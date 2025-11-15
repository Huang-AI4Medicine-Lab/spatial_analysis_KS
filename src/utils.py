#####

import os
# import openai

from datetime import datetime
from pathlib import Path
import pickle

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from gseapy import gsea, enrichr
import scipy.spatial as spatial
import seaborn as sns
import statsmodels
from matplotlib.collections import LineCollection
from scipy.stats import mannwhitneyu
from sklearn.cluster import DBSCAN, KMeans
from sklearn.neighbors import KernelDensity, KNeighborsClassifier, NearestNeighbors
from statsmodels.stats import multitest
from typing import Callable, List, Optional
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score
from datetime import datetime
from joblib import Parallel, delayed
from upsetplot import UpSet, from_memberships


plt.rcParams['svg.fonttype'] = 'none'


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, u):
        if self.parent[u] != u:
            self.parent[u] = self.find(self.parent[u])
        return self.parent[u]

    def union(self, u, v):
        root_u = self.find(u)
        root_v = self.find(v)
        if root_u != root_v:
            if self.rank[root_u] > self.rank[root_v]:
                self.parent[root_v] = root_u
            else:
                self.parent[root_u] = root_v
                if self.rank[root_u] == self.rank[root_v]:
                    self.rank[root_v] += 1

def plot_clusters(adata, title, column='core_id', ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 10))  # Change figsize to fit both x and y dimensions
    num_clusters = adata.obs[column].nunique()
    colormap = sns.color_palette("tab20", num_clusters)
    scatter = ax.scatter(adata.obs['x_centroid'], adata.obs['y_centroid'], c=adata.obs[column].astype(int), cmap=ListedColormap(colormap), s=10)
    ax.set_xlabel('x_centroid')
    ax.set_ylabel('y_centroid')
    ax.set_title(title)
    ax.invert_yaxis()  # Invert the y-axis
    if ax is None:  # Only add colorbar if ax was created internally
        fig.colorbar(scatter, ax=ax, label=column)

    # Annotate the plot with core IDs
    for core_id in np.unique(adata.obs[column]):
        if core_id != -1:  # Ignore noise points labeled as -1
            # Get the coordinates of the core
            core_coords = adata.obsm['spatial'][adata.obs[column] == core_id]  # Filter by core_id
            # Compute the centroid of the core for labeling
            centroid_x = np.mean(core_coords[:, 0])
            centroid_y = np.mean(core_coords[:, 1])
            ax.text(centroid_x, centroid_y, str(core_id), color='red', fontsize=16, ha='center', va='center')

def merge_close_clusters(adata, distance_threshold):
    # Calculate the centroids of each cluster
    core_centroids = adata.obs[adata.obs['core_id'] != -1].groupby('core_id')[['x_centroid', 'y_centroid']].mean()
    core_ids = core_centroids.index.to_list()
    n = len(core_centroids)

    # Compute pairwise distances between cluster centroids
    distances = pdist(core_centroids, metric='euclidean')
    distance_matrix = squareform(distances)
    
    # Create a DataFrame to store the distance matrix
    distance_df = pd.DataFrame(distance_matrix, index=core_ids, columns=core_ids)

    # Use Union-Find to merge clusters within the distance threshold
    uf = UnionFind(n)
    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i, j] < distance_threshold:
                uf.union(i, j)

    # Create the merge_map from the Union-Find structure
    merge_map = {core_ids[i]: core_ids[uf.find(i)] for i in range(n)}
    
    # Propagate the merge_map to ensure all clusters point to the final root
    final_merge_map = {}
    for core_id in merge_map:
        root = core_id
        while merge_map[root] != root:
            root = merge_map[root]
        final_merge_map[core_id] = root

    # Update core_id to merge close clusters
    adata.obs['core_id_merged'] = adata.obs['core_id'].map(lambda x: final_merge_map.get(x, x))

    # Ensure new labels starting from 0, preserve -1 for noise points
    unique_labels = np.sort(adata.obs['core_id_merged'].unique())
    label_map = {old_label: (new_label if old_label != -1 else -1) for new_label, old_label in enumerate(unique_labels)}
    adata.obs['core_id_merged'] = adata.obs['core_id_merged'].map(label_map).astype(int)

    return adata, distance_df

def combine_core_ids(adata, merge_tuples):
    # Create a mapping based on the merge_tuples
    merge_map = {}
    for target, source in merge_tuples:
        merge_map[source] = target
    
    # Update the core_id_merged column based on the merge_map
    def update_core_id(core_id):
        while core_id in merge_map:
            core_id = merge_map[core_id]
        return core_id
    
    adata.obs['core_id_merged'] = adata.obs['core_id_merged'].map(update_core_id)
    
    # Add core centroid coordinates to obs
    core_centroids_merged = adata.obs[adata.obs['core_id_merged'] != -1].groupby('core_id_merged')[['x_centroid', 'y_centroid']].mean()
    adata.obs['core_centroid_x'] = adata.obs['core_id_merged'].map(core_centroids_merged['x_centroid'])
    adata.obs['core_centroid_y'] = adata.obs['core_id_merged'].map(core_centroids_merged['y_centroid'])
    
    return adata

def load_xenium_TMA(sample, base_path, output_csv, save=True, return_adata=True, show=True, eps=250, min_samples=100, distance_threshold=500, min_cells_per_core=100):
    """
    Loads and processes Xenium TMA (Tissue Microarray) data for a given sample.

    Parameters:
    - sample (str): The sample identifier.
    - base_path (str): The base path where sample data is located.
    - output_csv (str): The path to save the resulting DataFrame as a CSV file.
    - save (bool): Whether to save the resulting DataFrame to a CSV file. Default is True.
    - return_adata (bool): Whether to return the processed AnnData object. Default is True.
    - show (bool): Whether to visualize the identified cores with core IDs. Default is True.
    - eps (float): The maximum distance between two samples for them to be considered as in the same neighborhood in DBSCAN. Default is 250.
    - min_samples (int): The number of samples (or total weight) in a neighborhood for a point to be considered as a core point in DBSCAN. Default is 100.
    - distance_threshold (float): The distance threshold for merging close clusters. Default is 500.
    - min_cells_per_core (int): The minimum number of cells required per core after merging clusters. Default is 100.

    Returns:
    - adata (AnnData): The processed AnnData object if return_adata is True.
    
    """
    # Read the 10x h5 file
    adata = sc.read_10x_h5(filename=f"{base_path}/{sample}/cell_feature_matrix.h5")

    # Read the cells CSV file
    df_ = pd.read_csv(f"{base_path}/{sample}/cells.csv.gz")

    # Update adata.obs_names and set the DataFrame index
    adata.obs_names = [f"{sample}_{r}" for r in range(len(adata))]
    df_.set_index(adata.obs_names, inplace=True)
    
    # Copy the DataFrame to adata.obs and add 'sample_id'
    adata.obs = df_.copy()
    adata.obs['sample_id'] = np.repeat(sample, len(adata))
    
    # Filter out cells with total_counts <= 0
    adata = adata[adata.obs.total_counts > 0]

    # Set the spatial coordinates in obsm
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
    
    print_with_time("Clustering cells")
    # Use DBSCAN to identify clusters (cores) based on spatial coordinates
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(adata.obsm["spatial"])
    adata.obs['core_id'] = db.labels_.astype(int)
    
    # Count cells per merged cluster
    cluster_counts = adata.obs['core_id'].value_counts()
    
    # Filter out clusters with fewer than min_cells_per_core cells
    clusters_to_keep = cluster_counts[cluster_counts >= 30].index
    adata = adata[adata.obs['core_id'].isin(clusters_to_keep)]
    
    # Order core_ids based on the centroid positions
    core_centroids = adata.obs[adata.obs['core_id'] != -1].groupby('core_id')[['x_centroid', 'y_centroid']].mean()
    core_centroids = core_centroids.reset_index()
    core_centroids_sorted = core_centroids.sort_values(by=['y_centroid', 'x_centroid']).reset_index(drop=True)
    core_centroids_sorted['new_core_id'] = range(len(core_centroids_sorted))

    core_id_map = dict(zip(core_centroids_sorted['core_id'], core_centroids_sorted['new_core_id']))
    core_id_map[-1] = -1  # Preserve -1 for noise points

    adata.obs['core_id'] = adata.obs['core_id'].replace(core_id_map).astype(int)
    
    print_with_time("Merging clusters")
    # Merge close clusters based on spatial proximity
    adata, distance_df = merge_close_clusters(adata, distance_threshold)
    
    # Count cells per merged cluster
    cluster_counts = adata.obs['core_id_merged'].value_counts()
    
    # Filter out clusters with fewer than min_cells_per_core cells
    clusters_to_keep = cluster_counts[cluster_counts >= min_cells_per_core].index
    adata = adata[adata.obs['core_id_merged'].isin(clusters_to_keep)]
    
    distance_df.to_csv(f'data/{sample}_distance_df_clusters.csv')
    adata.obs['TMA_core_id'] = f"TMA_{sample.split('_')[2]}_" + adata.obs['core_id_merged'].astype(str)
    adata.obs_names = [f"{sample}_{r}" for r in range(len(adata))]

    adata = adata[adata.obs.query("core_id_merged != -1").index]
    
    # Add core centroid coordinates to obs
    core_centroids_merged = adata.obs[adata.obs['core_id_merged'] != -1].groupby('core_id_merged')[['x_centroid', 'y_centroid']].mean()
    adata.obs['core_centroid_x'] = adata.obs['core_id_merged'].map(core_centroids_merged['x_centroid'])
    adata.obs['core_centroid_y'] = adata.obs['core_id_merged'].map(core_centroids_merged['y_centroid'])
    
    # Create a new DataFrame from adata.X with var_names as columns
    df = pd.DataFrame(adata.X.toarray(), columns=adata.var_names)

    # Set the index of the new DataFrame to match adata.obs.cell_id
    df.index = adata.obs.cell_id

    # Merge the new DataFrame with adata.obs
    df = pd.merge(df, adata.obs.copy(), left_index=True, right_index=True)
    
    if save:
        # Save the resulting DataFrame to a CSV file
        df.to_csv(output_csv)
    
    if show:
        print_with_time("Visualizing the identified cores with core IDs")
        fig, axs = plt.subplots(1, 2,  figsize=(9, 9))  # Adjust figsize as needed
        plot_clusters(adata, 'Initial Clusters', column='core_id', ax=axs[0])
        plot_clusters(adata, 'Refined Clusters', column='core_id_merged', ax=axs[1])
        
        plt.suptitle(f"{sample}; eps={eps}; min_samples={min_samples}")
        plt.tight_layout()
        
        plt.savefig(f"figures/cores/{sample}_before_after_refinement.png")
        plt.show()
    
    if return_adata:
        return adata

# Function to rotate the grid for a given TMA slide counterclockwise
def rotate_tma_grid_ccw(df, start_row):
    grid = df.iloc[start_row:start_row+3, 1:7].values
    rotated_grid = grid.T
    return rotated_grid


# Function to map unique IDs from tma_grids to the centroid locations in adata
def map_tma_grids_to_centroids(adata, tma_grids, y_threshold=1000, x_slack=100):
    for sample_id in adata.obs['sample_id'].unique():
        one_adata = adata[adata.obs['sample_id'] == sample_id]
        tma_slide_id = ('_').join(sample_id.split('_')[1:3])  # Extract TMA slide ID correctly

        if tma_slide_id in tma_grids:
            grid = tma_grids[tma_slide_id]
            points_set = set(zip(one_adata.obs['core_centroid_x'], one_adata.obs['core_centroid_y']))
            points = list(points_set)  # Convert set to list
            
            # Convert points to DataFrame
            points_df = pd.DataFrame(points, columns=['x', 'y'])
            
            # Assign rows based on y-coordinates (similar y points to the same row)
            points_df = points_df.sort_values(by='y', ascending=False).reset_index(drop=True)
            # Define a threshold for grouping similar y-coordinates using the y_threshold and identify rows
            points_df['row'] = (points_df['y'].diff().abs() > y_threshold).cumsum()
            
            # Ensure we only have 6 rows by reassigning rows if there are more than 6
            if points_df['row'].nunique() > 6:
                unique_rows = points_df['row'].unique()[:6]
                points_df = points_df[points_df['row'].isin(unique_rows)]

            points_df['row'] = points_df.groupby('row').ngroup()

            # Calculate the global minimum and maximum x-coordinates
            min_x = one_adata.obs['core_centroid_x'].min()
            max_x = one_adata.obs['core_centroid_x'].max()

            # Within each row, assign columns based on x-coordinates with slack
            result_df = pd.DataFrame(columns=['x', 'y', 'row', 'col'])
            for row_id, group in points_df.groupby('row'):
                group = group.sort_values(by='x').reset_index(drop=True)
                n_points = len(group)
                cols = [-1] * n_points  # Initialize with invalid column index

                if n_points == 1:
                    # Decide based on x-coordinate if the single point is more towards left, middle, or right
                    if group.iloc[0]['x'] < (min_x + max_x) / 3:
                        cols[0] = 0  # LEFT
                    elif group.iloc[0]['x'] > 2 * (min_x + max_x) / 3:
                        cols[0] = 2  # RIGHT
                    else:
                        cols[0] = 1  # MIDDLE
                elif n_points == 2:
                    if group.iloc[0]['x'] < (min_x + max_x) / 3:
                        cols[0] = 0  # LEFT
                        if group.iloc[1]['x'] > 2 * (min_x + max_x) / 3:
                            cols[1] = 2  # RIGHT
                        else:
                            cols[1] = 1  # MIDDLE
                    elif group.iloc[0]['x'] > 2 * (min_x + max_x) / 3:
                        cols[0] = 2  # RIGHT
                        if group.iloc[1]['x'] < (min_x + max_x) / 3:
                            cols[1] = 0  # LEFT
                        else:
                            cols[1] = 1  # MIDDLE
                    else:
                        cols[0] = 1  # MIDDLE
                        if group.iloc[1]['x'] < (min_x + max_x) / 3:
                            cols[1] = 0  # LEFT
                        else:
                            cols[1] = 2  # RIGHT
                elif n_points == 3:
                    cols = [0, 1, 2]  # LEFT, MIDDLE, RIGHT

                group['col'] = cols
                result_df = pd.concat([result_df, group], ignore_index=True)

            # Ensure we only have 3 columns
            result_df = result_df[result_df['col'] < 3]

            # Map each sorted point to the corresponding unique ID
            for _, row in result_df.iterrows():
                if row['row'] < 6 and row['col'] < 3:
                    unique_id = grid[int(row['row'])][int(row['col'])]
                    adata.obs.loc[(adata.obs['core_centroid_x'] == row['x']) & (adata.obs['core_centroid_y'] == row['y']), 'unique_id'] = unique_id
    
    return adata

# Plot the centroid locations with unique_id labels, with inverted y-axis
def plot_centroid_locations(adata, save=True):
    for sample in adata.obs['sample_id'].unique():
        one_adata = adata[adata.obs['sample_id'] == sample]
        points = set(zip(one_adata.obs['core_centroid_x'], one_adata.obs['core_centroid_y']))
        
        # Unpack the points into x and y coordinates
        x, y = zip(*points)

        # Create a scatter plot
        plt.figure(figsize=(5, 12))
        plt.scatter(x, y, color='blue')
        
        # Add labels
        for i, point in enumerate(points):
            unique_id = one_adata.obs.loc[(one_adata.obs['core_centroid_x'] == point[0]) & (one_adata.obs['core_centroid_y'] == point[1]), 'unique_id'].values[0]
            plt.text(point[0], point[1]-150, unique_id, fontsize=9, ha='center')

        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title(f'Scatter Plot of (x, y) Points for Sample {sample}')
        plt.gca().invert_yaxis()  # Invert the y-axis
        plt.grid(True)
        
        if save:
            plt.savefig(f"figures/{sample}_TMA_centroid_locations.png")
            
        plt.show()

def print_with_time(message):
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] {message}")
    
def calculate_auc(cell_type, cluster, clustername='leiden', verbose=False):
    try:
        meta_gene = marker_list[cell_type]
        meta_gene_expression = adata[:, meta_gene].X.sum(axis=1)
        labels = (adata.obs[clustername] == cluster).astype(int)
        
        if verbose:
            # Print or log intermediate values for debugging
            print(f"Calculating AUC for {cell_type}, {cluster}")
            print(f"Meta gene markers: {meta_gene}")
            print(f"Meta gene expression shape: {meta_gene_expression.shape}")
            print(f"Labels shape: {labels.shape}")
        
        # Calculate AUC for the cell type and cluster
        auc = roc_auc_score(labels, np.array(meta_gene_expression))
        
        # print(f"AUC calculated successfully: {auc}")
    except Exception as e:
        print_with_time(f"Error calculating AUC for {cell_type}, {cluster}: {e}")
        auc = np.nan
        
    return auc


def perform_auc_and_assignment(adata, marker_list, subclusters=None, auc_threshold=None, leiden_id_to_subcluster=None):
    cell_types = list(marker_list.keys())
    clusters = subclusters if subclusters is not None else adata.obs['leiden'].cat.categories
    
    if subclusters is not None and leiden_id_to_subcluster is None:
        raise ValueError("If subclusters are provided, leiden_id_to_subcluster parameter is required.")
    
    # Determine clustername based on conditions
    if subclusters is None:
        clustername = 'leiden'
    else:
        clustername = f"leiden_sub_{leiden_id_to_subcluster}"
    
    roauc_df = pd.DataFrame(0.0, index=cell_types, columns=clusters)

    # Parallel processing to calculate AUC values
    results = Parallel(n_jobs=-1)(delayed(calculate_auc)(cell_type, cluster, clustername) for cell_type in cell_types for cluster in clusters)

    # Reshape the results into the DataFrame
    result_idx = 0
    for cell_type in cell_types:
        for cluster in clusters:
            roauc_df.at[cell_type, cluster] = results[result_idx]
            result_idx += 1

    # Create a mapping from clusters to cell types based on maximum AUC
    cluster_to_cell_type = {}
    second_best_cell_type = {}
    third_best_cell_type = {}
    
    if not auc_threshold:
        auc_threshold = np.percentile(roauc_df.values.flatten(), 70)
        
    for cluster in clusters:
        sorted_auc = roauc_df[cluster].sort_values(ascending=False)
        best_cell_type = sorted_auc.index[0] if sorted_auc.iloc[0] >= auc_threshold else "unknown"
        cluster_to_cell_type[cluster] = best_cell_type
        second_best_cell_type[cluster] = sorted_auc.index[1] if len(sorted_auc) > 1 else "unknown"
        third_best_cell_type[cluster] = sorted_auc.index[2] if len(sorted_auc) > 2 else "unknown"
    
    
    if subclusters is None:
        # Map the cell types to the clusters in adata.obs
        adata.obs['cell_type'] = adata.obs[clustername].map(cluster_to_cell_type)

        # Add second and third best cell types to adata.obs
        adata.obs['second_best_cell_type'] = adata.obs[clustername].map(second_best_cell_type)
        adata.obs['third_best_cell_type'] = adata.obs[clustername].map(third_best_cell_type)
    
    else:
        # Ensure the categorical dtype for cell_type and add new categories
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
        adata.obs['second_best_cell_type'] = adata.obs['second_best_cell_type'].astype('category')
        adata.obs['third_best_cell_type'] = adata.obs['third_best_cell_type'].astype('category')

        # Map the cell types to the clusters in adata.obs
        for cluster in clusters:
            cell_type = cluster_to_cell_type[cluster]
            second_best = second_best_cell_type[cluster]
            third_best = third_best_cell_type[cluster]

            try:
                # Add new categories if they don't exist
                if second_best not in adata.obs['cell_type'].cat.categories:
                    adata.obs['cell_type'].cat.add_categories([second_best], inplace=True)
                if third_best not in adata.obs['cell_type'].cat.categories:
                    adata.obs['cell_type'].cat.add_categories([third_best], inplace=True)
                if second_best not in adata.obs['second_best_cell_type'].cat.categories:
                    adata.obs['second_best_cell_type'].cat.add_categories([second_best], inplace=True)
                if third_best not in adata.obs['third_best_cell_type'].cat.categories:
                    adata.obs['third_best_cell_type'].cat.add_categories([third_best], inplace=True)
            except:
                pass

            # Assign cell types
            adata.obs.loc[adata.obs[f'leiden_sub_{leiden_id_to_subcluster}'] == cluster, 'cell_type'] = cell_type
            adata.obs.loc[adata.obs[f'leiden_sub_{leiden_id_to_subcluster}'] == cluster, 'second_best_cell_type'] = second_best
            adata.obs.loc[adata.obs[f'leiden_sub_{leiden_id_to_subcluster}'] == cluster, 'third_best_cell_type'] = third_best

    return adata, roauc_df, cluster_to_cell_type

def plot_results(adata, roauc_df, suffix=''):
    # Plot clustermap of AUC scores
    sns.clustermap(roauc_df, cmap="vlag", standard_scale=1, figsize=(10,10), fmt=".2f")
    plt.show()

    # Plot heatmap of AUC scores
    plt.figure(figsize=(20, 10))
    sns.heatmap(roauc_df, annot=True, cmap='viridis', fmt=".2f")
    plt.title("Heatmap of Per-Cluster AUC Scores for All Available Cell Types")
    plt.xlabel("Clusters")
    plt.ylabel("Cell Types")
    
    plt.savefig(f"figures/AUC_heatmap_{suffix}.png")
            
    plt.show()

    # Plot the final UMAP
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": (300)}):
        sc.pl.umap(adata, color=['leiden', 'cell_type'], legend_loc='on data', size=0.5, save=f"_celltyped_{suffix}.png")

def combine_clusters(adata, cluster_to_cell_type):
    # Step 1: Create a dictionary where the keys are cell types and the values are lists of cluster IDs sharing that cell type
    clusters_by_cell_type = {}
    for cluster, cell_type in cluster_to_cell_type.items():
        if cell_type not in clusters_by_cell_type:
            clusters_by_cell_type[cell_type] = []
        clusters_by_cell_type[cell_type].append(cluster)
    
    # Step 2: For each cell type, find the minimum cluster ID by converting the cluster names to integers and then back to strings
    min_cluster_id_by_cell_type = {}
    for cell_type, clusters in clusters_by_cell_type.items():
        min_cluster_id = str(min(int(cluster) for cluster in clusters))
        min_cluster_id_by_cell_type[cell_type] = min_cluster_id
    
    # Step 3: Replace all cluster IDs sharing the same cell type with the minimum cluster ID
    adata_copy = adata.copy()
    adata_copy.obs['leiden'] = adata_copy.obs['leiden'].astype(str)
    for cell_type, min_cluster_id in min_cluster_id_by_cell_type.items():
        for cluster in clusters_by_cell_type[cell_type]:
            adata_copy.obs.loc[adata_copy.obs['leiden'] == cluster, 'leiden'] = min_cluster_id
    
    # Ensure the leiden_0.7 column is categorical with updated categories
    unique_clusters = sorted(adata_copy.obs['leiden'].unique(), key=int)
    adata_copy.obs['leiden'] = pd.Categorical(adata_copy.obs['leiden'], categories=unique_clusters)
    
    # Step 4: Add assert statements to check for correctness
    for cell_type, clusters in clusters_by_cell_type.items():
        min_cluster_id = min_cluster_id_by_cell_type[cell_type]
        for cluster in clusters:
            assert all(adata_copy.obs.loc[adata_copy.obs['leiden'] == cluster, 'leiden'] == min_cluster_id)
    
    assert adata_copy.obs['leiden'].nunique() < adata.obs['leiden'].nunique(), "The number of unique clusters should be less than the original"
    
    # Step 5: Create an updated cluster_to_cell_type mapping
    updated_cluster_to_cell_type = {}
    for cell_type, min_cluster_id in min_cluster_id_by_cell_type.items():
        updated_cluster_to_cell_type[min_cluster_id] = cell_type
    
    return adata_copy, updated_cluster_to_cell_type

def combine_clusters_and_update_mapping(adata, cluster_to_cell_type, cluster_name='leiden'):
    # Step 1: Create a dictionary where the keys are cell types and the values are lists of cluster IDs sharing that cell type
    clusters_by_cell_type = {}
    for cluster, cell_type in cluster_to_cell_type.items():
        if cell_type not in clusters_by_cell_type:
            clusters_by_cell_type[cell_type] = []
        clusters_by_cell_type[cell_type].append(cluster)
    
    # Step 2: For each cell type, find the minimum cluster ID by converting the cluster names to integers and then back to strings
    min_cluster_id_by_cell_type = {}
    for cell_type, clusters in clusters_by_cell_type.items():
        min_cluster_id = str(min(int(cluster) for cluster in clusters))
        min_cluster_id_by_cell_type[cell_type] = min_cluster_id
    
    # Step 3: Replace all cluster IDs sharing the same cell type with the minimum cluster ID
    adata_copy = adata.copy()
    adata_copy.obs[cluster_name] = adata_copy.obs[cluster_name].astype(str)
    for cell_type, min_cluster_id in min_cluster_id_by_cell_type.items():
        for cluster in clusters_by_cell_type[cell_type]:
            adata_copy.obs.loc[adata_copy.obs[cluster_name] == cluster, cluster_name] = min_cluster_id
    
    # Ensure the leiden_0.7 column is categorical with updated categories
    unique_clusters = sorted(adata_copy.obs[cluster_name].unique(), key=int)
    adata_copy.obs[cluster_name] = pd.Categorical(adata_copy.obs[cluster_name], categories=unique_clusters)
    
    # Step 4: Add assert statements to check for correctness
    for cell_type, clusters in clusters_by_cell_type.items():
        min_cluster_id = min_cluster_id_by_cell_type[cell_type]
        for cluster in clusters:
            assert all(adata_copy.obs.loc[adata_copy.obs[cluster_name] == cluster, cluster_name] == min_cluster_id)
    
    assert adata_copy.obs[cluster_name].nunique() < adata.obs[cluster_name].nunique(), "The number of unique clusters should be less than the original"
    
    # Step 5: Create an updated cluster_to_cell_type mapping
    updated_cluster_to_cell_type = {}
    for cell_type, min_cluster_id in min_cluster_id_by_cell_type.items():
        updated_cluster_to_cell_type[min_cluster_id] = cell_type
    
    return adata_copy, updated_cluster_to_cell_type


def subcluster(adata, leiden_id, resolution=1.0):

    """
    Subset an AnnData object based on the leiden column ID and perform leiden clustering on the subset cells
    with a specific resolution. Adds the new clusters back to the original AnnData object with a modified ID.
    
    Parameters:
    - adata: AnnData object
    - leiden_id: ID of the leiden cluster to subset
    - resolution: Resolution parameter for the leiden clustering
    
    Returns:
    - adata: Updated AnnData object with new subclusters added
    - new_subcluster_ids: List of new subcluster IDs
    """
    # Subset the AnnData object
    subset_adata = adata[adata.obs['leiden'] == str(leiden_id)].copy()
    
    # Perform leiden clustering on the subset cells
    sc.pp.neighbors(subset_adata)
    sc.tl.leiden(subset_adata, resolution=resolution, key_added='leiden_sub')
    
    # Create new subcluster IDs with _A, _B, etc.
    subcluster_mapping = {str(old): f"{leiden_id}_{chr(65 + i)}" for i, old in enumerate(subset_adata.obs['leiden_sub'].cat.categories)}
    subset_adata.obs['leiden_sub'] = subset_adata.obs['leiden_sub'].map(subcluster_mapping)
    
    # Add new subclusters back to the original AnnData object
    adata.obs[f'leiden_sub_{leiden_id}'] = adata.obs['leiden'].astype(str)
    adata.obs.loc[subset_adata.obs_names, f'leiden_sub_{leiden_id}'] = subset_adata.obs['leiden_sub']
    
    # Return the new subcluster IDs
    new_subcluster_ids = list(subcluster_mapping.values())
    
    return adata, new_subcluster_ids


# Function to clean cell type names
def clean_cell_type(cell_type):
    if " / " in cell_type:
        return cell_type.split(" / ")[0]
    return cell_type

# Function to calculate centroids and create plots for a given subset
def plot_sample_data(adata_subset, sample_id):
    # Calculate the centroid of spatial coordinates for each TMA core
    centroids = adata_subset.obs.groupby('sample_TMA_core_id')[['x_centroid', 'y_centroid']].mean().reset_index()

    # Normalize the centroid coordinates to fit into a 6x3 grid
    x_max, y_max = centroids['x_centroid'].max(), centroids['y_centroid'].max()
    centroids['row'] = np.floor(centroids['y_centroid'] / y_max * 6).astype(int)
    centroids['col'] = np.floor(centroids['x_centroid'] / x_max * 3).astype(int)

    # Ensure the indices are within bounds
    centroids['row'] = centroids['row'].clip(0, 5)
    centroids['col'] = centroids['col'].clip(0, 2)

    # Calculate cell type proportions for each TMA core
    cell_type_counts = adata_subset.obs.groupby('sample_TMA_core_id')['cell_type'].value_counts(normalize=True).unstack(fill_value=0)

    # Create a 6x3 subplot layout
    fig, axes = plt.subplots(6, 3, figsize=(16, 24), sharex=True, sharey=True)
    fig.tight_layout(pad=5.0)

    # Plot each TMA core in its corresponding subplot
    for idx, row in centroids.iterrows():
        r, c = row['row'], row['col']
        tma_core_id = row['sample_TMA_core_id']
        
        if tma_core_id in cell_type_counts.index:
            cell_type_prop = cell_type_counts.loc[tma_core_id]
            # Clean cell type names and get corresponding colors
            cell_types_cleaned = cell_type_prop.index.map(clean_cell_type)
            colors = [color_palette_subtype.get(cell_type, "#000000") for cell_type in cell_types_cleaned]
            cell_type_prop.plot(kind='bar', ax=axes[r, c], color=colors, width=1.3)
            axes[r, c].set_title(f"Core: {tma_core_id}", fontsize=10)
            axes[r, c].set_xlabel("Cell Type", fontsize=20)
            axes[r, c].set_ylabel("Proportion", fontsize=20)
            axes[r, c].tick_params(axis='x', labelsize=20)
            axes[r, c].tick_params(axis='y', labelsize=20)
        else:
            axes[r, c].set_visible(False)  # Hide axes without data

    # Adjust layout to make sure titles and labels are readable
    plt.subplots_adjust(hspace=0.2, wspace=0.2)
    plt.suptitle(f"Sample ID: {sample_id}", fontsize=20)
    plt.tight_layout()
    plt.savefig(f"figures/spatial_{sample_id}_celltype_percentage.png")
    plt.show()
    
def split_core_id(adata, core_id_to_split, eps=250, min_samples=100, sample=None):
    """
    Splits a specified core ID based on new DBSCAN clustering and reassigns all TMA_core_ids
    of the sample based on the x, y coordinate location starting from left to right.

    Parameters:
    - adata (AnnData): The AnnData object containing the data.
    - core_id_to_split (int): The core ID to split.
    - eps (float): The maximum distance between two samples for them to be considered as in the same neighborhood in DBSCAN. Default is 250.
    - min_samples (int): The number of samples (or total weight) in a neighborhood for a point to be considered as a core point in DBSCAN. Default is 100.
    - sample (str): The sample identifier. Required to update `TMA_core_id` in `adata.obs`.

    Returns:
    - adata (AnnData): The updated AnnData object.
    """
    # Ensure sample is provided
    if sample is None:
        raise ValueError("Sample identifier must be provided.")

    # Filter cells belonging to the core_id_to_split
    cells_to_split = adata.obs.query(f"core_id_merged == {core_id_to_split}").index

    # Perform DBSCAN clustering on these cells
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(adata[cells_to_split].obsm["spatial"])
    new_core_ids = db.labels_

    # Update core_id_merged for the cells to split
    adata.obs.loc[cells_to_split, 'core_id_merged'] = new_core_ids + adata.obs['core_id_merged'].max() + 1

    # Recalculate TMA_core_ids for the entire sample
    sample_cells = adata.obs.query(f"sample_id == '{sample}'").index
    core_centroids = adata.obs.loc[sample_cells].groupby('core_id_merged')[['x_centroid', 'y_centroid']].mean()
    core_centroids = core_centroids.reset_index()
    core_centroids_sorted = core_centroids.sort_values(by=['y_centroid', 'x_centroid']).reset_index(drop=True)
    core_centroids_sorted['new_core_id'] = range(len(core_centroids_sorted))

    core_id_map = dict(zip(core_centroids_sorted['core_id_merged'], core_centroids_sorted['new_core_id']))

    adata.obs.loc[sample_cells, 'core_id_merged'] = adata.obs.loc[sample_cells, 'core_id_merged'].map(core_id_map).astype(int)

    # Update TMA_core_id
    adata.obs['TMA_core_id'] = f"TMA_{sample.split('_')[2]}_" + adata.obs['core_id_merged'].astype(str)

    return adata




def calculate_neighborhood_cell_composition(df, n_neighbors=200, x='CenterX_global_px_zero', y='CenterY_global_px_zero', ctype_col='cell_type'):
    """
    This function calculates the composition of cell types in the neighborhoods 
    of each cell in a given DataFrame.
    
    Parameters:
    - df: DataFrame containing cell data, including x and y coordinates and cell types.
    - n_neighbors: Number of nearest neighbors to consider for each cell. Default is 200.
    
    Returns:
    - df: DataFrame with added columns representing the composition of cell types 
          in the neighborhoods of each cell.
    """

    # Extracting coordinates and cell types from the DataFrame
    coords = df[[x, y]].values
    cell_types = df[ctype_col].values
    
    # Obtaining unique cell types and sorting them
    unique_cell_types = sorted(df[ctype_col].unique())

    # Initializing a NearestNeighbors object and fitting it to the data
    neigh = NearestNeighbors(n_neighbors=n_neighbors)
    neigh.fit(coords)
    
    # Finding the indices of nearest neighbors for each point
    _, neighbors_indices = neigh.kneighbors(coords)

    # Mapping cell types to indices for faster processing
    cell_type_to_index = {cell_type: i for i, cell_type in enumerate(unique_cell_types)}
    cell_type_indices = np.vectorize(cell_type_to_index.get)(cell_types)

    # Initializing an array to hold the counts of cell types in neighborhoods
    cell_composition_counts = np.zeros((len(df), len(unique_cell_types)))

    # Printing progress information
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Calculating Neighborhood Cell Composition")

    # Counting the occurrences of each cell type in the neighborhoods
    for i, neighbors in enumerate(neighbors_indices):
        # Updating progress every 30000 iterations
        if i % 30000 == 0:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Processing cell ID: {i}")
        neighbor_types = cell_type_indices[neighbors]
        for neighbor_type in neighbor_types:
            cell_composition_counts[i, neighbor_type] += 1

    # Creating a DataFrame from the counts
    cell_composition_df = pd.DataFrame(cell_composition_counts, columns=[f'n_{ct}' for ct in unique_cell_types], index=df.index)

    # Joining the new DataFrame with the original one
    df = pd.merge(df, cell_composition_df, left_index=True, right_index=True)

    # Printing completion message
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Completed Neighborhood Cell Composition")

    return df



def remove_unnamed_columns(df):
    unnamed_columns = [col for col in df.columns if col.startswith('Unnamed:')]
    df = df.drop(columns=unnamed_columns, axis=1)
    return df


def mpl_default():
    """
    Some function to restore default values
    """
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.style.use('default')

def scatter_df(df, x='x', y='y', c='k', s=[1], figsize=(14,14), cmap=None, 
               dpi = 140, show_axis=True, ax=None, title=None, marker='o', alpha=1.0, label="",
              flip_y=False, flip_x=False, hue=None, vmin=None, vmax=None):
    
    dfc = df.copy()
        
    with mpl.rc_context({'figure.dpi':dpi, 'savefig.dpi':dpi, "figure.figsize":figsize}):
        
        if ax is None:
            f, ax = plt.subplots(figsize=figsize)
            
        # plt.rcParams['figure.dpi'] = dpi
        # plt.rcParams['savefig.dpi'] = dpi
        # plt.rcParams["figure.figsize"] = figsize
        
        if flip_y:
            dfc[y] = -dfc[y]

        if flip_x:
            dfc[x] = -dfc[x]
            
        if hue:
            ax.scatter(dfc[x].to_numpy(), dfc[y].to_numpy(), c=dfc[hue], s=s, cmap=cmap, marker=marker, alpha=alpha, label=label, vmin=None, vmax=None)
        else:
            ax.scatter(dfc[x].to_numpy(), dfc[y].to_numpy(), c=c, s=s, cmap=cmap, marker=marker, alpha=alpha, label=label, vmin=None, vmax=None)
        
        if not show_axis:
            plt.axis('off')
        if title:
            plt.title(title)

        plt.tight_layout()
    return ax
    
    
    
def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)
    

def spatial_neighborhood_density(cell_locations_df, x='x', y='y', gene=None, radius = 200.0, leafsize=10000, 
                                 figsize=(14,14), fontsize=9, cmap='Reds', dpi=100, show_axis=True, 
                                 return_results=True, kde=False, kde_levels=10, nlargest= 500, min_samples=10, cluster_color='red',s=1, marker='o'):
    """
    This function calculates the density of individual cells in the spatial dataset and optionally returns the density vector [:return(neighbors, frequency, ax)].
    """
    tree = spatial.KDTree(cell_locations_df[[x,y]].to_numpy(), leafsize=5000)
    neighbors = tree.query_ball_tree(tree, radius)
    frequency = np.array([len(i) for i in neighbors])
    # density = frequency/radius**2
    plt.rcParams["font.size"]=fontsize
    
    if gene:
        gene_expression = cell_locations_df[gene]
        ax = scatter_df(cell_locations_df,  x=x, y=y, figsize=figsize, c=np.multiply(frequency, gene_expression), cmap=cmap, dpi=dpi, s=s, marker=marker,show_axis=show_axis)
        if return_results:
            return(neighbors, frequency, gene_expression)
    else:
        ax = scatter_df(cell_locations_df, x=x, y=y, figsize=figsize, c=frequency, cmap=cmap, dpi=dpi, marker=marker, show_axis=show_axis, s=s)
        
        if kde:
            kdeplot_out = sns.kdeplot(data=cell_locations_df, x=x, y=y, ax=ax, levels=kde_levels, shade=True, cmap='viridis', alpha=0.35, c=frequency)
            kdedata = []
            for i in kdeplot_out.get_children():
                if i.__class__.__name__ == 'LineCollection':
                    kdedata.append(i.get_paths())
        
        if return_results:
            return(neighbors, frequency, ax)
        
        
def spatial_KDE_plot(cell_locations_df, x='x', y='y', scatter=False, figsize=(14,14), fontsize=9, cmap='viridis', alpha_dots=0.5, dpi=100, show_axis=True, title=None,
                                 return_results=True, kde_levels=10, nlargest= 500, min_samples=10, cluster_color='red',s=1, marker='o', alpha=0.5):
    """
    This function calculates the density of individual cells in the spatial dataset and optionally returns the density vector.
    """
    
    if scatter:
        ax = scatter_df(cell_locations_df, figsize=figsize, c='k', dpi=dpi, marker=marker, show_axis=show_axis, s=0.3, alpha=alpha_dots)
        kdeplot_out = sns.kdeplot(data=cell_locations_df, x=x, y=y, ax=ax, levels=kde_levels, shade=True, cmap=cmap, alpha=alpha)
    else:
        kdeplot_out = sns.kdeplot(data=cell_locations_df, x=x, y=y, levels=kde_levels, shade=True, cmap=cmap, alpha=alpha)
    
    kdedata = []
    for i in kdeplot_out.get_children():
        if i.__class__.__name__ == 'LineCollection':
            kdedata.append(i.get_paths())

    if not show_axis:
        plt.axis('off')
    if title:
        plt.title(title)
        
    if return_results:
        return(kdedata)


# openai.api_key = os.getenv("OPENAI_API_KEY")
# def ask_gpt(query=None,  model="gpt-3.5-turbo", system_role="You are an expert in coding", temperature=1, max_tokens=256, top_p=1, 
#             frequency_penalty=0, presence_penalty=0, return_full=False):
    
#     response = openai.ChatCompletion.create(
#       model=model,
#       messages=[
#         {
#           "role": "system",
#           "content": system_role
#         },
#         {
#           "role": "user",
#           "content": query
#         }
#       ],
#       temperature=temperature,
#       max_tokens=max_tokens,
#       top_p=top_p,
#       frequency_penalty=frequency_penalty,
#       presence_penalty=presence_penalty
#     )

#     print(response.choices[0].message.content)

#     if return_full:
#         return(response)
    

def run_DE_high_low(explanation_df_norm, annotation_label_dict, region, gene_list,
                          percentage_high, percentage_low):
    explanation_df_norm_q = explanation_df_norm.query(f"annotation == {annotation_label_dict[region]} and cell_type == 17")

    column_masks_high = {}
    column_masks_low = {}    

    explanation_df_sum = explanation_df_norm_q[attr_gene_list].sum(axis=1)
    
    # Define the maximum value and percentages for high and low thresholds
    max_value = explanation_df_sum.max()

    # Calculate threshold values based on percentages
    threshold_value_high = max_value * percentage_high
    threshold_value_low = max_value * percentage_low

    row_mask_high = explanation_df_sum >= threshold_value_high
    row_mask_low = explanation_df_sum < threshold_value_low

    explanation_df_norm_high = explanation_df_norm_q[row_mask_high]
    explanation_df_norm_low = explanation_df_norm_q[row_mask_low]
    
    print(len(explanation_df_norm_high), len(explanation_df_norm_low))
    above_threshold_df = explanation_df_norm_high[gene_list]
    below_threshold_df = explanation_df_norm_low[gene_list]
    
    adata_above_threshold = sc.AnnData(X=above_threshold_df.values)
    adata_below_threshold = sc.AnnData(X=below_threshold_df.values)

    adata_above_threshold.var_names = above_threshold_df.columns
    adata_below_threshold.var_names = below_threshold_df.columns

    adata_above_threshold.obs['dataset'] = 'above'
    adata_below_threshold.obs['dataset'] = 'below'

    combined_adata = adata_above_threshold.concatenate(adata_below_threshold)

    sc.tl.rank_genes_groups(combined_adata, groupby='dataset', method='wilcoxon')
    result_df = sc.get.rank_genes_groups_df(combined_adata, group='above')
    
    filtered_result_df = result_df[result_df['pvals_adj'] < 0.01]
    filtered_result_df.index = filtered_result_df.names
    
    top_upregulated = filtered_result_df.nlargest(20, 'logfoldchanges')
    top_downregulated = filtered_result_df.nsmallest(20, 'logfoldchanges')

    top_genes = pd.concat([top_upregulated, top_downregulated])
    top_genes.sort_values(by='logfoldchanges', inplace=True, ascending=False)
    top_genes = top_genes[~top_genes.index.duplicated(keep='first')]
    
    top_genes['pvals_adj'] = top_genes['pvals_adj'] + 1e-235
    names = top_genes.index
    logfoldchanges = top_genes['logfoldchanges']
    pvals_adj = -np.log10(top_genes['pvals_adj'])

    fig, axes = plt.subplots(1, 1, figsize=(20, 2))
    sns.heatmap(
        top_genes[['logfoldchanges']].T,
        cmap="RdBu_r",
        annot=False,
        fmt=".2f",
        center=0,
        cbar=True,
        linewidths=0.5,
        annot_kws={"size": 10},
        cbar_kws={"label": "logfoldchanges"},
    )
    plt.tight_layout()
    plt.title(f'Top DE genes (high vs. low attr) of {region}')
    plt.show()

    """# top_genes_sorted = top_genes[top_genes['logfoldchanges'] > 0].sort_values(by='logfoldchanges', ascending=False)
    top_genes_sorted = top_genes.sort_values(by='logfoldchanges', ascending=False)
    gene_list_ = top_genes_sorted.index.tolist()

    enr = gp.enrichr(gene_list=gene_list_,
                     gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human'],
                     organism='human',
                     outdir='gsea_results',
                    )

    # enr.res2d.to_csv(f'./results/GCN-subgraph-GSEA_MSigDB_Hallmark_2020_KEGG_2021_Human_top_genes_{region}_high_vs_low.csv')
    gsea_results = enr.res2d

    gsea_results_sorted = gsea_results.sort_values(by='Combined Score', ascending=False)
    top_results = gsea_results_sorted.head(15)

    plt.figure(figsize=(2, 4), dpi=200)
    plt.barh(top_results['Term'], top_results['Combined Score'])
    plt.xlabel('Normalized Enrichment Score')
    plt.ylabel('Gene Set Term')
    plt.title(f'Top Gene Set Enrichment Results for {region}')
    plt.gca().invert_yaxis()
    plt.show()"""
    
    ranked_result_df = filtered_result_df['logfoldchanges']
    ranked_result_df.index = filtered_result_df['names']
    ranked_result_df = pd.DataFrame(ranked_result_df)

    # # run prerank
    # # enrichr libraries are supported by prerank module. Just provide the name
    # # use 4 process to acceralate the permutation speed
    pre_res = gp.prerank(rnk=ranked_result_df,
                         gene_sets=['MSigDB_Hallmark_2021','KEGG_2021_Human'],
                         threads=4,
                         min_size=5,
                         max_size=1000,
                         permutation_num=1000, # reduce number to speed up testing
                         outdir=None, # don't write to disk
                         seed=6,
                         verbose=True, # see what's going on behind the scenes
                        )
    
    gsea_results_sorted = pre_res.res2d.sort_values(by='NES')
    gsea_results_sorted = gsea_results_sorted[gsea_results_sorted['NOM p-val'] < 0.01]

    # Bar plot for NES and Terms
    plt.figure(figsize=(8, 6))
    plt.barh(gsea_results_sorted['Term'], gsea_results_sorted['NES'], color='skyblue')
    plt.xlabel('NES (Normalized Enrichment Score)')
    plt.title(f'Enrichment Scores for Gene Sets in {region} high vs. low')
    plt.grid(axis='x')  # Add gridlines along the x-axis
    plt.tight_layout()
    plt.show()
    

def run_DE_A_vs_B(explanation_df_norm_A, explanation_df_norm_B, region_A, region_B, annotation_label_dict, gene_list,
                          percentage_high, percentage_low, pval_de=0.01, pval_gsea=0.1, gene_sets=['MSigDB_Hallmark_2021','KEGG_2021_Human'],):

    explanation_df_norm_A_q = explanation_df_norm_A.query(f"annotation == {annotation_label_dict[region_A]} and cell_type == 17")
    explanation_df_norm_B_q = explanation_df_norm_B.query(f"annotation == {annotation_label_dict[region_B]} and cell_type == 17")

    explanation_df_sum_A = explanation_df_norm_A_q[attr_gene_list].sum(axis=1)
    explanation_df_sum_B = explanation_df_norm_B_q[attr_gene_list].sum(axis=1)

    
    # Define the maximum value and percentages for high and low thresholds
    max_value_A = explanation_df_sum_A.max()
    max_value_B = explanation_df_sum_B.max()

    # Calculate threshold values based on percentages
    threshold_value_high_A = max_value_A * percentage_high
    threshold_value_low_A = max_value_A * percentage_low
    
    threshold_value_high_B = max_value_B * percentage_high
    threshold_value_low_B = max_value_B * percentage_low
    
    row_mask_high_A = explanation_df_sum_A >= threshold_value_high_A
    row_mask_low_A = explanation_df_sum_A < threshold_value_low_A
    
    row_mask_high_B = explanation_df_sum_B >= threshold_value_high_B
    row_mask_low_B = explanation_df_sum_B < threshold_value_low_B
    
    explanation_df_norm_high_A = explanation_df_norm_A_q[row_mask_high_A]
    explanation_df_norm_low_A = explanation_df_norm_A_q[row_mask_low_A]
    
    explanation_df_norm_high_B = explanation_df_norm_B_q[row_mask_high_B]
    explanation_df_norm_low_B = explanation_df_norm_B_q[row_mask_low_B]
    
    print(len(explanation_df_norm_high_A), len(explanation_df_norm_low_A))
    print(len(explanation_df_norm_high_B), len(explanation_df_norm_low_B))
    
    above_threshold_df_A = explanation_df_norm_high_A[gene_list]
    below_threshold_df_A = explanation_df_norm_low_A[gene_list]
    
    above_threshold_df_B = explanation_df_norm_high_B[gene_list]
    below_threshold_df_B = explanation_df_norm_low_B[gene_list]
    
    
    adata_above_threshold_A = sc.AnnData(X=above_threshold_df_A.values)
    adata_below_threshold_A = sc.AnnData(X=below_threshold_df_A.values)

    adata_above_threshold_B = sc.AnnData(X=above_threshold_df_B.values)
    adata_below_threshold_B = sc.AnnData(X=below_threshold_df_B.values)
    
    
    adata_above_threshold_A.var_names = above_threshold_df_A.columns
    adata_below_threshold_A.var_names = below_threshold_df_A.columns
    
    adata_above_threshold_B.var_names = above_threshold_df_B.columns
    adata_below_threshold_B.var_names = below_threshold_df_B.columns
    
    adata_above_threshold_A.obs['dataset'] = 'above_A'
    adata_below_threshold_A.obs['dataset'] = 'below_A'
    
    adata_above_threshold_B.obs['dataset'] = 'above_B'
    adata_below_threshold_B.obs['dataset'] = 'below_B'

    combined_adata = adata_above_threshold_A.concatenate(adata_above_threshold_B)

    sc.tl.rank_genes_groups(combined_adata, groupby='dataset', method='wilcoxon')
    result_df = sc.get.rank_genes_groups_df(combined_adata, group='above_A')
    
    filtered_result_df = result_df[result_df['pvals_adj'] < pval_de]
    filtered_result_df.index = filtered_result_df.names
    
    top_upregulated = filtered_result_df.nlargest(20, 'logfoldchanges')
    top_downregulated = filtered_result_df.nsmallest(20, 'logfoldchanges')

    top_genes = pd.concat([top_upregulated, top_downregulated])
    top_genes.sort_values(by='logfoldchanges', inplace=True, ascending=False)
    top_genes = top_genes[~top_genes.index.duplicated(keep='first')]
    
    top_genes['pvals_adj'] = top_genes['pvals_adj'] + 1e-235
    names = top_genes.index
    logfoldchanges = top_genes['logfoldchanges']
    pvals_adj = -np.log10(top_genes['pvals_adj'])

    fig, axes = plt.subplots(1, 1, figsize=(15, 2))
    sns.heatmap(
        top_genes[['logfoldchanges']].T,
        cmap="RdBu_r",
        annot=False,
        fmt=".2f",
        center=0,
        cbar=True,
        linewidths=0.5,
        annot_kws={"size": 10},
        cbar_kws={"label": "logfoldchanges"},
    )
    plt.tight_layout()
    plt.title(f'Top DE genes')
    plt.show()

    ranked_result_df = filtered_result_df['logfoldchanges']
    ranked_result_df.index = filtered_result_df['names']
    ranked_result_df = pd.DataFrame(ranked_result_df)

    # # run prerank
    # # enrichr libraries are supported by prerank module. Just provide the name
    # # use 4 process to acceralate the permutation speed
    pre_res = gp.prerank(rnk=ranked_result_df,
                         gene_sets=gene_sets,
                         threads=4,
                         min_size=5,
                         max_size=1000,
                         permutation_num=1000, # reduce number to speed up testing
                         outdir=None, # don't write to disk
                         seed=6,
                         verbose=True, # see what's going on behind the scenes
                        )
    
    gsea_results_sorted = pre_res.res2d.sort_values(by='NES')
    gsea_results_sorted = gsea_results_sorted[gsea_results_sorted['NOM p-val'] < pval_gsea]

    # Bar plot for NES and Terms
    plt.figure(figsize=(8, 6))
    plt.barh(gsea_results_sorted['Term'], gsea_results_sorted['NES'], color='skyblue')
    plt.xlabel('NES (Normalized Enrichment Score)')
    plt.title(f'Enrichment Scores for Gene Sets')
    plt.grid(axis='x')  # Add gridlines along the x-axis
    plt.tight_layout()
    plt.show()

from upsetplot import UpSet, from_memberships
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_upset(adata, genes, interactions=False, title='UpSet Plot', min_subset_size=400, max_degree=6, show_percentages=True, 
               element_size=50, color_degrees=True, color_cat=None, color_no_cat=None, color_min_val=None, save_path=None, figsize=(20, 20), 
               dpi=200, save_upset_data_path=None, return_data=False, orientation='horizontal', sort_by='degree', sort_categories_by='cardinality', 
               subset_size='auto', sum_over=None, max_subset_size=None, max_subset_rank=None, min_degree=None, facecolor='auto', other_dots_color=0.18, 
               shading_color=0.05, with_lines=True, intersection_plot_elements=6, totals_plot_elements=2, show_counts='', include_empty_subsets=False):
    """
    Function to plot UpSet plot for an AnnData object.

    Parameters:
        adata (AnnData): AnnData object containing the data.
        genes (list): List of genes to consider for the plot.
        title (str): Title of the plot.
        min_subset_size (int): Minimum subset size to consider for the plot.
        max_degree (int): Maximum degree of intersections to consider.
        show_percentages (bool): Whether to show percentages in the plot.
        element_size (int): Size of elements in the plot.
        orientation (str): Orientation of the plot.
        sort_by (str): Criteria to sort subsets by.
        sort_categories_by (str): Criteria to sort categories by.
        subset_size (str): Criteria to calculate the size of subsets.
        sum_over (str or None): Field to sum over for subset sizes.
        max_subset_size (int or str): Maximum size of a subset to be shown.
        max_subset_rank (int): Limit to the top N ranked subsets.
        min_degree (int): Minimum degree of a subset to be shown.
        facecolor (str or float): Color for bar charts and active dots.
        other_dots_color (str or float): Color for shading of inactive dots.
        shading_color (str or float): Color for shading of odd rows.
        with_lines (bool): Whether to show lines joining dots in the matrix.
        intersection_plot_elements (int): Elements in the intersections plot.
        totals_plot_elements (int): Elements in the totals plot.
        show_counts (bool or str): Whether to label the intersection size bars.
        include_empty_subsets (bool): Whether to show all possible category combinations.
    """
    # Filter the AnnData object for positive cells

    # Ensure genes are in var_names
    genes_present = [gene for gene in genes if gene in adata.var_names]
    
    # Extract the expression data for relevant genes
    gene_expression = adata[:, genes_present].X.toarray()
    print(f"Shape of gene expression matrix: {gene_expression.shape}")  # Debug print
    
    # Determine positive cells and create a binary matrix
    binary_matrix = np.zeros((adata.shape[0], len(genes_present)), dtype=int)
    for i, gene in enumerate(genes_present):
        binary_matrix[:, i] = gene_expression[:, i] > 0

    print(f"Binary matrix shape: {binary_matrix.shape}")  # Debug print
    
    # Create a DataFrame from the binary matrix
    binary_df = pd.DataFrame(binary_matrix, columns=genes_present, index=adata.obs.index)

    # Transform binary DataFrame to count the occurrences of each combination
    binary_df['count'] = 1
    grouped_df = binary_df.groupby(genes_present).size().reset_index(name='count')
    print(f"Grouped DataFrame:\n{grouped_df.head()}")  # Debug print
    # Filter the grouped DataFrame for rows where both KSHV.ORF50 and KSHV.ORF57 are expressed
    orf50_orf57_positive = grouped_df[(grouped_df['KSHV.ORF50'] == 1) & (grouped_df['KSHV.ORF57'] == 1)]
    
    # Display the filtered DataFrame
    print(f"orf50_orf57_positive:\n{orf50_orf57_positive['count'].sum()}")  # Debug print print(orf50_orf57_positive)
    
    # Display the filtered DataFrame
    print(f"orf50_orf57_positive:\n{orf50_orf57_positive}")  # Debug print print(orf50_orf57_positive)

    # Convert to the required format for UpSet plot
    memberships = [set(genes_present[i] for i, val in enumerate(row) if val) for row in grouped_df[genes_present].values]

    # Generate the UpSet data
    upset_data = from_memberships(memberships, data=grouped_df['count'])
    print(f"UpSet data:\n{upset_data}")  # Debug print
    
    if isinstance(interactions, pd.core.indexes.multi.MultiIndex):
        try:
            upset_data = upset_data.reindex(interactions, fill_value=0)
            # upset_data.reindex(new_indices)
        except Exception as e:
            print(f"Interactions are incorrect. Please check the MultiIndex format and retry. Error: {e}")

    # Create and plot the UpSet plot
    upset_plot = UpSet(
        upset_data,
        orientation=orientation,
        sort_by=sort_by,
        sort_categories_by=sort_categories_by,
        subset_size=subset_size,
        sum_over=sum_over,
        min_subset_size=min_subset_size,
        max_subset_size=max_subset_size,
        max_subset_rank=max_subset_rank,
        min_degree=min_degree,
        max_degree=max_degree,
        facecolor=facecolor,
        other_dots_color=other_dots_color,
        shading_color=shading_color,
        with_lines=with_lines,
        element_size=element_size,
        intersection_plot_elements=intersection_plot_elements,
        totals_plot_elements=totals_plot_elements,
        show_counts=show_counts,
        show_percentages=show_percentages,
        include_empty_subsets=include_empty_subsets
    )

    if return_data:
        return upset_plot
        
    if color_degrees:
        upset_plot.style_subsets(min_degree=1, facecolor="k")
        upset_plot.style_subsets(min_degree=2, facecolor="blue")
        upset_plot.style_subsets(min_degree=3, facecolor="green")
        upset_plot.style_subsets(min_degree=4, facecolor="red")
        upset_plot.style_subsets(min_degree=5, facecolor="purple")

    if isinstance(color_cat, list):
        upset_plot.style_subsets(present=color_cat, facecolor="blue", label=None)

    if isinstance(color_no_cat, list):
        upset_plot.style_subsets(absent=color_no_cat, facecolor="green", label=None)
        
    if isinstance(color_min_val, int):
        upset_plot.style_subsets(min_subset_size=color_min_val, facecolor="red", label=None)

    plt.figure(figsize=figsize, dpi=dpi)
    upset_plot.plot()
    plt.title(title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    
    plt.show()

