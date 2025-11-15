import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import seaborn as sns
import datetime
from pathlib import Path
from shapely.ops import unary_union
from scipy.spatial import KDTree
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib import gridspec
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, box
from tqdm import tqdm
import concurrent.futures



# Path to the current script
BASE_DIR = Path(__file__).resolve().parent

# Function to create hexagonal grid cells
def create_hexagonal_grid(gdf, hex_size):
    hexagons = []
    sqrt3 = np.sqrt(3)

    # Get the bounding box of the polygon
    minx, miny, maxx, maxy = gdf.total_bounds
    bbox = box(minx, miny, maxx, maxy)
    # Generate hexagons within the bounding box
    for row in np.arange(miny, maxy, hex_size * sqrt3):
        for col in np.arange(minx, maxx, hex_size * 1.5):
            # Adjust y-coordinate for staggered rows
            y_offset = (col % 2) * (hex_size * sqrt3 / 2)
            hexagon = Polygon([
                (col + hex_size * np.cos(np.pi / 3 * i), row + y_offset + hex_size * np.sin(np.pi / 3 * i))
                for i in range(6)
            ])
            intersects = bbox.intersects(hexagon)
            if intersects:
                hexagons.append(hexagon)
    return hexagons

def assign_niche(hex_gdf, category_column, core_adata_gdf):
    core_adata_gdf['centroid'] = core_adata_gdf.geometry.centroid
    centroid_gdf = core_adata_gdf[['cell_id', 'centroid', category_column]].copy()
    centroid_gdf.rename(columns={'centroid': 'geometry'}, inplace=True)    
    cells_in_polygons = gpd.sjoin(centroid_gdf, hex_gdf, how="inner", predicate="within")
    cluster_counts = cells_in_polygons.groupby(['index_right', category_column], observed=True).size().reset_index(name='count')
    majority_cluster = cluster_counts.loc[cluster_counts.groupby('index_right')['count'].idxmax()]
    hex_gdf[category_column] = hex_gdf.index.map(majority_cluster.set_index('index_right')[category_column])
    hex_gdf.dropna(subset=[category_column], inplace=True)    
    return hex_gdf

def count_vessels(vessels, hex_gdf):
    vessels_in_polygons = gpd.sjoin(vessels, hex_gdf, how="inner", predicate="intersects")
    cluster_counts = vessels_in_polygons.groupby(['index_right'], observed=True).size().reset_index(name='count')
    hex_gdf['vessel_count'] = hex_gdf.index.map(cluster_counts.set_index('index_right')['count'])
    hex_gdf['vessel_count'] = hex_gdf['vessel_count'].fillna(0)
    return hex_gdf

def get_adata(path_block, category_column, vessels, cell_boundaries):
    vessels_subset = vessels[vessels['path_block_core'] == path_block].copy()
    core_adata_gdf = cell_boundaries[cell_boundaries['path_block_core'] == path_block].copy()

    # Process vessel geometries
    vessels_subset['geometry'] = vessels_subset['geometry'].apply(
        lambda geom: unary_union([g.buffer(5) for g in geom.geoms]).buffer(-5)
        if geom.geom_type in ['MultiPolygon', 'GeometryCollection'] else geom.buffer(5).buffer(-5)
    )
    return core_adata_gdf, vessels_subset

def process_path_block(path_block, vessels, cell_boundaries, output_dir):
    category_column = 'niche_with_tumor_proximity'
    hex_size = 42  # Size of the hexagon (radius)
    output_folder = output_dir/f'hexgrid_c{category_column}_s{hex_size}'
    output_folder.mkdir(parents=True, exist_ok=True)

    try:
        core_adata_gdf, vessels_subset = get_adata(path_block, category_column, vessels, cell_boundaries)
        # Create hexagonal grid for the polygon
        hexagonal_cells = create_hexagonal_grid(core_adata_gdf, hex_size)

        # Create a new GeoDataFrame for the hexagonal grid
        hex_gdf = gpd.GeoDataFrame(geometry=hexagonal_cells)
        hex_gdf = assign_niche(hex_gdf, category_column, core_adata_gdf)
        hex_gdf = count_vessels(vessels_subset, hex_gdf)
        hex_gdf['path_block_core'] = path_block
        hex_gdf['Stage'] = core_adata_gdf['Stage'].unique()[0]
        with open(output_folder/f'{path_block}.pkl', 'wb') as f:
            pickle.dump(hex_gdf, f)
    except Exception as e:
        print(f"Error processing {path_block}: {e}")

def compile_all_vessel_counts(input_folder, output_dir):
    output_dir = Path(output_dir)
    all_dfs = []
    for pkl_file in input_folder.iterdir():
        if pkl_file.suffix == '.pkl':
            df = pd.read_pickle(pkl_file)
            all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    return combined_df

def process_vessel_counts(vessel_counts):
    hex_size = 42
    area = (3 * np.sqrt(3) / 2) * (hex_size ** 2)
    category_column = 'niche_with_tumor_proximity'
    vessel_counts['normalized_count'] = vessel_counts['vessel_count']/area
    vessel_counts_per_core_per_niche = vessel_counts.groupby(['path_block_core', category_column, 'Stage'], observed=False)['normalized_count'].mean().reset_index()
    vessel_counts_per_core_per_niche['normalized_count'].fillna(0, inplace=True)
    return vessel_counts_per_core_per_niche

def process_all_cores():
    adata_path = BASE_DIR/'adata.h5ad'
    vessels_path = BASE_DIR/'vessels.pkl'
    cell_boundaries_path = BASE_DIR/'cell_boundaries.pkl'
    adata = sc.read_h5ad(adata_path)
    vessels = pd.read_pickle(vessels_path)
    cell_boundaries = pd.read_pickle(cell_boundaries_path)
    for path_block in tqdm(adata.obs['path_block_core'].unique().tolist()):
        output_dir = Path(BASE_DIR/'hexgrid_output')
        output_dir.mkdir(parents=True, exist_ok=True)
        process_path_block(path_block, vessels, cell_boundaries, output_dir)

def normalize_counts():
    input_folder = Path(BASE_DIR/'hexgrid_output/hexgrid_cniche_with_tumor_proximity_s42')
    output_folder = Path(BASE_DIR/'hexgrid_output')
    combined_vessel_counts = compile_all_vessel_counts(input_folder, output_folder)
    processed_counts = process_vessel_counts(combined_vessel_counts)
    processed_counts.to_pickle(output_folder/'vessel_counts.pkl')


if __name__ == "__main__":
    # process_all_cores()
    normalize_counts()
