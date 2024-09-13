#afshin.ghahramania@des.qld.gov.au

"""
This script is designed to process Digital Elevation Model (DEM) data by merging multiple DEM files into a single raster
 and then using a shapefile to clip this raster to a specific area. The process begins by loading the DEM files and checking
their Coordinate Reference System (CRS). It then ensures the shapefile is in the same CRS as the DEMs. If necessary,
the shapefile is reprojected to match the CRS of the DEMs. After merging the DEMs into a single raster, the script
checks if any of the shapefile's geometries overlap with the raster bounds. If they do, it applies a mask to clip
the raster based on the shapefile geometries. The final clipped raster is saved as a new TIFF file. This script
is useful for geographic data processing where rasters need to be aligned and clipped according to specific geographic features.
"""




import rasterio
from rasterio.merge import merge
import rasterio.mask
from rasterio.io import MemoryFile
import glob
import os
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import box


# Function to open and read a DEM file
def read_dem_file(dem_file):
    src = rasterio.open(dem_file)
    return src


# Step 1: Load DEMs and Determine Their Common CRS
# Specify the directory where your DEM files are stored
dem_files_directory = r'D:\Elvis_data\Data\QLD Government\DEM\1 Metre'
dem_files = glob.glob(os.path.join(dem_files_directory, '*.tif'))

# Assuming all DEMs share the same CRS, check the CRS of the first DEM
with rasterio.open(dem_files[0]) as first_dem:
    dem_crs = first_dem.crs
print(f"DEM CRS: {dem_crs}")

# Step 2: Load Shapefile and Check Its CRS
shapefile_path = r'D:\Elvis_data\Grand_Canyont.shp'
shapefile = gpd.read_file(shapefile_path)
print(f"Original Shapefile CRS: {shapefile.crs}")

# Step 3: Reproject Shapefile to Match DEM CRS, If Necessary
if shapefile.crs != dem_crs:
    shapefile = shapefile.to_crs(dem_crs)
    print("Shapefile has been reprojected to match the DEM's CRS.")
else:
    print("Shapefile is already in the DEM's CRS.")

# Step 4: Merge DEM Files
with ThreadPoolExecutor(max_workers=4) as executor:
    src_files_to_mosaic = list(
        tqdm(executor.map(read_dem_file, dem_files), total=len(dem_files), desc='Reading DEM files in parallel'))

print('Merging DEM files...')
mosaic, out_trans = merge(src_files_to_mosaic)

# Update the metadata for creating an in-memory file with the merged DEM
out_meta = src_files_to_mosaic[0].meta.copy()
out_meta.update({"driver": "GTiff", "height": mosaic.shape[1], "width": mosaic.shape[2], "transform": out_trans})

# Step 5: Apply Mask Using Reprojected Shapefile Geometries
geometries = [geom for geom in shapefile.geometry]

with MemoryFile() as memfile:
    with memfile.open(**out_meta) as tmp_dst:
        tmp_dst.write(mosaic)
        raster_polygon = box(*tmp_dst.bounds)

        # Check if any geometry falls within the raster bounds
        any_within_bounds = any(raster_polygon.intersects(geom) for geom in geometries)
        print("Any geometry falls within raster bounds:", any_within_bounds)

        if not any_within_bounds:
            print("No geometries overlap with the raster; skipping mask operation.")
            exit(1)

        # Apply the mask using geometries
        masked_mosaic, masked_transform = rasterio.mask.mask(tmp_dst, geometries, crop=True)

# Update the metadata to reflect the new dimensions and transform after masking
out_meta.update({"transform": masked_transform, "height": masked_mosaic.shape[1], "width": masked_mosaic.shape[2]})

# Define path for the merged and masked file
output_file = r'D:\Elvis_data\Data\QLD Government\DEM\merged_dem_1m.tif'

# Write the masked mosaic to a new TIFF file
with rasterio.open(output_file, "w", **out_meta) as dest:
    dest.write(masked_mosaic)

print('Writing merged and masked file completed successfully.')

# Close the source files
for src in src_files_to_mosaic:
    src.close()

print("Merging and masking completed successfully.")
