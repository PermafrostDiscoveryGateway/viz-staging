"""
An updated & simplified version of the prepare_footprints.py script for delta
"""

from pathlib import Path
import json
import geopandas as gpd
import warnings

# Where is all the original input data, before any staging or processing
input_dir = '/scratch/bbki/kastanday/maple_data_xsede_bridges2/glacier_water_cleaned_shp/'
# Where are the original footprints?
footprint_dir = '/scratch/bbki/kastanday/maple_data_xsede_bridges2/outputs/footprints/footprints_new'
# Where should we save the prepared footprints?
footprint_staged_dir = '/scratch/bbki/thiessenbock/pdg/staged_footprints'
# Where should we save the JSON list of input files missing footprints?
missing_footprints_path = '/scratch/bbki/thiessenbock/pdg/files-missing-footprints.json'


def get_fp_path(input_path):
    """Given the input path, get the staged footprint path"""
    # Remove the base directory from the input path
    input_path_rel = str(input_path).removeprefix(input_dir).removeprefix('/')
    fp_path = Path(footprint_staged_dir).joinpath(
        input_path_rel).with_suffix('.gpkg')
    return fp_path


def dir_to_filename(dir):
    """Change the scene directory from input path to file name for
    footprints"""
    return 'selection_' + dir.removesuffix('_iwp') + '.shp'


def id_to_name(id):
    """Convert the scene ID from input filename to 'Name' attribute in
    footprint shapefiles"""
    return '_'.join(id.split('_')[:-2])


def contains_multipolygons(gdf):
    types = list(set(gdf.geom_type))
    return 'MultiPolygon' in types


# Find input files
print('Finding input files')
input_file_list = Path(input_dir).rglob('*.shp')
print('Input files found')

inputs_missing_footprints = []
for i, path in enumerate(input_file_list):
    print(f'🏁 Looking for footprint for file {i+1}')
    path_parts = path.parts
    scene_id = path_parts[-2]
    scene_dir = path_parts[-3]
    region = path_parts[-4]
    # footprints follow format FOOTPRINT_DIR/REGION/SHORT_ID.shp
    footprint_file = Path(footprint_dir).joinpath(
        region).joinpath(dir_to_filename(scene_dir))
    if not (footprint_file.exists()):
        print('🤷‍♀️ No footprint file found')
        inputs_missing_footprints.append(str(path))
    else:
        # each footprint shapefile contains footprints for multiple inputs
        all_fps = gpd.read_file(footprint_file)
        fp = all_fps[all_fps['Name'] == id_to_name(scene_id)]
        if len(fp) == 0:
            print('🤷‍♀️ No footprint ID found')
            inputs_missing_footprints.append(str(path))
        else:
            print('✅ Found a footprint')
            if contains_multipolygons(fp):
                fp = fp.explode()
            # Create the path where the footprint should be saved
            footprint_staged_path = get_fp_path(path)
            # Make any necessary parent directories
            footprint_staged_path.parents[0].mkdir(parents=True, exist_ok=True)
            # Save the single-polygon vector file to the filepath
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FutureWarning)
                fp.to_file(footprint_staged_path)
    # For testing
    # if i == 10:
    #     break

# Save a list of all the files that do not have footprints
with open(missing_footprints_path, 'w') as f:
    json.dump(inputs_missing_footprints, f)
