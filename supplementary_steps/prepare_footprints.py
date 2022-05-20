"""
This python script takes footprint files for the Ice Wedge Polygon (IWP)
project and prepares them for use in the deduplication step of the the PDG
viz-workflow. For details on deduplication, see:
https://github.com/PermafrostDiscoveryGateway/viz-staging/blob/develop/docs/deduplication.md

Although this process is very specific to the IWP project, it is possible to
configure certain parameters to work with other projects, should they have a
similar footprint file naming convention.

INPUT DATA:
------------

The footprint filenames that are prepared by this script are provided with the
following naming convention: `selection_<directory_id>.shp`. Where
`directory_id` is gives a partial name for the subdirectory where the
associated IWP input data can be found. The matching IWP subdirectory name
excludes the `selection_` prefix, but includes an the `_iwp` suffix, e.g.
`<directory_id>_iwp`. This script searches recursively for IWP input
subdirectories that correspond to the footprint file name as described.

The input footprint files contain one feature for each footprint. Each of these
features correspond to a particular IWP input file, and each feature has a
'Name' property that is a partial match to the IWP input file name. The IWP
input file name is this 'Name' (file ID) property plus a variable two part
suffix following the format `<file_id>_suffix1_suffix2.shp`, e.g.
`<file_id>_u16rf3413_pansh.shp`. This script matches the file_ids in the
footprint file to the IWP input file names.

The second part of the file ID (when ID is split by underscore), is the
Date & Time that the IWP image was captured.

Here is an example:

/footprint_input_dir/
    |
    |- selection_<directory_id>.shp
        | --------------------------
        | geometry | Name       | ...
        | ...      | <file_id1> | ... <- **FOOTPRINT 1**
        | ...      | <file_id2> | ... <- **FOOTPRINT 2**
        ...

... MATCHES TO ...

/IWP_input_dir/
    |
    | subdir1/
        | subdir2/
            ...
            | <directory_id>_iwp/
                |- <file_id1>_suffix1_suffix2.shp <- **IWP INPUT 1**
                |- <file_id2>_suffix1_suffix2.shp <- **IWP INPUT 2**
                ...

OUTPUT DATA:
-------------
The script matches footprint features to IWP input files as described above,
then re-saves each footprint as a single file with the same name and directory
structure as the IWP input files, as required by the PDG viz-workflow. It also
parses the date from the file ID and saves it as a new property, 'Date', in the
footprint file.

For details on the required format of footprints, see documentation here:
https://github.com/PermafrostDiscoveryGateway/viz-staging/blob/develop/docs/footprints.md.
Details of footprint matches or lack thereof are saved to several JSON files.
"""

import os
import warnings
import json
import geopandas as gpd
from pdgstaging import TileStager

# Some options that could change on a per-project basis
options = {
    'base_dir': '/home/pdg/data/ice-wedge-polygon-data/version 01',
    'subdir_suffix': '_iwp',
    'footprint_file_prefix': 'selection_',
    'dir_vector_input_data': '',
    'dir_footprints_in': 'footprints/original_footprints',
    'dir_footprints_out': 'footprints/staged_footprints',
    'ext_footprints_in': '.shp',
    'ext_footprints_out': '.gpkg',
    'prop_file_id': 'Name',
    'prop_date': 'Date',
    # Directories to skip when searching for IWP input data that matches
    # footprint file names
    'dirs_skip': [
        'web_tiles_1TB_testrun',
        'high_ice/russia/russia',
        'high_ice/medium_ice/alaska_m',
        'high_ice/medium_ice/russia_m',
        'footprints'
    ],
    # Where to save records of footprint matches and failures
    'filename_unmatched_footprints':
        'footprints/footprint_files_unmatched_to_subdirs.json',
    'filename_matched_footprints':
        'footprints/footprint_files_matched_to_subdirs.json',
    'filename_multimatch_footprints':
        'footprints/footprint_files_multimatch_to_subdirs.json',
    'filename_unmatched_footprint_features':
        'footprints/footprint_features_unmatched_to_files.json',
    'filename_matched_footprint_features':
        'footprints/footprint_features_matched_to_files.json'
}

# Add base directory to options
path_opts = [
    'dir_vector_input_data',
    'dir_footprints_in',
    'dir_footprints_out',
    'filename_unmatched_footprints',
    'filename_matched_footprints',
    'filename_multimatch_footprints',
    'filename_unmatched_footprint_features',
    'filename_matched_footprint_features']

for o in path_opts:
    options[o] = os.path.join(options['base_dir'], options[o])

for i in range(0, len(options['dirs_skip'])):
    options['dirs_skip'][i] = os.path.join(
        options['base_dir'], options['dirs_skip'][i])


def get_base_name(path):
    """
    Get the base name of a file, without the extension
    """
    return os.path.basename(path).split('.')[0]


def date_from_id(string):
    """
    Parse date from IWP file name
    """
    # Split string by underscore
    parts = string.split('_')
    # Return the second part, represents date & time
    return parts[1]


def id_from_input_path(input):
    """
        Get just the IWP file 'ID' code from the full path name that
        includes a two-part suffix
    """
    input = get_base_name(input)
    parts = input.split('_')
    parts = parts[:-2]
    input = '_'.join(parts)
    return input


def subdir_from_footprint_path(footprint_path):
    """
        Get the IWP subdirectory name from the footprint file path. The
        sub-directory of the input file is the the filename of the matching
        footprint, minus the prefix 'selection_', plus the suffix '_iwp'
    """
    subdir = get_base_name(footprint_path)
    subdir = subdir.removeprefix(options['footprint_file_prefix'])
    subdir += options['subdir_suffix']
    return subdir


def get_input_subdirs(rootdir='.', dirs_skip=[], recursive=True):
    """
    List all subdirectories in a directory that include the '_iwp' suffix,
    excluding the ones contained within in the dirs_skip list.
    """
    dirs = []
    for item in os.scandir(rootdir):
        if item.is_dir() and (item.path not in dirs_skip):
            if item.name.endswith(options['subdir_suffix']):
                dirs.append(item.path)
            if recursive:
                subdirs = get_input_subdirs(item.path, dirs_skip, recursive)
                if subdirs:
                    dirs.extend(subdirs)
    return dirs


# Create a tile stager with the location of the input IWP data files, and the
# *DESIRED* location of the footprint files. dir_footprints is where we will
# save the footprints that are prepared for staging.
tileStager = TileStager({
    'dir_input': options['dir_vector_input_data'],
    'dir_footprints': options['dir_footprints_out'],
    'ext_input': options['ext_footprints_in'],
    'ext_footprints': options['ext_footprints_out']
})

# To help create and parse paths
pathManager = tileStager.tiles

# To get options from the configuration for the workflow
config = tileStager.config

# Add the directory where the footprints are currently stored
pathManager.add_base_dir(
    name='footprints_original',
    dir_path=options['dir_footprints_in'],
    ext=options['ext_footprints_in'],
)

# Get the paths to all of the original footprint files
footprint_paths = pathManager.get_filenames_from_dir('footprints_original')

dir_input = config.get('dir_input')
ext_input = config.get('ext_input')

# Get a list of all possible input data sub-directories
all_input_subdirs = get_input_subdirs(dir_input, options['dirs_skip'])

# For record keeping
directory_matches = []
file_matches = []

for footprint_path in footprint_paths:

    # Get the dir that contains the corresponding input files
    subdir_name = subdir_from_footprint_path(footprint_path)

    # Find the dir in all_input_subdirs that matches the subdir_name
    subdir_input = [dir for dir in all_input_subdirs if dir.split(
        os.sep)[-1] == subdir_name]

    # Record the match (or lack thereof)
    directory_match = {
        'original_footprint_path': footprint_path,
        'status': 'no_match',
        'match': None
    }
    if len(subdir_input) > 1:
        directory_match['status'] = 'multiple_matches'
        directory_match['match'] = subdir_input
    elif len(subdir_input) == 1:
        subdir_input = subdir_input[0]
        directory_match['status'] = 'matched'
        directory_match['match'] = subdir_input

    directory_matches.append(directory_match)

    # Stop here if not 1 and only 1 match
    if directory_match['status'] != 'matched':
        continue

    # Get the paths & IDs of input files in the matched sub-directory
    pathManager.add_base_dir(subdir_input, subdir_input, ext_input)
    subset_iwp_files = pathManager.get_filenames_from_dir(subdir_input)
    subset_ids = [id_from_input_path(i) for i in subset_iwp_files]

    # Read the footprint GDF, split into 1 GeoDataFrame per row
    fp = gpd.read_file(footprint_path)
    fp_list = [v for k, v in fp.groupby('Name', as_index=False)]

    # For each GDF (row)
    for fp_row in fp_list:

        # Set the date
        fp_row.reset_index(drop=True, inplace=True)
        fp_id = fp_row.Name[0]
        date = date_from_id(fp_id)
        fp_row['Date'] = date

        file_match = {
            'original_footprint_path': footprint_path,
            'footprint_id': fp_id,
            'date': date,
            'status': 'no_match',
            'match': None
        }

        # Find the corresponding IWP file
        try:
            # Match footprint ID (from GDF) to IWP ID (from filename)
            id_index = subset_ids.index(fp_id)
            matching_file = subset_iwp_files[id_index]

            # Create the file path expected by the viz-workflow
            footprint_path_staging = config.footprint_path_from_input(
                matching_file, check_exists=False)

            # Make any necessary parent directories
            pathManager.create_dirs(footprint_path_staging)

            # Save the single-polygon vector file to the filepath
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FutureWarning)
                fp_row.to_file(footprint_path_staging)

            # Record the match
            file_match['status'] = 'matched'
            file_match['match'] = matching_file
            file_matches.append(file_match)

        except ValueError:
            file_matches.append(file_match)

# Filter records and map to output filenames
record_map = {
    options['filename_unmatched_footprints']:
        [d for d in directory_matches if d['status'] == 'no_match'],
    options['filename_matched_footprints']:
        [d for d in directory_matches if d['status'] == 'matched'],
    options['filename_multimatch_footprints']:
        [d for d in directory_matches if d['status'] == 'multiple_matches'],
    options['filename_unmatched_footprint_features']:
        [f for f in file_matches if f['status'] == 'no_match'],
    options['filename_matched_footprint_features']:
        [f for f in file_matches if f['status'] == 'matched']
}

# Write all of the records to the output files as JSON
for filename, records in record_map.items():
    with open(filename, 'w') as f:
        json.dump(records, f, indent=2)
