import os
import json

"""
This script prepares a list of input paths for the Ice Wedge Polygon project
(IWP). It takes a list of main input directories and a list of secondary,
preferred input directories. For every path that it finds in the main
directories, it looks for a file with the exact file name in the preferred
directories. If it finds one, it uses the preferred directory path. If it
doesn't find one, it uses the main directory path. This is required for the IWP
project because there exists a 'water_clipped' directory, with shapefiles that
are identical to those in the other directories, except that the polygons are
removed where there is a water body. In the workflow, we want to use the water
masked version of the data over the unmasked version.
"""

filename_path_list = 'input_paths.json'

# Will be prepended to all fo the main & preferred input dir paths
base_dir = '/home/pdg/data/ice-wedge-polygon-data/version 01/'

# Directories that contain all the data
main_dirs = [
    'alaska',
    'high_ice'
]

# Directories that contain better/newer/more desired versions of the same data
# in main dir
preferred_dirs = [
    'water_clipped'
]
ext = '.shp'

# Prepend the base dir
main_dirs = [os.path.join(base_dir, d) for d in main_dirs]
preferred_dirs = [os.path.join(base_dir, d) for d in preferred_dirs]


def get_file_paths(dir, ext):
    """
    Recursively find all files in the given directory with the given extension.
    """
    files = []
    for root, dirs, filenames in os.walk(dir):
        for filename in filenames:
            if filename.endswith(ext):
                files.append(os.path.join(root, filename))
    return files


# Recursively find all files in the main directories with the correct extension
# and add them to a list
main_files = []
for main_dir in main_dirs:
    files = get_file_paths(main_dir, ext)
    main_files.extend(files)

# Extract just the filename and ext, not the path, from the main files
main_file_basenames = [os.path.basename(f) for f in main_files]

print(f'Found {len(main_files)} files with ext {ext} in main directories.')

# Repeat the same process for the preferred directories
preferred_files = []
for preferred_dir in preferred_dirs:
    files = get_file_paths(preferred_dir, ext)
    preferred_files.extend(files)
preferred_files_basenames = [os.path.basename(f) for f in preferred_files]

print(
    f'Found {len(preferred_files)} files with ext {ext} in preferred '
    f'directories.')

# When a main file basename exists in the list of preferred files, replace
# the main file path with the preferred file path
final_path_list = main_files.copy()
replacement_count = 0
for main_file_basename in main_file_basenames:
    if main_file_basename in preferred_files_basenames:
        ind = main_file_basenames.index(main_file_basename)
        ind_p = preferred_files_basenames.index(main_file_basename)
        final_path_list[ind] = preferred_files[ind_p]
        # print(f'Replced {main_files[ind]} with {preferred_files[ind_p]}')
        replacement_count += 1

# Write the list of paths to a JSON file
with open(filename_path_list, 'w') as f:
    json.dump(final_path_list, f)

print(f'Wrote {len(final_path_list)} paths to {filename_path_list}.')
print(f'Replaced {replacement_count} paths with the preferred version.')
# Did all the preferred files exist in the main directories?
if len(preferred_files) != replacement_count:
    print(
        f'WARNING: There were {len(preferred_files)} preferred files, but '
        f'only {replacement_count} main file paths were replaced.')


# #############################################################################

# which files in the preferred files list do not have a corresponding file in
# the main files list?

# no_match = [
#   f for f in preferred_files_basenames if f not in main_file_basenames]

# no_match_full_path = []
# for f in preferred_files_basenames:
#     if f not in main_file_basenames:
#         ind_p = preferred_files_basenames.index(f)
#         no_match_full_path.append(preferred_files[ind_p])

# [f for f in no_match_full_path if 'canada' not in f]
