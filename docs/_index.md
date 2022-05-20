#  PDG Staging

## Configuration

The staging (and rasterization) process is configured with a python `dict` or a `.json` file. The `dict` or `json` path is passed to the `TileStager` class.

Configuration options are detailed in the `ConfigManager` class, see `help(pdgstaging.ConfigManager)` for details and `pdgstaging.ConfigManager.defaults` for default config values.

For more in depth explanation of some of the configuration options, see:
- [deduplication](deduplication.md)
- [footprints for deduplication](footprints.md)
- [output path structure](tile_path_structure.md)

## Supplementary steps

Some projects will require additional processing for the visualization workflow. Scripts
that are more or less project-specific are stored in the `supplementary_steps` of this repo.

The `prepare_footprints.py` file prepares the footprint files from the Ice
Wedge Polygon project to be used to deduplicate input vector files. See the
documentation in that file for a detailed overview.