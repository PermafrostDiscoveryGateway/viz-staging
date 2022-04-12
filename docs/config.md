# Configuration

The staging (and rasterization) process is configured with a python `dict` or a `.json` file. The `dict` or `json` path is passed to the `TileStager` class.

Configuration options are detailed in the `ConfigManager` class, see `help(pdgstaging.ConfigManager)` for details and `pdgstaging.ConfigManager.defaults` for default config values.

For more in depth explanation of some of the configuration options, see:
- [deduplication](deduplication.md)
- [output path structure](tile_path_structure.md)
