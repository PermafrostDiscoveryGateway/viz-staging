import logging
import logging.config
import argparse
from stager import TileStager

# Set up logging (TODO: move to config file)

# log_dict = {
#     'version': 1,
#     'disable_existing_loggers': False,
#     'formatters': {
#         'standard': {
#             'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
#         },
#     },
#     'handlers': {
#         'default': {
#             'level': 'INFO',
#             'formatter': 'standard',
#             'class': 'logging.StreamHandler',
#         },
#         'file_handler': {
#             'level': 'INFO',
#             'filename': 'viz-staging.log',
#             'class': 'logging.FileHandler',
#             'formatter': 'standard'
#         }
#     },
#     'loggers': {
#         '': {
#             'handlers': ['file_handler'],
#             'level': 'INFO',
#             'propagate': True
#         },
#     }
# }

# logging.config.dictConfig(log_dict)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Stages vector files for the PDG viz tiling pipeline'
    )
    parser.add_argument(
        '-l', '--logging-config',
        help='Path to logging configuration file',
        default='logging.yaml'
    )
    parser.add_argument(
        '-i', '--input',
        help='Input shapefile',
        required=True
    )
    parser.add_argument(
        '-o', '--output',
        help='Output directory',
        required=True
    )
    parser.add_argument(
        '-ie', '--inputext',
        help='Output file extension',
        default='.shp'
    )
    parser.add_argument(
        '-oe', '--outputext',
        help='Output file extension',
        default='.shp'
    )
    parser.add_argument(
        '-c', '--crs',
        help='Input shapefile CRS'
    )
    parser.add_argument(
        '-s', '--simplify',
        help='Simplify tolerance',
        default=0.000001
    )
    parser.add_argument(
        '-t', '--tms',
        help='Tile Matrix Set identifier',
        default='WorldCRS84Quad'
    )
    parser.add_argument(
        '-z', '--zoom',
        help='Zoom level',
        default=14
    )
    parser.add_argument(
        '-f', '--format',
        help='Tile path format',
        default=['z', 'x', 'y']
    )
    parser.add_argument(
        '-m', '--summary',
        help='Summary file',
        default='staging-summary.csv'
    )

    args = parser.parse_args()

    stager = TileStager(
        input_dir=args.input,
        input_ext=args.inputext,
        output_dir=args.output,
        output_ext=args.outputext,
        input_crs=args.crs,
        simplify_tolerance=args.simplify,
        tms_identifier=args.tms,
        zoom_level=args.zoom,
        tile_path_structure=args.format,
        summary_path=args.summary
    )

    stager.stage_all()

    logging.info('Done')
    logging.shutdown()
    exit(0)
