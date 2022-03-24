import logging
import logging.config
import argparse
from .TileStager import TileStager

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
        '-c', '--config',
        help='Path to configuration JSON file',
        default='config.json',
        type=str
    )

    args = parser.parse_args()

    stager = TileStager(args.config)
    stager.stage_all()

    logging.info('Done')
    logging.shutdown()
    exit(0)
