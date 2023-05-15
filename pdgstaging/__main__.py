import logging
import logging.config
import argparse
from pdgstaging import TileStager

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
