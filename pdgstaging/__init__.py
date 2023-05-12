from .ConfigManager import ConfigManager
from .Grid import TMSGrid
from .TilePathManager import TilePathManager
from .TileStager import TileStager
from .Deduplicator import *

# # configure logger
# import logging
# logger = logging.getLogger("logger")
# # prevent logging statements from being printed to terminal
# logger.propagate = False
# # set up new handler
# handler = logging.FileHandler("/tmp/log.log")
# formatter = logging.Formatter(logging.BASIC_FORMAT)
# handler.setFormatter(formatter)
# logger.addHandler(handler)
# logger.setLevel(logging.INFO)