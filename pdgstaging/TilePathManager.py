import os
import re
import morecantile
import geopandas as gpd
from shapely.geometry import box


class TilePathManager():
    """
        The Tile & Path Manager class handles creating and parsing paths to
        geospatial data that have been tiled in accordance with a TileMatrixSet
        (tms). Using the class, paths to tiled files follow a standardized
        structure that incorporate the tile's x, y and z indices. Also includes
        helper methods for managing paths in general (e.g. creating the base
        directory for a path for a file, if the dir doesn't exist), and for
        working with the morecantile library in general (e.g. get the bounding
        box for a tile as a dict).
    """

    def __init__(self,
                 tms_id='WorldCRS84Quad',
                 path_structure=('style', 'tms', 'z', 'x', 'y'),
                 base_dirs={}
                 ):
        """
            Create a TilePathManager object.

            Parameters
            ----------
            tms_id : str, optional
                The Tile Matrix Set ID. Must be one of the IDs used by the
                morecantile library. Defaults to 'WorldCRS84Quad'.
            path_structure : list or tuple, optional
                A list of strings that represent the directory structure of
                last segment of the path that uses the tms (TileMatrixSet),
                style (layer/statistic), x index (TileCol), y index (TileRow),
                and z index (TileMatrix) of the tile. By default, the path will
                be in the format of
                {TileMatrixSet}/{Style}/{TileMatrix}/{TileCol}/{TileRow}.ext,
                configured as ('style', 'tms', 'z', 'x', 'y'). The x, y, z
                indices must always be last in the list.
            base_dirs : dict, optional
                A dictionary of base directories that can be referenced in some
                of the methods in this class. The keys are the names of the
                base directories and the values are dictionaries with the keys
                'path' and 'ext'. The 'path' key is the directory path and the
                'ext' key is the file extension including the '.', e.g.:
                {'base_dir_name': {'path': '/path/to/dir', 'ext': '.tif'}}
        """

        isValid, error_message, error_type = self.validate_tms_id(tms_id)
        if not isValid:
            raise error_type(error_message)

        isValid, error_message, error_type = self.validate_path_structure(
            path_structure)
        if not isValid:
            raise error_type(error_message)

        self.tms_id = tms_id
        self.tms = morecantile.tms.get(tms_id)
        self.crs = self.tms.crs
        self.tile_type = morecantile.commons.Tile

        self.path_structure = path_structure
        self.base_dirs = base_dirs or {}

    def validate_tms_id(self, tms_id):
        """
            Check if the TMS is valid

            Parameters
            ----------
            tms_id : str
                The Tile Matrix Set ID.

            Returns
            -------
            (isValid, error_message, error_type) : tuple
                Returns a tuple with a boolean indicating whether the ID is
                valid, a string with a human-readable error message, and an
                error type that may be raised.
        """
        error_message = ''
        error_type = None
        isValid = False
        if not isinstance(tms_id, str):
            error_message = 'Tile Matrix Set ID must be a string'
            error_type = TypeError
        elif tms_id not in morecantile.tms.list():
            error_message = f'{tms_id} is not a valid tms_id. ID must ' + \
                f'be one of the following: {morecantile.tms.list()}'
            error_type = ValueError
        else:
            isValid = True
        return (isValid, error_message, error_type)

    def validate_path_structure(self, path_structure):
        """
            Check if the path structure is valid. A path structure is valid if
            it is a tuple containing only the strings 'style', 'tms', 'z', 'x',
            and 'y'. And 'x', 'y', and 'z' are the last three elements in the
            tuple (in any order).

            Parameters
            ----------
            path_structure : tuple
                A tuple of strings that represent the directory structure of a
                tile path.

            Returns
            -------
            (isValid, error_message, error_type) : tuple
                Returns a tuple with a boolean indicating whether the path
                structure is valid, a string with a human-readable error
                message, and an error type that may be raised.
        """
        error_message = ''
        error_type = None
        isValid = False

        required_keys = ['style', 'tms', 'z', 'x', 'y']
        last_3_keys = ['x', 'y', 'z']

        if not isinstance(path_structure, (list, tuple)):
            error_message = 'Path structure must be a tuple or a list'
            error_type = TypeError
        elif not all(isinstance(x, str) for x in path_structure):
            error_message = 'Path structure must be a tuple of strings'
            error_type = TypeError
        elif not all(x in required_keys for x in path_structure):
            error_message = 'Path structure must contain the following ' + \
                f'keys: {required_keys}'
            error_type = ValueError
        elif [x for x in path_structure if x not in required_keys]:
            error_message = 'Path structure contains invalid keys.'
            error_type = ValueError
        elif any(path_structure.count(x) > 1 for x in path_structure):
            error_message = 'Path structure contains duplicate keys'
            error_type = ValueError
        elif not all(x in path_structure[-3:] for x in last_3_keys):
            error_message = 'The last three keys must be: ' + \
                f'{last_3_keys}'
            error_type = ValueError
        else:
            isValid = True
        return (isValid, error_message, error_type)

    def add_base_dir(self, name, dir_path, ext):
        """
            Add a base directory to the Path Manager.

            Parameters
            ----------
            name : str
                The name of the base directory.
            dir_path : str
                The directory path.
            ext : str
                The extension for the files in this directory.
        """

        if not isinstance(name, str):
            raise TypeError('Name must be a string')

        if not isinstance(dir_path, str):
            raise TypeError('Path must be a string')

        if name in self.base_dirs:
            raise ValueError('Base directory name already exists')

        self.base_dirs[name] = {
            'path': dir_path,
            'ext': ext
        }

    def get_base_dir(self, name):
        """
            Get the base directory for the given name.

            Parameters
            ----------
            name : str
                The name of the base directory.

            Returns
            -------
            dict
                A dictionary of the base directory: {'path': '/path/to/dir',
                'ext': '.tif'}
        """
        try:
            return self.base_dirs[name]
        except KeyError:
            return None

    def get_filenames_from_dir(self, base_dir, z=None):
        """
            Given the name of a base directory that has been saved to this Path
            Manager, get all of the files within that directory that have the
            associated extension.

            Parameters
            ----------
            base_dir : str
                The name of the base directory.

            z : int, optional
                The zoom level to get the files for. If not provided, files
                from all zoom levels will be returned.

            Returns
            -------
            paths : list
                A list of paths to the files in the base directory.
        """
        base_dir = self.get_base_dir(base_dir)
        if base_dir is None:
            raise ValueError('Base directory does not exist')

        dir_path = base_dir['path']
        ext = base_dir['ext']

        paths = []
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                if file.endswith(ext):
                    fullpath = os.path.join(root, file)
                    if z is not None:
                        tile = self.dict_from_path(fullpath)
                        if tile['z'] == z:
                            paths.append(fullpath)
                    else:
                        paths.append(fullpath)

        return paths

    def tile_from_path(self, path):
        """
            Create a Tile object from a path that follows the standardized
            directory structure set in this Path Manager.

            Parameters
            ----------
            path : str
                Path to a tile file

            Returns
            -------
            tile : morecantile.Tile
                Tile object
        """
        if not isinstance(path, str):
            raise TypeError('Path must be a string')
        path_dict = self.dict_from_path(path)
        return morecantile.Tile(path_dict['x'], path_dict['y'], path_dict['z'])

    @staticmethod
    def tile_from_string(tile_str=None):
        """
            Parse a morecantile tile that has been cast as a string in the
            format 'Tile(x=6, y=10, z=4)'. Get the x, y, and z values and
            return a morecantile.Tile object.

            Parameters
            ----------
            tile_str : str
                A string in the format 'Tile(x=6, y=10, z=4)' as used by the
                morecantile library.

            Returns
            -------
            tile : morecantile.Tile
                A morecantile.Tile object
        """
        if not isinstance(tile_str, str):
            raise TypeError('Tile must be a string')
        regex = re.compile(r'(?<=x=)\d+|(?<=y=)\d+|(?<=z=)\d+')
        x, y, z = [int(i) for i in regex.findall(tile_str)]
        tile = morecantile.Tile(x, y, z)
        return tile

    def tile_from_path_or_str(self, path_or_str):
        """
            Create a Tile object from a path that follows the standardized
            directory structure set in this Path Manager or a or a serialized
            tile string in the format "Tile(x=6, y=10, z=4)". The method will
            try to determine if the input string is a path or a tile, and then
            parse as needed.

            Parameters
            ----------
            path_or_str : str
                Path to a tile file or seralized tile string

            Returns
            -------
            tile : morecantile.Tile
                Tile object
        """

        if not isinstance(path_or_str, str):
            raise TypeError('Path or tile string must be a string')

        if path_or_str.startswith('Tile(') and path_or_str.endswith(')'):
            return self.tile_from_str(path_or_str)
        else:
            return self.tile_from_path(path_or_str)

    def dict_from_path(self, path):
        """
            Given a path that follows the structure specified in the init
            method, return a dictionary of the path components.

            Parameters
            ----------
            path : str
                Path to a tile file

            Returns
            -------
            dict
                A dictionary of the path components: {'style': '', 'tms': '',
                'z': '', 'x': '', 'y': '', 'ext'='', 'base_dir': ''}
        """

        path_dict = {}
        path_items = path.split(os.sep)

        # Remove the extension from the last item in the path, save it in the
        # path_dict
        ext = os.path.splitext(path_items[-1])[1]
        path_dict['ext'] = ext
        path_items[-1] = path_items[-1].replace(ext, '')

        # Assume the last three items are the x, y, and z indices, in the order
        # specified in the path_structure
        index_values = path_items[-3:]
        index_keys = self.path_structure[-3:]
        for i in range(len(index_keys)):
            path_dict[index_keys[i]] = int(index_values[i])

        path_items = path_items[:-3]

        if(len(path_items) == 1):
            # Style and base_dir are optional, so if they are not present,
            # assume the remaining item in the path_items list is the tms_id
            path_dict['tms'] = path_items[0]
        else:
            # The last two items in the path_items list are the style and tms,
            # in the order specified in the path_structure.
            other_values = path_items[-2:]
            other_keys = self.path_structure[:2]
            for i in range(len(other_keys)):
                path_dict[other_keys[i]] = other_values[i]
            path_items = path_items[:-2]

            # Any remaining items in the path_items list are the base
            # directory. Add the base_dir to the path_dict as a combined string
            # (if multiple parts)
            if len(path_items) > 0:
                path_dict['base_dir'] = os.path.join(*path_items)

        return path_dict

    def path_from_dict(self, dict):
        """
            Given a dictionary of ALL path components, return the path to the
            tile file.

            Parameters
            ----------
            dict : dict
                A dictionary of all the path components: {'style': '', 'tms':
                '', 'z': '', 'x': '', 'y': '', 'base_dir': '', 'ext': ''}

            Returns
            -------
            str
                Path to a tile file
        """

        path_elements_in_order = [dict['base_dir']]
        for key in self.path_structure:
            if key != 'ext':
                path_el = str(dict[key])
                path_elements_in_order.append(path_el)
        # Add extension to last element in the list
        path_elements_in_order[-1] = path_elements_in_order[-1] + dict['ext']

        return os.path.join(*path_elements_in_order)

    def path_from_tile(self, tile=None, base_dir=None, style=None):
        """
            Given a tile, create a path that follows the standardized directory
            structure set in this Path Manager.

            Parameters
            -----------
            tile : morecantile.Tile or str, required
                Tile object to create a path for. A string in the format
                'Tile(x=6, y=10, z=4)' or a standard tile path (e.g. one that
                has a different base dir/ext) is also accepted.
            base_dir : str or dict, required
                The name of the base directory to use, or a dictionary of the
                base directory: {'path': '/path/to/dir', 'ext': '.tif'}.
                Provide a base directory even if no base directory needs to be
                prepended to the path, so that the file extension can be added.
            style : str, optional
                The style, layer name, or statistic that the tile belongs to.
                Defaults to None.

            Returns
            -------
            path : str
                A string with a path to the tile file.
        """

        if isinstance(tile, str):
            tile = self.tile_from_path_or_str(tile)

        if (tile is None):
            raise TypeError('A tile is required to create a path')

        if (base_dir is None):
            raise TypeError('A base directory is required to create a path')

        if (isinstance(base_dir, str)):
            base_dir = self.get_base_dir(base_dir)
            if base_dir is None:
                raise ValueError('Base directory does not exist')

        if (style is None):
            style = ''

        path_dict = {
            'style': style,
            'tms': self.tms_id,
            'z': tile.z,
            'x': tile.x,
            'y': tile.y,
            'ext': base_dir['ext'],
            'base_dir': base_dir['path']
        }

        return self.path_from_dict(path_dict)

    def get_child_tiles(self, tile):
        """
            Get the 1-zoom-level-down children tiles that a given tile
            comprises. Assumes every parent tile comprises 4 children tiles.

            Parameters
            ----------
            tile : morecantile.Tile or str, required
                A tile object from the morecantile library. A string in the
                format 'Tile(x=6, y=10, z=4)' or a standardized path to a tile
                is also accepted.

            Returns
            -------
            list of morecantile.Tile
                A list of morecantile.Tile objects.
        """

        if isinstance(tile, str):
            tile = self.tile_from_path_or_str(tile)

        if (tile is None):
            raise TypeError('A tile is required to get child tiles')

        x2 = tile.x * 2
        y2 = tile.y * 2

        child_z = tile.z + 1
        child_x = (x2, (x2 + 1))
        child_y = (y2, (y2 + 1))

        tiles = []

        for x in child_x:
            for y in child_y:
                tile = morecantile.Tile(x, y, child_z)
                tiles.append(tile)

        return tiles

    def get_parent_tile(self, tile):
        """
            Get the 1-zoom-level-up parent tile that includes a given child
            tile.

            Parameters
            ----------
            tile : morecantile.Tile or str, required
                A tile object from the morecantile library. A string in the
                format 'Tile(x=6, y=10, z=4)' or a path to a standardized tile
                is also accepted.

            Returns
            -------
            morecantile.Tile
                The parent tile
        """

        if isinstance(tile, str):
            tile = self.tile_from_path_or_str(tile)

        if tile is None:
            raise TypeError('A tile is required to get the parent tile')

        parent_z = tile.z - 1
        parent_x = tile.x // 2
        parent_y = tile.y // 2

        parent_tile = morecantile.Tile(parent_x, parent_y, parent_z)

        return parent_tile

    def get_child_paths(self, tile, base_dir=None, style=None):
        """
            For a given tile, get the paths to the tiles that are the children
            of the tile.

            Parameters
            ----------
            tile : morecantile.Tile or str
                The parent tile to get the children tile paths for. A string in
                the format 'Tile(x=6, y=10, z=4)' or a path to a standardized
                tile is also accepted.
            base_dir : str or dict, required
                The name of the base directory where the child tiles are
                stored, or a dictionary of the base directory: {'path':
                '/path/to/dir', 'ext': '.tif'}. Provide a base directory even
                if no base directory needs to be prepended to the path, so that
                the file extension can be added.
            style : str, optional
                The style, layer name, or statistic that the children tiles
                belongs to. Defaults to None.

            Returns
            -------
            list of str
                A list of paths to the child tiles.
        """

        if isinstance(tile, str):
            tile = self.tile_from_path_or_str(tile)

        if (tile is None):
            raise TypeError('A tile is required to create a path')

        child_tiles = self.get_child_tiles(tile)
        paths = []
        for t in child_tiles:
            paths.append(self.path_from_tile(t, base_dir, style))
        return paths

    def get_parent_path(self, tile, base_dir=None, style=None):
        """
            For a given tile, get the path to the parent tile.

            Parameters
            ----------
            tile : morecantile.Tile or str
                The child tile to get the parent tile path for. A string in the
                format 'Tile(x=6, y=10, z=4)' or a path to a standardized tile
                is also accepted.
            base_dir : str or dict, required
                The name of the base directory where the parent tile is stored,
                or a dictionary of the base directory: {'path': '/path/to/dir',
                'ext': '.tif'}. Provide a base directory even if no base
                directory needs to be prepended to the path, so that the file
                extension can be added.
            style : str, optional
                The style, layer name, or statistic that the parent tile
                belongs to. Defaults to None.

            Returns
            -------
            path : str
                A string with a path to the parent tile file.
        """
        if isinstance(tile, str):
            tile = self.tile_from_path_or_str(tile)

        if tile is None:
            raise ValueError('A tile is required to get the parent tile path')

        parent_tile = self.get_parent_tile(tile)
        return self.path_from_tile(parent_tile, base_dir, style)

    @staticmethod
    def create_dirs(paths):
        """
            For a given path or paths, create the directories that the path
            points to if they do not already exist.

            Parameters
            ----------
            paths : str or list, required
                A string or list of strings with paths to directories.

            Returns
            -------
            None
        """

        if isinstance(paths, str):
            paths = [paths]

        for path in paths:
            dir = os.path.dirname(path)
            if not os.path.exists(dir):
                os.makedirs(dir)

    @staticmethod
    def remove_nonexistent_paths(paths):
        """
            Helper method that takes a path or list of paths, and removes any
            for which files do not exist.

            Parameters
            ----------
            paths : str or list, required
                A list of strings with paths to files.

            Returns
            -------
            list or str or None
                If a list of paths was passed, then the list of paths that
                exist. If a string was passed, then the path that exists. If no
                paths exist, then None is returned.
        """

        if isinstance(paths, str):
            paths = [paths]

        if isinstance(paths, list):
            paths = [path for path in paths if os.path.exists(path)]
            if paths:
                return paths
            else:
                return None
        else:
            return None

    @staticmethod
    def tile(x, y, z):
        """
            Helper method that takes x, y, and z indices and returns a
            morecantile.Tile object.
        """
        return morecantile.Tile(x, y, z)

    def get_bounding_box(self, tile=None, as_dict=True):
        """
            Get the bounding box for a given tile.

            Parameters
            ----------
            tile : morecantile.Tile or str
                The tile to get the bounding box for. A string in the format
                'Tile(x=6, y=10, z=4)' is also accepted, as is a path to a tile
                file that follows the standardized path structure set in this
                class.
            as_dict : bool, optional
                If True, the bounding box is returned as a dictionary with keys
                'top', 'left', 'bottom', and 'right'. If False, the bounding
                box is returned as a morecantile.BoundingBox object. Defaults
                to True.

            Returns
            -------
            dict or morecantile.BoundingBox
                The bounding box for the tile, in the unit of the CRS of the
                TMS set on this class.
        """

        if isinstance(tile, str):
            tile = self.tile_from_path_or_str(tile)

        if tile is None:
            raise ValueError('A tile is required to get the bounding box')

        bounding_box = self.tms.bounds(tile)

        if(as_dict):
            bounding_box = {'top': bounding_box.top,
                            'left': bounding_box.left,
                            'bottom': bounding_box.bottom,
                            'right': bounding_box.right}

        return bounding_box

    def get_total_bounding_box(self, dir, z=None):
        """
            Get the total bounding box for all tiles in a directory.

            Parameters
            ----------
            dir : str
                The name of the directory to get the total bounding box for.
            z : int
                The zoom level to get the total bounding box for. If None, the
                total bounding box for all zoom levels is returned.

            Returns
            -------
            list
                The total bounding box for the tiles in the directory, in the
                unit of the CRS of the TMS set on this class.
        """

        tile_paths = self.get_filenames_from_dir(dir, z=z)
        tiles = [self.tile_from_path(path) for path in tile_paths]
        # get the total bounds for all the tiles
        bounds = [self.get_bounding_box(tile) for tile in tiles]
        # get the total bounds for all the tiles
        polygons = [
            box(b['left'], b['bottom'], b['right'], b['top']) for b in bounds]
        bounds_gs = gpd.GeoSeries(polygons, crs=self.tms.crs)
        total_bounds = list(bounds_gs.total_bounds)

        return total_bounds
