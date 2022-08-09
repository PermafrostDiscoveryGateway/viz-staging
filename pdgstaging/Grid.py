"""
This module comprises models that represent grids in geographic space.
"""

from geopandas import GeoDataFrame
from shapely.geometry import box
from numpy import linspace
from morecantile import tms as morecantile_tms
import warnings


class Grid():
    """
    The Grid class represents a lattice of rows and columns, even spaced along
    projected or geographic coordinates. Each row and column has a unique
    index, and each cell is a square. Similar to a 'fishnet' in ArcGIS.
    """

    ROW_IND_NAME = 'grid_row_index'
    COL_IND_NAME = 'grid_column_index'

    def __init__(
        self,
        bounds=[-180, -90, 180, 90],
        nrows=10,
        ncols=10,
        first_row_i=0,
        first_col_i=0,
        crs=4326
    ):
        """
        Initialize a grid.

        Parameters
        ----------
        bounds : list of floats
            The bounds of the grid.
        nrows : int
            The number of rows in the grid.
        ncols : int
            The number of columns in the grid.
        first_row_i : int
            The index of the first row.
        first_col_i : int
            The index of the first column.
        crs : int
            The coordinate reference system of the grid.
        """
        self.bounds = bounds
        self.nrows = nrows
        self.ncols = ncols
        self.first_row_i = first_row_i
        self.first_col_i = first_col_i
        self.crs = crs
        self.__build_grid__()

    def __build_grid__(self):
        """
        Build the grid.
        """

        if hasattr(self, '_gdf_cells'):
            del self._gdf_cells

        self.minx, self.miny, self.maxx, self.maxy = \
            minx, miny, maxx, maxy = self.bounds

        self.height = self.maxy - self.miny
        self.width = self.maxx - self.minx
        self.area = self.height * self.width

        self.cell_height = self.height / self.nrows
        self.cell_width = self.width / self.ncols
        self.cell_area = self.cell_height * self.cell_width

        nrows = self.nrows
        ncols = self.ncols

        # Calculate the number of rows and columns. Reverse the rows so that
        # the first row is the bottom row.
        self.row_fences = rf = linspace(miny, maxy, nrows + 1)[::-1]
        self.col_fences = cf = linspace(minx, maxx, ncols + 1)

        self.row_indices = ri = []
        row_geoms = []
        self.col_indices = ci = []
        col_geoms = []

        # Make nrow row geometries.
        for i in range(nrows):
            ri.append(self.first_row_i + i)
            row_geoms.append(
                box(minx, rf[i], maxx, rf[i + 1]))
        # Make ncol column geometries.
        for i in range(ncols):
            ci.append(self.first_col_i + i)
            col_geoms.append(
                box(cf[i], miny, cf[i + 1], maxy))

        # Make two geodataframes. One for the rows and one for the columns.
        self.gdf_rows = GeoDataFrame(
            {self.ROW_IND_NAME: ri, 'geometry': row_geoms}, crs=self.crs)
        self.gdf_cols = GeoDataFrame(
            {self.COL_IND_NAME: ci, 'geometry': col_geoms}, crs=self.crs)

    @property
    def gdf_cells(self):
        """
        Return a GeoDataFrame with one geometry for each cell in the grid,
        along with the associated row and column indices.
        """
        if not hasattr(self, '_gdf_cells'):
            self._gdf_cells = self.gdf_cols.overlay(
                self.gdf_rows,
                how='union').reset_index(drop=True)
        return self._gdf_cells

    def overlay(self, gdf, how='intersection', **kwargs):
        """
        Spatially superimpose the grid onto a GeoDataFrame of polygons. Return
        the GeoDataFrame with new polygon shapes based on places where the
        polygons overlap with the grid. The resulting GeoDataFrame will also
        have a new MultiIndex comprising the grid's row and column indices.

        This method uses GeoPanda's overlay method, but rather than overlaying
        the grid cells onto the GeoDataFrame, it overlays the rows, then the
        columns. Performing two overlay operations (rows & columns) is faster
        than one (cells), especially with higher row and column counts, but the
        resulting geometries are the same.

        All set-operations except that are supported by GeoPandas.overlay are
        supported by this overlay method, except for 'difference':
          - 'intersection' (the default) will essentially slice the input
            polygons by the grid cells. Wherever a polygon overlaps with a grid
            cell, it will be divided along the grid lines and the resulting
            polygons will be returned.
          - 'union' will return both the sliced polygons from the
            'intersection' operation as well as the grid cells geometries with
            the areas where the polygons overlap removed.
          - 'identity' gives the same result as 'intersection'.
          - 'symmetric difference' will return the grid cell geometries with
            the areas where the polygons overlap removed.

        For more information about spatial set-operations and the overlay
        method, see:
          - https://geopandas.org/en/stable/docs/user_guide/set_operations.html
          - https://geopandas.org/en/stable/docs/reference/api/geopandas.overlay.html
          - https://shapely.readthedocs.io/en/stable/manual.html#set-theoretic-methods

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to overlay the grid onto. The GeoDataFrame must
            be in the same coordinate reference system as the grid, and must
            comprise only polygon geometries.
        how : 'intersection', 'union', 'identity', 'symmetric difference'
            The set-operation to perform. Default is 'intersection'.
        **kwargs : keyword arguments
            Keyword arguments to pass to the GeoPandas.overlay method.

        Returns
        -------
        GeoDataFrame
            The GeoDataFrame with new polygon shapes based on places where the
            polygons overlap with the grid. The resulting GeoDataFrame will
            also have a new MultiIndex comprising the grid's row and column
            indices.
        """

        # Don't modify the original GeoDataFrame.
        gdf_c = gdf.copy()

        # Check validity of inputs
        if how == 'difference':
            raise NotImplementedError(
                'The difference operation is not supported by the grid '
                'overlay method. Use the symmetric_difference operation '
                'instead.')

        supported_methods = [
            'intersection',
            'union',
            'identity',
            'symmetric_difference']
        if how not in supported_methods:
            raise ValueError(
                f'The operation "{how}" is not supported by the grid overlay '
                'method. The supported operations are: {supported_methods}')

        self.__check_all_polygons__(gdf_c)

        if self.crs != gdf_c.crs:
            if gdf_c.crs is None:
                raise ValueError(
                    'The GeoDataFrame requires a coordinate reference system.')
            warnings.warn(
                'The CRS of the GeoDataFrame does not match the CRS of the '
                'grid. The resulting GeoDataFrame will be re-projected to the '
                'CRS of the grid.')
            gdf_c.to_crs(self.crs, inplace=True)

        # Perform the overlay operation.
        if how != 'symmetric_difference':
            gdf_c_cols = gdf_c.overlay(self.gdf_cols, how, **kwargs)
            gdf_c_rows_col = gdf_c_cols.overlay(self.gdf_rows, how, **kwargs)
        else:
            gdf_c_rows_col = gdf_c.overlay(self.gdf_cells, how=how, **kwargs)

        # Set the new MultiIndex.
        gdf_c_rows_col.set_index(
            [self.ROW_IND_NAME, self.COL_IND_NAME], inplace=True)

        return gdf_c_rows_col

    def overlay_cells(self, gdf, how=None, **kwargs):
        gdf_c = gdf.copy()
        if(how != 'difference'):
            gdf_c_rows_col = gdf_c.overlay(self.gdf_cells, how, **kwargs)
            gdf_c_rows_col.set_index(
                [self.ROW_IND_NAME, self.COL_IND_NAME], inplace=True)
        else:
            gdf_c_rows_col = self.gdf_cells.overlay(gdf_c, how, **kwargs)
        return gdf_c_rows_col

    def __check_crs_match__(self, gdf):
        if self.crs != gdf.crs:
            raise ValueError(
                'The grid and the GeoDataFrame do not have the same '
                'coordinate reference system.')

    @staticmethod
    def __check_all_polygons__(gdf):
        are_polys = gdf.geometry.type == 'Polygon'
        if not are_polys.all():
            raise ValueError(
                'The GeoDataFrame must contain only polygons.')

    def plot(self):
        """
        Plot the grid with matplotlib. This is useful for debugging.
        """
        grid_line_color = '#a1a8ab'
        rows = self.gdf_rows
        cols = self.gdf_cols
        ax = rows.plot(edgecolor=grid_line_color, facecolor='none')
        cols.plot(ax=ax, edgecolor=grid_line_color, facecolor='none')
        for i, row in rows.iterrows():
            ax.text(
                row.geometry.bounds[0],
                row.geometry.centroid.y,
                'ROW ' + str(row[self.ROW_IND_NAME]),
                horizontalalignment='left',
                verticalalignment='center')
        for i, col in cols.iterrows():
            ax.text(
                col.geometry.centroid.x,
                col.geometry.bounds[3],
                'COL ' + str(col[self.COL_IND_NAME]),
                horizontalalignment='center',
                verticalalignment='top')
        return ax


class TMSGrid(Grid):
    """
    The TMSGrid class represents a grid that follows one of the OGC
    TileMatrixSets for a specific bounding box and a single z-level.
    Uses morecantile: https://developmentseed.org/morecantile/
    """

    def __init__(self, tms_id, z, bounds):
        """
        Initialize a TMSGrid.

        Parameters
        ----------
        tms_id : str
            The ID of the TileMatrixSet. See morecantile.tms.list() for a
            list of available TileMatrixSets.
        z : int
            The z-level of the grid.
        bounds : list of floats
            The minimum bounds for the grid. All tiles that these bounds
            cover will be included in the grid, and the subsequent grid
            bounds will be expanded to include all tiles.
        """

        self.z = z
        self.tms_id = tms_id
        self.tms = tms = morecantile_tms.get(tms_id)
        self.original_bounds = west, south, east, north = bounds

        LL_EPSILON = 1e-11

        if west > east:
            bbox_west = (tms.bbox.left, south, east, north)
            bbox_east = (west, south, tms.bbox.right, north)
            bboxes = [bbox_west, bbox_east]
        else:
            bboxes = [(west, south, east, north)]

        x = set()
        y = set()

        top_limits = []
        bottom_limits = []
        left_limits = []
        right_limits = []

        for w, s, e, n in bboxes:

            # Clamp bounding values.
            w = max(tms.bbox.left, w)
            s = max(tms.bbox.bottom, s)
            e = min(tms.bbox.right, e)
            n = min(tms.bbox.top, n)

            ul_tile = tms.tile(w + LL_EPSILON, n - LL_EPSILON, z)
            lr_tile = tms.tile(e - LL_EPSILON, s + LL_EPSILON, z)

            ul_tile_bounds = tms.bounds(ul_tile)
            lr_tile_bounds = tms.bounds(lr_tile)

            top_limits.append(ul_tile_bounds.top)
            bottom_limits.append(lr_tile_bounds.bottom)
            left_limits.append(ul_tile_bounds.left)
            right_limits.append(lr_tile_bounds.right)

            bbox_xs = list(range(ul_tile.x, lr_tile.x + 1))
            bbox_ys = list(range(ul_tile.y, lr_tile.y + 1))
            x.update(bbox_xs)
            y.update(bbox_ys)

        top_limit = max(top_limits)
        bottom_limit = min(bottom_limits)
        left_limit = min(left_limits)
        right_limit = max(right_limits)

        super().__init__(
            bounds=[left_limit, bottom_limit, right_limit, top_limit],
            nrows=len(y),
            ncols=len(x),
            first_row_i=min(y),
            first_col_i=min(x),
            crs=tms.crs
        )
