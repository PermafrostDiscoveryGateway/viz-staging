"""
This module comprises models that represent grids in geographic space.
"""

from geopandas import GeoDataFrame
from shapely.geometry import box
from numpy import linspace
from morecantile import tms as morecantile_tms
import matplotlib.pyplot as plt


class Grid():
    """
    The Grid class represents a lattice of rows and columns, even spaced along
    projected or geographic coordinates. Each row and column has a unique
    index, and each cell is a square.
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

    def show(self):
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
        plt.show()


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
