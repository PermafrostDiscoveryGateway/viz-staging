"""
This module comprises models that represent grids in geographic space.
"""

from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import box
from numpy import linspace, searchsorted, array, logical_and, logical_or
from morecantile import tms as morecantile_tms
from uuid import uuid4
from warnings import warn


class Grid():
    """
    The Grid class represents a lattice of rows and columns, evenly spaced
    along projected or geographic coordinates. Rows and columns are each
    assigned a unique integer index, starting at 0 or some other specified
    integer and increasing consecutively in steps of 1 (i.e. rows that are
    adjacent to each other have index n and n+1).

    A Grid is similar to a 'fishnet' in ArcGIS.
    """

    ROW_IND_NAME = 'grid_row_index'
    COL_IND_NAME = 'grid_column_index'
    ID = 'GRID_' + str(uuid4())

    def __init__(
        self,
        bounds=[-180, -90, 180, 90],
        crs=4326,
        nrows=10,
        ncols=10,
        first_row_i=0,
        first_col_i=0,
        left_to_right=True,
        top_to_bottom=True
    ):
        """
        Initialize a grid.

        Parameters
        ----------
        bounds : list of floats
            The bounds of the grid.
        crs : int
            The coordinate reference system of the grid.
        nrows : int
            The number of rows in the grid.
        ncols : int
            The number of columns in the grid.
        first_row_i : int
            The index of the first row.
        first_col_i : int
            The index of the first column.
        left_to_right : bool
            When True (default) the grid columns are numbered from left to
            right, otherwise they are numbered from right to left.
        top_to_bottom : bool
            When True (default) the grid rows are numbered from top (North) to
            bottom (South), otherwise they are numbered from bottom (South) to
            top (North).
        """
        self.bounds = bounds
        self.nrows = nrows
        self.ncols = ncols
        self.first_row_i = first_row_i
        self.first_col_i = first_col_i
        self.crs = crs
        self.left_to_right = left_to_right
        self.top_to_bottom = top_to_bottom
        self.__build_grid__()

    def __build_grid__(self):
        """
        Build the grid.
        """

        if hasattr(self, '_gdf_cells'):
            del self._gdf_cells

        # Get min and max values for the x and y coordinates.
        self.minx, self.miny, self.maxx, self.maxy = \
            minx, miny, maxx, maxy = self.bounds

        # Calculate height, width, area of the entire grid
        self.height = self.maxy - self.miny
        self.width = self.maxx - self.minx
        self.area = self.height * self.width

        # Calculate height, width, area of each grid cell
        self.cell_height = self.height / self.nrows
        self.cell_width = self.width / self.ncols
        self.cell_area = self.cell_height * self.cell_width

        # Get number of rows and columns
        nrows = self.nrows
        ncols = self.ncols
        self.ncells = nrows * ncols

        # Calculate the location of the x and y grid lines (i.e. fences)
        self.row_fences = rf = linspace(miny, maxy, nrows + 1)
        self.col_fences = cf = linspace(minx, maxx, ncols + 1)

        # Calculate the row and column indices
        first_row_i = self.first_row_i
        first_col_i = self.first_col_i

        last_row_i = first_row_i + nrows
        last_col_i = first_col_i + ncols

        self.row_indices = ri = list(range(first_row_i, last_row_i))
        self.col_indices = ci = list(range(first_col_i, last_col_i))

        if not self.left_to_right:
            self.col_indices = ci = ci[::-1]
        if self.top_to_bottom:
            self.row_indices = ri = ri[::-1]

        # Make row and column geodataframes
        r_geoms = [box(minx, rf[i], maxx, rf[i + 1]) for i in range(nrows)]
        c_geoms = [box(cf[i], miny, cf[i + 1], maxy) for i in range(ncols)]
        self.gdf_rows = GeoDataFrame(
            {self.ROW_IND_NAME: ci, 'geometry': r_geoms}, crs=self.crs)
        self.gdf_cols = GeoDataFrame(
            {self.COL_IND_NAME: ci, 'geometry': c_geoms}, crs=self.crs)
        # TODO ^ set row & col indices as GDF index?

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
        Spatially superimpose the grid onto a GeoDataFrame of polygons, and
        return the GeoDataFrame with new polygon shapes based on places where
        the polygons overlap with the grid. The resulting GeoDataFrame will
        also have a new MultiIndex comprising the grid's row and column
        indices. If the goal is not to create new geometries, then use the
        Grid.sjoin method.

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
            The GeoDataFrame to overlay the grid onto. The GeoDataFrame must be
            in the same coordinate reference system as the grid, and must
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

        self.__all_geom_type__(gdf_c, 'Polygon', True)
        gdf_c = self.__check_crs_match__(gdf_c)

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

    def sjoin(self, gdf, how='left', predicate='intersects'):
        """
        Spatially join the grid and a GeoDataFrame. Return either the
        GeoDataFrame with the associated row & column indices for each
        geometry, or return the grid cell geometries with properties from the
        GeoDataFrame. This method does not alter geometries in anyway, it only
        adds additional properties. If the goal is to create new geometries,
        then use the Grid.overlay method.

        This method uses GeoPandas's sjoin method, but certain circumstances,
        it uses other methods which are sped up by taking advantage of the
        grid's regular and known structure.

        All binary predicates that are supported by GeoPandas.sjoin are
        supported by this sjoin method. If using the 'intersection' method with
        only polygon geometries or only point geometries, then this method will
        be multitudes faster than the GeoPandas version (especially with large
        nrows & ncols)

        For more information about sjoin binary predicates, see: -
        https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.sjoin.html
        - https://geopandas.org/en/stable/docs/user_guide/mergingdata.html

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to join the grid onto. The GeoDataFrame must be in
            the same coordinate reference system as the grid.
        how : 'left', 'right', 'inner'
            The type of join:
              - left: use keys from the GeoDataFrame; retain only GeoDataFrame
                geometry column
              - right: use keys from the Grid; retain only Grid cell geometry
                column
              - inner: use intersection of keys from both; retain only
                GeoDataFrame geometry column. In this case, 'inner' will give
                the same result as 'left'.
        predicate : 'intersects', 'contains', 'within', 'crosses', 'overlaps',
            'touches', 'equals', 'disjoint' The binary predicate to use for the
            sjoin. Default is 'intersects'.

        Returns
        -------
        GeoDataFrame
            The input GeoDataFrame with row and column indices added (when
            how='left' or 'inner'), or the grid cell GeoDataFrame with
            properties from the input GeoDataFrame added (when how='right').
        """

        self.__check_join_type__(how)
        self.__check_predicate__(predicate)

        # If the sjoin falls under certain conditions, then use the faster
        # version of sjoin that takes advantage of the grid's uniform
        # structure.
        if predicate == 'intersects':
            if self.__all_geom_type__(gdf, 'Point'):
                return self.__sjoin_points__(gdf, how)
            if self.__all_geom_type__(gdf, 'Polygon'):
                return self.__sjoin_polygons__(gdf, how)
        # Otherwise, use the built-in geopandas sjoin method with the GDF
        # of grid cells (slowest)
        else:
            return self.__sjoin_cells__(gdf, how, predicate)

    def __sjoin_cells__(self, gdf, how='left', predicate='intersects'):
        """
            Spatially join the grid and a GeoDataFrame. This method is the
            slowest of the sjoin method and is called by the Grid.sjoin method
            when when the predicate is not 'intersects' and when the input
            GeoDataFrame does not contains only Point or only Polygon
            geometries.

            See the Grid.sjoin method for more information.
        """
        self.__check_join_type__(how)
        # works but is the slowest method.
        cells = self.gdf_cells
        lsuffix = self.ID + '_left'
        rsuffix = self.ID + '_right'
        joined = gdf.sjoin(
            cells,
            how=how,
            predicate=predicate,
            lsuffix=lsuffix,
            rsuffix=rsuffix)
        # drop all columns starting with the GRID ID
        joined.drop(columns=joined.filter(regex=self.ID), inplace=True)
        return joined

    def __sjoin_polygons__(self, gdf, how='left'):
        """
            Spatially join the grid and a GeoDataFrame. This method is called
            by the Grid.sjoin method when the input GeoDataFrame contains only
            Polygon geometries, and the predicate is 'intersects'.

            See the Grid.sjoin method for more information.
        """

        # TODO: adapt to different predicates (other than intersects)

        self.__check_join_type__(how)

        # works for intersection only so far.
        gdf_c = gdf.copy()
        gdf_c = self.__check_crs_match__(gdf_c)

        # Names for the row and column index columns.
        ri = self.ROW_IND_NAME
        ci = self.COL_IND_NAME
        idc = self.ID

        # give a temporary ID to each geometry
        gdf_c[idc] = range(len(gdf_c))

        # Slice then merge
        index_id_map = self.overlay(gdf_c, 'intersection') \
            .groupby([ri, ci, idc], as_index=False) \
            .size() \
            .drop(columns='size')

        polygons_with_indices = gdf_c \
            .merge(index_id_map, on=idc, how='left') \
            .drop(columns=idc) \
            .reset_index(drop=True)

        if how == 'left' or how == 'inner':
            return polygons_with_indices

        # if how is 'right'
        cells = self.gdf_cells.copy()
        poly_info = polygons_with_indices.drop(columns='geometry')
        poly_info = DataFrame(poly_info)
        cells_with_poly_info = cells \
            .merge(poly_info, on=[ri, ci], how='left')
        return cells_with_poly_info

    def __sjoin_points__(self, gdf, how='left'):
        """
            Spatially join the grid and a GeoDataFrame. This method is called
            by the Grid.sjoin method when when the input GeoDataFrame contains
            only Point geometries, and the predicate is 'intersects'.

            See the Grid.sjoin method for more information.
        """

        # TODO: adapt to different predicates (other than intersects)

        self.__check_join_type__(how)

        gdf_c = gdf.copy()

        # Names for the row and column index columns.
        ri = self.ROW_IND_NAME
        ci = self.COL_IND_NAME
        idc = self.ID

        # make a temporary ID number
        # gdf_c[idc] = range(len(gdf_c))
        x = array(gdf_c.geometry.x)
        y = array(gdf_c.geometry.y)

        # Add an ID column to the GeoDataFrame.
        gdf_c[idc] = ids = array(range(len(gdf_c)))

        row_ind = searchsorted(self.row_fences, y) - 1
        col_ind = searchsorted(self.col_fences, x) - 1

        # Don't count points outside the grid.
        inside_rows = ~logical_or(row_ind < 0, row_ind > self.nrows)
        inside_cols = ~logical_or(col_ind < 0, col_ind > self.ncols)
        inside_grid = logical_and(inside_rows, inside_cols)
        row_ind = row_ind[inside_grid]
        col_ind = col_ind[inside_grid]
        ids = ids[inside_grid]

        # Convert from 0 to nrows and 0 to ncols to first_row_i to first_row_i
        # + nrows and first_col_i to first_col_i + ncols
        col_ind = [self.col_indices[i] for i in col_ind]
        row_ind = [self.row_indices[i] for i in row_ind]

        index_id_map = DataFrame({ri: row_ind, ci: col_ind, idc: ids})

        points_with_indices = gdf_c \
            .merge(index_id_map, on=idc, how='left') \
            .drop(columns=idc)

        if how == 'left' or how == 'inner':
            return points_with_indices

        # if how is 'right'
        cells = self.gdf_cells.copy()
        point_info = points_with_indices.drop(columns='geometry')
        point_info = DataFrame(point_info)
        cells_with_point_info = cells \
            .merge(point_info, on=[ri, ci], how='left')
        return cells_with_point_info

    def __check_crs_match__(self, gdf):
        """
            Check that a GeoDataFrame has a CRS and that it is the same CRS as
            the grid. If it's not, reproject it with a warning. Return the GDF.
        """
        if self.crs != gdf.crs:
            if gdf.crs is None:
                raise ValueError(
                    'The GeoDataFrame requires a coordinate reference system.')
            warn('The CRS of the GeoDataFrame does not match the CRS of the '
                 'grid. The resulting GeoDataFrame will be re-projected to the'
                 ' CRS of the grid.')
            gdf.to_crs(self.crs, inplace=True)
        return gdf

    def __check_join_type__(self, how):
        """
            Check if a spatial join type is valid. Raise an error if it is not.
            Used for sjoin.
        """
        supported_types = ['left', 'right', 'inner']
        if how not in supported_types:
            raise ValueError(
                f'The join type "{how}" is not allowed for a spatial join. '
                'Allowed join types include: {supported_types}')

    def __check_predicate__(self, predicate):
        """
            Check if a spatial predicate is valid. Raise an error if it is not.
            Used for sjoin.
        """
        supported_predicates = ['intersects', 'contains', 'within',
                                'touches', 'crosses', 'overlaps']
        if predicate not in supported_predicates:
            raise ValueError(
                f'The predicate "{predicate}" is not supported by the grid '
                'sjoin method. The supported predicates are: '
                f'{supported_predicates}')

    @staticmethod
    def __all_geom_type__(gdf, geom_type='Polygon', raise_error=False):
        """
            Check if all geometries in a GeoDataFrame are of the same type.
            Optionally raise an error if they are not.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to check.
            geom_type : str, optional
                The geometry type to check for. Default is 'Polygon'.
            raise_error : bool, optional
                If True, raise an error if the geometries are not of the same
                type. Default is False.
        """
        are_all_geom_type = (gdf.geom_type == geom_type).all()
        if raise_error and not are_all_geom_type:
            raise ValueError(
                f'The GeoDataFrame must contain only {geom_type} geometries.')
        return are_all_geom_type

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
