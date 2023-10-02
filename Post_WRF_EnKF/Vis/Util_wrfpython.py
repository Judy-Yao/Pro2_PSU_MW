import numpy as np

#@set_interp_metadata("2dxy")
#def interp2dxy(field3d, xy, meta=True):
#    """Return a cross section for a three-dimensional field.
#
#    The returned array will hold the vertical cross section data along the
#    line described by *xy*.
#
#    This method differs from :meth:`wrf.vertcross` in that it will return
#    all vertical levels found in *field3d*.  :meth:`wrf.vertcross` includes
#    an additional interpolation to set the output to a fixed number of
#    vertical levels.  Also, a :class:`numpy.ma.MaskedArray` is not created
#    and this routine should be considered as low-level access to the underlying
#    Fortran routine.
#
#    See Also:
#
#        :meth:`wrf.xy`, :meth:`wrf.vertcross`
#
#    Args:
#
#        field3d (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The
#            array to interpolate with at least three dimensions, whose
#            rightmost dimensions are nz x ny x nx.
#
#        xy (:class:`xarray.DataArray` or :class:`numpy.ndarray`): An array
#            of one less dimension than *field3d*, whose rightmost dimensions
#            are nxy x 2. This array holds the x,y pairs of a line across the
#            model domain. The requested vertical cross section will be
#            extracted from *field3d* along this line.
#
#        meta (:obj:`bool`, optional): Set to False to disable metadata and
#            return :class:`numpy.ndarray` instead of
#            :class:`xarray.DataArray`.  Default is True.
#
#    Warning:
#
#        The input arrays must not contain any missing/fill values or
#        :data:`numpy.nan` values.
#
#
#    Returns:
#
#        :class:`xarray.DataArray` or :class:`numpy.ndarray`: An array
#        containing the vertical cross section along the line *xy*.  The
#        returned dimensions will be the same as *xy*, but with the rightmost
#        dimensions being nz x nxy. If xarray is enabled and the *meta*
#        parameter is True, then the result will be a :class:`xarray.DataArray`
#        object.  Otherwise, the result will
#        be a :class:`numpy.ndarray` object with no metadata.
#
#    Examples:
#        Example 1:  Calculate the vertical cross section for RH for a diagonal
#        line from the lower left to the upper right of the domain.
#
#        .. code-block:: python
#
#            from wrf import getvar, xy, interp2dxy
#            from netCDF4 import Dataset
#
#            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
#
#            rh = getvar(wrfnc, "rh")
#            start = (0, 0)
#            end = (-1, -1)
#            xy_line = xy(rh, start_point=start, end_point=end)
#
#            vert_cross = interp2dxy(rh, xy_line)
#
#    """
#    return _interp2dxy(field3d, xy)


#@set_interp_metadata("xy")
def xy_edit(field, pivot_point=None, angle=None, start_point=None, end_point=None):
    """Return the x,y points for a line within a two-dimensional grid.

    This function is primarily used to obtain the x,y points when making a
    cross section.

    Args:

        field (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A
            field with at least two dimensions.

        pivot_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries,
            in the form of [x, y] (or [west_east, south_north]), which
            indicates the x,y location through which the plane will pass.
            Must also specify `angle`.

        angle (:obj:`float`, optional): Only valid for cross sections where
            a plane will be plotted through
            a given point on the model domain. 0.0 represents a S-N cross
            section.  90.0 is a W-E cross section.

        start_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries, in the form of
            [x, y] (or [west_east, south_north]), which indicates the start
            x,y location through which the plane will pass.

        end_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries, in the form of
            [x, y] (or [west_east, south_north]), which indicates the end x,y
            location through which the plane will pass.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: An array of
        x,y points, which has shape num_points x 2.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    Examples:
        Example 1: Using Pivot Point and Angle

        .. code-block:: python

            from wrf import getvar, xy
            from netCDF4 import Dataset

            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            field = wrf.getvar(wrfnc, "slp")

            # Use the center of the grid
            pivot = (field.shape[-1]/2.0, field.shape[-2]/2.0)

            # West-East
            angle = 90.0

            xy_points = xy(field, pivot_point=pivot, angle=angle)

        Example 2: Using Start Point and End Point

        .. code-block:: python

            from wrf import getvar, xy
            from netCDF4 import Dataset

            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            field = wrf.getvar(wrfnc, "slp")

            # Make a diagonal of lower left to upper right
            start = (0, 0)
            end = (-1, -1)

            xy_points = xy(field, start_point=start, end_point=end)


    """
    return get_xy_edit(field, pivot_point, angle, start_point, end_point)


def get_xy_edit(var, pivot_point=None, angle=None,
           start_point=None, end_point=None):
    """Return the x,y points for the horizontal cross section line.

    Args:

        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A variable
            that contains a :attr:`shape` attribute.

        pivot_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries,
            in the form of [x, y] (or [west_east, south_north]), which
            indicates the x,y location through which the plane will pass.
            Must also specify `angle`.

        angle (:obj:`float`, optional): Only valid for cross sections where
            a plane will be plotted through
            a given point on the model domain. 0.0 represents a S-N cross
            section.  90.0 is a W-E cross section.

        start_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries, in the form of
            [x, y] (or [west_east, south_north]), which indicates the start
            x,y location through which the plane will pass.

        end_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries, in the form of
            [x, y] (or [west_east, south_north]), which indicates the end x,y
            location through which the plane will pass.

    Returns:

        :class:`np.ndarray`: A two-dimensional array with the left index
        representing each point along the line, and the rightmost dimension
        having two values for the x and y coordinates [0=X, 1=Y].

    """
    if pivot_point is not None:
        pos_pivot = to_positive_idxs(var.shape[-2:], pivot_point)
    else:
        pos_pivot = pivot_point

    if start_point is not None:
        pos_start = to_positive_idxs(var.shape[-2:], start_point)
    else:
        pos_start = start_point

    if end_point is not None:
        pos_end = to_positive_idxs(var.shape[-2:], end_point)
    else:
        pos_end = start_point

    xdim = var.shape[-1]
    ydim = var.shape[-2]

    xy = _calc_xy_edit(xdim, ydim, pos_pivot, angle, pos_start, pos_end)

    return xy

def _calc_xy_edit(xdim, ydim, pivot_point=None, angle=None,
             start_point=None, end_point=None):
    """Return the x,y points for the horizontal cross section line.

    Args:

        xdim (:obj:`int`): The x-dimension size.

        ydim (:obj:`int`): The y-dimension size.

        pivot_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries,
            in the form of [x, y] (or [west_east, south_north]), which
            indicates the x,y location through which the plane will pass.
            Must also specify `angle`.

        angle (:obj:`float`, optional): Only valid for cross sections where
            a plane will be plotted through
            a given point on the model domain. 0.0 represents a S-N cross
            section.  90.0 is a W-E cross section.

        start_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries, in the form of
            [x, y] (or [west_east, south_north]), which indicates the start
            x,y location through which the plane will pass.

        end_point (:obj:`tuple` or :obj:`list`, optional): A
            :obj:`tuple` or :obj:`list` with two entries, in the form of
            [x, y] (or [west_east, south_north]), which indicates the end x,y
            location through which the plane will pass.

    Returns:

        :class:`np.ndarray`: A two-dimensional array with the left index
        representing each point along the line, and the rightmost dimension
        having two values for the x and y coordinates [0=X, 1=Y].

    """
    # Have a pivot point with an angle to find cross section
    if pivot_point is not None and angle is not None:
        xp = pivot_point[-2]
        yp = pivot_point[-1]

        if xp >= xdim or yp >= ydim:
            raise ValueError("pivot point {} is outside of domain "
                             "with shape {}".format(pivot_point,
                                                    (xdim, ydim)))

        if (angle > 315.0 or angle < 45.0
                or ((angle > 135.0) and (angle < 225.0))):

            slope = -(360.-angle)/45.
            if(angle < 45.):
                slope = angle/45.
            if(angle > 135.):
                slope = (angle-180.)/45.

            intercept = xp - yp*slope

            # find intersections with domain boundaries
            y0 = 0.
            x0 = y0*slope + intercept

            if(x0 < 0.):  # intersect outside of left boundary
                x0 = 0.
                y0 = (x0 - intercept)/slope
            if(x0 > xdim-1):  # intersect outside of right boundary
                x0 = xdim-1
                y0 = (x0 - intercept)/slope
            y1 = ydim-1.  # need to make sure this will be a float?
            x1 = y1*slope + intercept

            if(x1 < 0.):  # intersect outside of left boundary
                x1 = 0.
                y1 = (x1 - intercept)/slope

            if(x1 > xdim-1):  # intersect outside of right boundary
                x1 = xdim-1
                y1 = (x1 - intercept)/slope
        else:
            #  y = x*slope + intercept
            slope = (90.-angle)/45.
            if (angle > 225.):
                slope = (270.-angle)/45.
            intercept = yp - xp*slope

            # Find intersections with domain boundaries
            x0 = 0.
            y0 = x0*slope + intercept

            if (y0 < 0.):  # intersect outside of bottom boundary
                y0 = 0.
                x0 = (y0 - intercept)/slope

            if (y0 > ydim-1):  # intersect outside of top boundary
                y0 = ydim-1
                x0 = (y0 - intercept)/slope

            x1 = xdim-1.  # need to make sure this will be a float?
            y1 = x1*slope + intercept

            if (y1 < 0.):  # intersect outside of bottom boundary
                y1 = 0.
                x1 = (y1 - intercept)/slope

            if (y1 > ydim-1):  # intersect outside of top boundary
                y1 = ydim-1
                x1 = (y1 - intercept)/slope
    elif start_point is not None and end_point is not None:
        x0 = start_point[-2]
        y0 = start_point[-1]
        x1 = end_point[-2]
        y1 = end_point[-1]

        if x0 >= xdim or y0 >= ydim:
            raise ValueError("start_point {} is outside of domain "
                             "with shape {}".format(start_point, (xdim, ydim)))

        if x1 >= xdim or y1 >= ydim:
            raise ValueError("end_point {} is outside of domain "
                             "with shape {}".format(end_point, (xdim, ydim)))
    else:
        raise ValueError("invalid start/end or pivot/angle arguments")

    dx = x1 - x0
    dy = y1 - y0
    distance = (dx*dx + dy*dy)**0.5
    npts = int(distance) + 1

    xy = np.zeros((npts, 2), "float")

    dx = dx/(npts-1)
    dy = dy/(npts-1)

    for i in range(npts):
        xy[i, 0] = x0 + i*dx
        xy[i, 1] = y0 + i*dy

    return xy

def to_positive_idxs(shape, coord):
    """Return the positive index values.

    This function converts negative index values to positive index values.

    Args:

        shape (indexable sequence): The array shape.

        coord (indexable sequence): The coordinate pair for x and y.

    Returns:

        :obj:`list`: The coordinate values with all positive indexes.

    """
    if (coord[-2] >= 0 and coord[-1] >= 0):
        return coord

    return [x if (x >= 0) else shape[-i-1]+x for (i, x) in enumerate(coord)]
