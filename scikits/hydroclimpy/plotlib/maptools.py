"""
===================
Maptools extensions
===================

Tools for plotting maps.

Adding parallels and meridians
==============================

.. autofunction:: get_default_parallels
.. autofunction:: get_default_meridians
.. autofunction:: draw_default_parallels
.. autofunction:: draw_default_meridians


Adding a colorbar
=================

.. autofunction:: draw_colorbar


OGR/matplotlib conversions
==========================

.. autofunction:: polygon_to_geometry
.. autofunction:: geometry_to_vertices


:class:`mpl_toolkits.basemap.Basemap` extensions
================================================

.. autofunction:: fillwaterbodies


:author: Pierre GF Gerard-Marchant
:contact: pierregmcode_at_gmail_dot_com
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 257 $"
__date__ = '$Date: 2007-08-26 19:49:21 -0400 (Sun, 26 Aug 2007) $'

import warnings

import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as pyplot
import matplotlib.delaunay as delaunay
from matplotlib.patches import Polygon
import mpl_toolkits.basemap as mtb
from mpl_toolkits.basemap import Basemap

try:
    try:
        import osgeo.ogr as ogr
    except ImportError:
        import ogr
    has_osgeo = True
except ImportError, msg:
    warnmsg = "%s\nUnable to load the OGR/matplotlib interface in %s"
    warnings.warn(warnmsg % (msg, __file__), ImportWarning)
    has_osgeo = False

def get_default_parallels(basemap, resol=1):
    """
    Returns an array of default parallels at the given resolution
    for the current basemap.

    Parameters
    ----------
    basemap : Basemap object
        Basemap on which to draw.
    resol : {float}, optional
        Space between two parallels (in degrees).
    """
    return np.arange(np.floor(basemap.llcrnrlat / 5.) * 5.,
                     np.ceil(basemap.urcrnrlat / 5.) * 5., float(resol))

def get_default_meridians(basemap, resol=1):
    """
    Returns an array of default meridians at the given resolution for the current
    basemap.

    Parameters
    ----------
    basemap : Basemap object
        Basemap on which to draw.
    resol : {float}, optional
        Space between two meridians (in degrees).
    """
    return np.arange(np.floor(basemap.llcrnrlon / 5.) * 5.,
                     np.ceil(basemap.urcrnrlon / 5.) * 5., float(resol))

def draw_default_parallels(basemap, parallels=None, resol=1., ax=None, **kwargs):
    """
    Draws parallels on the given basemap and axis.

    Parameters
    ----------
    basemap : Basemap object
        Basemap on which to draw.
    parallels : {None, sequence}, optional
        List of meridians to draw. 
        If None, default parallels are determined from the current basemap 
        at the given resolution.
    resol : {float}, optional
        Space between two parallels (in degrees).
    ax : {None, Axe}, optional
        Axe on which to draw the parallels.
        If None, the current axe is selected.
    """
    if parallels is None:
        parallels = get_default_parallels(basemap, resol=resol)
    if ax is None:
        ax = pyplot.gca()
    labels = [int(ax.is_first_col()), int(ax.is_last_col()), 0, 0]
    kwargs.update(dict(labels=labels, ax=ax))
    basemap.drawparallels(parallels, **kwargs)
    return
Basemap.draw_default_parallels = draw_default_parallels


def draw_default_meridians(basemap, meridians=None, resol=1., ax=None, **kwargs):
    """
    Draws meridians on the given basemap and axis.

    Parameters
    ----------
    basemap : Basemap object
        Basemap on which to draw.
    meridians : {None, sequence}, optional
        List of meridians to draw. 
        If None, default meridians are determined from the current basemap 
        at the given resolution.
    resol : {float}, optional
        Space between two parallels (in degrees).
    ax : {None, Axe}, optional
        Axe on which to draw the parallels.
        If None, the current axe is selected.
    """
    if meridians is None:
        meridians = get_default_meridians(basemap, resol=resol)
    if ax is None:
        ax = pyplot.gca()
    labels = [0, 0, int(ax.is_first_row()), int(ax.is_last_row())]
    kwargs.update(dict(labels=labels, ax=ax))
    basemap.drawmeridians(meridians, **kwargs)
    return
Basemap.draw_default_meridians = draw_default_meridians



def draw_colorbar(mappable, fig=None, axislist=None, orientation='horizontal',
                  width=0.025, offset=0.05, **kwargs):
    """
    Draws a colorbar for the current mappable object.

    In addition to the parameters listed below, the function accepts all the 
    optional parameters of :func:`matplotlib.pylab.colorbar`.

    Parameters
    ----------
    mappable :
        A mappable object (image, contours...)
    fig : {None, Figure}, optional
        Figure on which to draw the mappable.
        If None, uses the current figure.
    axislist: {None, sequence}, optional
        List of axes to take as reference for drawing.
        If None, uses the current list of axes.
    orientation : {'horizontal', 'horiz', 'h', 'vertical', 'vert', 'v'}, optional
        Orientation of the colorbar.
    width : {0.025, float}, optional
        Width/height of the colorbar, in axes units.
    offset : {0.05, optional}
        Offest from the right/bottom.
    """
    if fig is None:
        fig = pyplot.gcf()
    if axislist is None:
        axislist = fig.axes
    #
    positions = np.array([ax.get_position().get_points() for ax in axislist])
    x0 = positions[:, 0, 0].min()
    x1 = positions[:, 1, 0].max()
    y0 = positions[:, 0, 1].min()
    y1 = positions[:, 1, 1].max()
    #
    if orientation in ('vertical', 'vert', 'v'):
        orientation = 'vertical'
        cbar_ax = pyplot.axes([x1 + offset, y0, width, y1 - y0])
    elif orientation in ('horizontal', 'horiz', 'h'):
        orientation = 'horizontal'
        cbar_ax = pyplot.axes([x0, y0 - offset, x1 - x0, width])
    #
    norm = mappable.norm
    magord = min(np.floor(np.log10(max(np.abs(norm.vmin), 1))),
                 np.floor(np.log10(max(np.abs(norm.vmax), 1))))
    format = kwargs.pop('format', None)
    if format is None:
        format = "%%.%if" % max(1 - magord, 0)
    #
    cbar = pyplot.colorbar(mappable, cbar_ax,
                           format=format, orientation=orientation,
                           **kwargs
                           )
    return cbar


#####--------------------------------------------------------------------------
#---- --- OGR conversion ---
#####--------------------------------------------------------------------------

def polygon_to_geometry(polygon):
    """
    Creates a new ogr.Geometry object from a matplolib.Polygon.
    """
    if not isinstance(polygon, Polygon):
        raise TypeError, "The input data should be a valid Polygon object!"
    listvert = ["%s %s" % (x, y)  for (x, y) in polygon.get_verts()]
    listvert.append(listvert[0])
    return ogr.CreateGeometryFromWkt("POLYGON ((%s))" % ','.join(listvert))

#-----------------------------------------------------------------------------
#---- OGR to matplotlib ---
#-----------------------------------------------------------------------------

def _geometry_to_vertices(geometry):
    """
    Creates a list of vertices (x,y) from the current geometry.
    """
    verts = []
    for nl in range(geometry.GetGeometryCount()):
        line = geometry.GetGeometryRef(nl)
        points = [(line.GetX(i), line.GetY(i)) for i in range(line.GetPointCount())]
        verts.append(points)
    return verts


def geometry_to_vertices(geometry):
    """
    Creates lists of vertices (x,y) from the current geometry.
    """
    if not isinstance(geometry, ogr.Geometry):
        raise TypeError, "The input data should be a valid ogr.Geometry object!"
    vertices = []
    #
    if geometry.GetGeometryType() == ogr.wkbPolygon:
        vertices.extend(_geometry_to_vertices(geometry))
    elif geometry.GetGeometryType() == ogr.wkbMultiPolygon:
        for np in range(geometry.GetGeometryCount()):
            poly = geometry.GetGeometryRef(np)
            vertices.extend(_geometry_to_vertices(poly))
    return vertices

#####--------------------------------------------------------------------------
#---- --- Add-ons
#####--------------------------------------------------------------------------

def fillwaterbodies(basemap, color='blue', inlands=True, ax=None, zorder=None):
    """
    Fills the water bodies with color.
    If inlands is True, inland water bodies are also filled.

    Parameters
    ----------
    basemap : Basemap
        The basemap on which to fill.
    color : {'blue', string, RGBA tuple}, optional
        Filling color
    inlands : {True, False}, optional
        Whether inland water bodies should be filled.
    ax : {None, Axes instance}, optional
        Axe on which to plot. If None, the current axe is selected.
    zorder : {None, integer}, optional
        zorder of the water bodies polygons.

    """
    if not isinstance(basemap, Basemap):
        raise TypeError, "The input basemap should be a valid Basemap object!"
    #
    if ax is None and basemap.ax is None:
        try:
            ax = pyplot.gca()
        except:
            import matplotlib.pyplot as pyplot
            ax = pyplot.gca()
        basemap.ax = ax
    #
    coastpolygons = basemap.coastpolygons
    coastpolygontypes = basemap.coastpolygontypes
    (llx, lly) = (basemap.llcrnrx, basemap.llcrnry)
    (urx, ury) = (basemap.urcrnrx, basemap.urcrnry)
    background = Polygon([(llx, lly), (urx, lly), (urx, ury), (llx, ury)])
    #
    ogr_background = polygon_to_geometry(background)
    inland_polygons = []
    #
    for (poly, polytype) in zip(coastpolygons, coastpolygontypes):
        if polytype != 2:
            verts = ["%s %s" % (x, y) for (x, y) in zip(*poly)]
            ogr_poly = ogr.CreateGeometryFromWkt("POLYGON ((%s))" % ','.join(verts))
            ogr_background = ogr_background.Difference(ogr_poly)
        else:
            inland_polygons.append(Polygon(zip(*poly),
                                           facecolor=color,
                                           edgecolor=color,
                                           linewidth=0))
    #
    background = geometry_to_vertices(ogr_background)
    for xy in background:
        patch = Polygon(xy, facecolor=color, edgecolor=color, linewidth=0)
        if zorder is not None:
            patch.set_zorder(zorder)
        ax.add_patch(patch)
    #        
    if inlands:
        for poly in inland_polygons:
            if zorder is not None:
                poly.set_zorder(zorder)
            ax.add_patch(poly)
    basemap.set_axes_limits(ax=ax)
    return



def extrapolate_data(dataset, basemap, gridsize_x, gridsize_y,
                     maskoceans=False):
    """
    Extrapolate `dataset` on a grid of size `(gridsize_x, gridsize_y)`
    based on `basemap`.

    A regular grid of the user-defined size is created from the basemap.
    The dataset coordinates are then Delaunay triangulated, and the corresponding
    data extrapolated on the regular grid using the nearest-neighbor method

    Parameters
    ----------
    dataset : ndarray
        A structured ndarray, w/ fields ['lon', 'lat', 'data']
    basemap : Basemap
        The projection basemap
    gridsize_x : int
        Number of cells in the x direction ('lon')
    gridsize_y : int
        Number of cells in the x direction ('lat')
    maskoceans : 
    """
    # Get the grid
    (glon, glat, gx, gy) = basemap.makegrid(gridsize_x, gridsize_y,
                                            returnxy=True)
    # Transforms the lon/lat of the dataset in basemap units
    (llon, llat) = basemap(dataset['lon'], dataset['lat'])
    # Triangulate the dataset
    triangul = delaunay.Triangulation(llon, llat)
    # Define an extrapolator (using natural neighbors)...
    # ... and extrapolate the data along the grid...
    extrapolator = triangul.nn_extrapolator(dataset['data'])
    extrapolated = ma.fix_invalid(extrapolator(gx, gy))
    if maskoceans:
        extrapolated = mtb.maskoceans(glon, glat, extrapolated)
    return (extrapolated, gx, gy)

