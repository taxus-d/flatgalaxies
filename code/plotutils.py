# -*- coding: utf-8 -*-
"""
    code.plotutils
    ~~~~~~~~~~~~~~~~~~~~~~

    scatter-density plots and utilities
    two backends are supported:
        - plotly (`scatter_density_plotly`)
        - matplotlib (`scatter_density_plot`)

    :copyright: (c) 2019 by taxus-d.
    :license: MIT, see LICENSE for more details.
"""

from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from scipy.spatial import cKDTree
from scipy.interpolate import griddata

import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib as mpl


import plotly.graph_objects as go

from joblib import Memory

cachedir = "./cached/"
memory = Memory(cachedir, verbose=0)


def _colorize_z_none(xx, yy, modepars, plotargs):
    zlabel = "dummy"
    zz = np.zeros_like(xx)
    return zz, zlabel, plotargs


def _colorize_z_hist(xx, yy, modepars, plotargs):
    grid = np.vstack([xx, yy])
    bins = modepars.get('bins', (15, 15))
    h, x_e, y_e = np.histogram2d(xx, yy, bins=bins)

    stepx = x_e[1]-x_e[0]
    x_e = np.pad(x_e, (1, 1), 'constant', constant_values=(x_e[0]-stepx, x_e[-1]+stepx))
    stepy = y_e[1]-y_e[0]
    y_e = np.pad(y_e, (1, 1), 'constant', constant_values=(y_e[0]-stepy, y_e[-1]+stepy))

    # padding helps to avoid extrapolation, as it introduces nans in zz
    h = np.pad(h, ((1, 1), (1, 1)), 'constant', constant_values=((0, 0), (0, 0)))

    xc = (x_e[1:] + x_e[:-1])/2
    yc = (y_e[1:] + y_e[:-1])/2

    zz = interpn((xc, yc), h, grid.T, method = "linear")
    zlabel = "interp hist2d"

    return zz, zlabel, plotargs


def _colorize_z_kde(xx, yy, modepars, plotargs):
    grid = np.vstack([xx, yy])
    kernel = gaussian_kde(grid)
    zz = kernel(grid)
    zlabel = "kde"
    return zz, zlabel, plotargs


def _colorize_z_near(xx, yy, modepars, plotargs):
    scalingx = (np.max(xx) - np.min(xx))
    scalingy = (np.max(yy) - np.min(yy))
    xxn = (xx - np.min(xx)) / scalingx
    yyn = (yy - np.min(yy)) / scalingy

    gridn = np.vstack([xxn, yyn])
    tree = cKDTree(gridn.T)
    searchradius = modepars.get('searchradius', 0.05)
    njobs = modepars.get('njobs', 2)

    nb = tree.query_ball_point(gridn.T, searchradius, n_jobs=njobs)
    zz = np.fromiter(map(len, nb), dtype=int)

    circle = mpath.Path.unit_circle()
    verts = np.copy(circle.vertices)

    # stretching correction
    if 'ax' in plotargs.keys():
        xp1, yp1, xp2, yp2 = plotargs['ax'].figure.bbox.bounds
        ap = (xp2 - xp1) * searchradius
        bp = (yp2 - yp1) * searchradius

        verts[:, 0] *= ap/bp
        verts[:, 1] *= 1
        ell_marker = mpath.Path(verts, circle.codes)

        area_p = ap*bp*np.pi

        if plotargs.get('s', None) in ['fair']:
            plotargs['s'] = area_p
        plotargs['marker'] = plotargs.pop('marker', ell_marker)

    zlabel = r"# neighbours"

    return zz, zlabel, plotargs


colorize_z_type = {
    'none': _colorize_z_none,
    'hist': _colorize_z_hist,
    'kde': _colorize_z_kde,
    'near': _colorize_z_near,
}


def sort_by_zorder(xx, yy, zz, sortp = True):
    idx = range(len(zz))
    if sortp:
        idx = np.argsort(zz)
        xx = np.array(xx)[idx]
        yy = np.array(yy)[idx]
        zz = zz[idx]
    return idx, xx, yy, zz


def scatter_density_plot(
    xx,
    yy,
    xrange,
    yrange,
    xlabel,
    ylabel,
    mode="hist",
    modepars={},
    sort=True,
    contours=False,
    pointlabels=None,
    pointlabellim=np.inf,
    ax=None,
    fig=plt,
    **kwargs
):
    ax = plt.gca() if ax is None else ax

    plotargs = {}
    plotargs['s'] = kwargs.pop('s', mpl.rcParams['lines.markersize']**2)
    plotargs['marker'] = kwargs.pop('marker', 'o')
    plotargs['ax'] = ax

    zz, zlabel, plotargs = colorize_z_type[mode](xx, yy, modepars, plotargs)

    idx, xx, yy, zz = sort_by_zorder(xx, yy, zz, sort)

    plotargs.pop('ax', None)
    mapping = ax.scatter(xx, yy, c=zz, **plotargs, **kwargs)

    if contours:
        grid_x, grid_y = np.mgrid[
            xrange[0]:xrange[1]:100j,
            yrange[0]:yrange[1]:100j]

        grid_z = griddata((xx, yy), zz, (grid_x, grid_y), method="cubic", rescale=True)
        ax.contour(grid_x, grid_y, grid_z, colors="lightgray", linewidths=1)

    pointradius = np.sqrt(plotargs['s']/np.pi)

    if pointlabels is not None:
        pointlabels = np.array(pointlabels)[idx]
        for i in range(len(xx)):
            if zz[i] < pointlabellim:
                ax.annotate(pointlabels[i], (xx[i], yy[i]),
                            (pointradius*1.1, pointradius*1.1),
                            textcoords='offset points', fontsize='6')

    ax.set_xlim(*xrange)
    ax.set_ylim(*yrange)

    ax.set(xlabel=xlabel, ylabel=ylabel)
    fig.colorbar(mapping, ax=ax, extend='max', label=zlabel)
    return ax


def seaborn_to_plotly(scl):
    ''' converts a seaborn color palette to a plotly colorscale '''
    return [[
        float(i)/float(len(scl)-1),
        'rgb' + str((scl[i][0]*255,
                     scl[i][1]*255,
                     scl[i][2]*255))] for i in range(len(scl))]


#  ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
#   'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
#   'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
#   'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
#   'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
#   'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
#   'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
#   'orrd', 'oryel', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg',
#   'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor',
#   'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy',
#   'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral',
#   'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose',
#   'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'twilight',
#   'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd']

def scatter_density_plotly(
    xx,
    yy,
    xrange,
    yrange,
    xlabel,
    ylabel,
    logx=False,
    logy=False,
    mode="hist",
    modepars={},
    subplotpos={},
    subplotlyt={},
    fig=None,
    sort=True,
    contours=False,
    pointlabels=None,
    scattertype=go.Scattergl,
    **kwargs
):
    fig = go.Figure() if fig is None else fig
    plotargs = {}
    zz, zlabel, plotargs = colorize_z_type[mode](xx, yy, modepars, plotargs)
    idx, xx, yy, zz = sort_by_zorder(xx, yy, zz, sort)

    alpha = kwargs.get("alpha", 1)
    points = scattertype(
        x=xx,
        y=yy,
        mode="markers",
        marker=dict(
            color=zz,
            showscale=True,
            colorscale="matter",
            reversescale=True,
            opacity=alpha,
            colorbar={"title": zlabel},
        ),
        showlegend=False,
    )
    subplotindex = (subplotpos.get("row", 1) - 1) * subplotlyt.get("rows", 1)
    subplotindex += subplotpos.get("col", 1)
    axisind = "" if subplotindex == 1 else str(subplotindex)

    xrange, xtype = (
        ((np.log10(xrange[0]), np.log10(xrange[1])), "log")
        if logx
        else (xrange, "linear")
    )
    yrange, ytype = (
        ((np.log10(yrange[0]), np.log10(yrange[1])), "log")
        if logy
        else (yrange, "linear")
    )
    layout = {
        "clickmode": "select",
        "xaxis" + axisind: {"title": {"text": xlabel}, "range": xrange, "type": xtype},
        "yaxis" + axisind: {"title": {"text": ylabel}, "range": yrange, "type": ytype},
    }
    if pointlabels is not None:
        points.hoverinfo = "text"
        points.hovertext = list(np.array(pointlabels)[idx])

    if contours:
        grid_x, grid_y = np.mgrid[
            np.min(xx) : np.max(xx) : 30j,
            np.min(yy) : np.max(yy) : 30j
        ]
        grid_z = griddata((xx, yy), zz, (grid_x, grid_y), method="cubic", rescale=True)
        contour = go.Contour(
            x=grid_x[:, 0],
            y=grid_y[0, :],
            z=grid_z.T,
            contours_coloring="lines",
            showlegend=False,
            showscale=False,
            line={"width": 1, "color": "lightgray"},
            hoverinfo="none",
        )
        fig.add_trace(contour, **subplotpos)

    fig.add_trace(points, **subplotpos)
    fig.update_layout(layout)

    return fig

