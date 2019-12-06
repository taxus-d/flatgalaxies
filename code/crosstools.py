# -*- coding: utf-8 -*-
"""
    code.crosstools
    ~~~~~~~~~~~~~~~

    Various utilities to fetch data with caching

    :copyright: (c) 2019 by taxus-d.
    :license: MIT, see LICENSE for more details.
"""

import numpy as np
from joblib import Memory

from astropy.visualization import PercentileInterval, AsinhStretch, LogStretch, LinearStretch
from scipy.ndimage import rotate as rotim
from matplotlib.patches import Ellipse, Circle
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS

import astropy.utils

import warnings

from .panstarrs import geturl

from argparse import Namespace
PANPARS = Namespace()
PANPARS.scale = 4
PANPARS.defaultvalue = -999.0

astropy.utils.data.Conf.remote_timeout = 100

cachedir = "./cached/"
memory = Memory(cachedir, verbose=0)


def clean_cache():
    memory.clear()


@memory.cache
def cone_search_getobjs(jobs, query, **kwargs):
    results = jobs.quick(query, task_name="galaxy-like cone search")
    df = results.to_pandas()
    df[df == PANPARS.defaultvalue] = np.nan
    return df


def cone_galaxy_search(jobs, template, pos, size, filt):
    """
    pos -- in degrees, size -- in arcmin
    """
    query = template.format(ra=pos[0], dec=pos[1], s=size, f=filt)

    return cone_search_getobjs(jobs, query)


def inellipse(pos, center, theta, a, b):
    """
    A test whether a point lies inside rotated ellipse
    """
    c = np.cos(np.radians(theta))
    s = np.sin(np.radians(theta))
    r = b/a
    x, y = pos
    x0, y0 = center

    return (
        (x-x0)**2*(c**2 + s**2/r**2) +
        (y-y0)**2*(c**2/r**2 + s**2) +
        2*(x-x0)*(y-y0)*c*s*(1./r**2 - 1.) < (a/3600.)**2
    )


def place_ellipse(a, b, pos, theta, color, label, ls='-', ax=None,
                  labels_cache=set()):
    """
    Draw ellipse on matplotlib axes object
    """
    ax = plt.gca() if ax is None else ax
    if label in labels_cache:
        label = '__nolegend__'
    else:
        labels_cache |= {label}

    el = Ellipse(
        (0, 0),
        width=2*a, height=2*b,
        facecolor='none', edgecolor=color, label=label, ls=ls
    )
    transf = transforms.Affine2D()\
        .rotate_deg(theta) \
        .translate(*pos)

    el.set_transform(transf + ax.transData)
    ax.add_patch(el)


@memory.cache(ignore=['pos'])
def getfits(pos, size, name, filt):
    fitsurl = geturl(*pos, size=size, filters=filt, format="fits")
    print(fitsurl)
    hdu = fits.open(fitsurl[0])[0]
    return hdu.copy()


def plot_panstarrs(center, size, name, filt, df, ref=None,
                   image=True,
                   median=True, kron=True, petrosian=True,
                   exp = False, sersic=True, voculer=False,
                   transform=LinearStretch() + PercentileInterval(99.5),
                   subplotindex=111, fig=None):
    """
    A complex routine to plot get fits image at {center} with {size}
    and draw various fitted ellipses on it

    Parameters
    ----------
    center: `tuple`
        (ra, dec) in degrees
    size: `float`
        size of image in pixels

    Returns
    -------
    ax: `matplotlib.pyplot.Axes`
    """

    if fig is None:
        fig = plt.figure(figsize=(7, 7))

    if image:
        hdu = getfits(center, size, name, filt)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wcs = WCS(hdu.header, fix=False)

        # print(wcs)
        im = hdu.data
        # nan to zeros
        im[np.isnan(im)] = 0.0

        # set contrast to something reasonable
        im_sane = transform(im)

    else:
        # generate a fake WCS
        wcs = WCS(naxis=2)
        cdelt = 1/3600/PANPARS.scale
        wcs.wcs.crpix = [int(size/2), int(size/2)]
        wcs.wcs.cdelt = np.array([cdelt, cdelt])
        wcs.wcs.crval = center
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.pc = [[-1, 0], [0, 1]]
        wcs.wcs.latpole = -30

    ax = plt.subplot(subplotindex, projection=wcs)
    ax.set_title(name)

    if image:
        ax.imshow(im_sane, cmap="bone", origin="lower")
    else:
        ax.imshow(np.zeros((size, size)))

    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())

    labels_cache = set()
    if ref is not None:
        ref_pixcoords = wcs.all_world2pix([(ref['ra'], ref['dec'])], 0)
        ax.scatter(*ref_pixcoords.T, color='gray', marker='x')
        place_ellipse(ref['a']*PANPARS.scale, ref['b']*PANPARS.scale, ref_pixcoords[0], ref['PA'],
                      'gray', 'Reference', ax=ax, ls='--', labels_cache=labels_cache)

    if len(df) > 0:
        panstarrs_src = wcs.all_world2pix(df[['raMean', 'decMean']].values, 0)
        ax.scatter(*panstarrs_src.T, color='yellow', marker='*')

    for i in range(len(df)):
        row = df.iloc[i:i+1]

        number = row.index[0]
        row = row.iloc[0]
        x, y = panstarrs_src[i]
        name = row['objName']

        ax.annotate(number, xy=(x, y),
                    xytext=(10, 0), textcoords="offset points",
                    va="center", ha="left",
                    bbox=dict(boxstyle="round", fc="w", alpha=0.8))

        if median:
            a, b = row[[f"{filt}GalMajor", f"{filt}GalMinor"]]*PANPARS.scale/2
            PA = row[f"{filt}GalPhi"]
            # n = row[f"{filt}GalIndex"]
            place_ellipse(a, b, (x, y), PA, 'red', 'sectormedian', ax=ax,
                          labels_cache=labels_cache)
        if kron:
            kronrad = row[f"{filt}KronRad"]*PANPARS.scale
            place_ellipse(kronrad, kronrad, (x, y), PA, 'yellow', 'Kron', ax=ax,
                          ls='--',
                          labels_cache=labels_cache)
        if sersic:
            sera = row[f"{filt}SerRadius"]*PANPARS.scale
            serab = row[f"{filt}SerAb"]
            serb = sera * serab
            serPA = row[f"{filt}SerPhi"]
            place_ellipse(sera, serb, (x, y), serPA, 'orange', 'Sersic', ax=ax,
                          labels_cache=labels_cache)
        if exp:
            expa = row[f"{filt}ExpRadius"]*PANPARS.scale
            expab = row[f"{filt}ExpAb"]
            expb = expa * expab
            expPA = row[f"{filt}ExpPhi"]
            place_ellipse(expa, expb, (x, y), expPA, 'green', 'Exp', ax=ax,
                          labels_cache=labels_cache)
        if voculer:
            deva = row[f"{filt}ExpRadius"]*PANPARS.scale
            devab = row[f"{filt}ExpAb"]
            devb = deva * devab
            devPA = row[f"{filt}ExpPhi"]
            place_ellipse(deva, devb, (x, y), devPA, 'blue', 'Voculer', ax=ax,
                          labels_cache=labels_cache)
        if petrosian:
            petrorad = row[f"{filt}petRadius"]
            place_ellipse(petrorad, petrorad, (x, y), 0, 'yellowgreen', 'Petro', ax=ax,
                          ls='--', labels_cache=labels_cache)
        # if eff: place_ellipse(reff, reff, (x,y), PA, 'green', '$R_e$ (L/2)', ax=ax, ls='--') 
    plt.legend()
    return ax

    
def ref_from_rfgc(sample):
    """
    rename columns from RFGC catalog
    """
    ref = dict(
        ra = sample['RAJ2000'],
        dec = sample['DEJ2000'],
        a = sample['aO'],
        b = sample['bO'],
        PA = sample['PA']
    )
    return ref


def show_galaxy_rfgc(
    sample, filt, df=None, sortby=None, zoom=1,
    jobs=None, template=None, **kwargs
):
    ref = ref_from_rfgc(sample)
    size = 2 * int(ref["a"] * PANPARS.scale)
    size_arcmin = size / PANPARS.scale / 60 / 2  # radius

    name = "RFGC " + str(int(sample["RFGC"])) + " in " + filt

    if df is None and jobs is not None and template is not None:
        df = cone_galaxy_search(
            jobs, template, (ref["ra"], ref["dec"]), size_arcmin, filt
        )
        print(df.T)
    if sortby is not None:
        df = df.sort_values(by=sortby, ascending=False)
    df.index = range(1, len(df) + 1)

    plot_panstarrs(
        (ref["ra"], ref["dec"]),
        int(size / zoom),
        name,
        filt,
        df=df,
        ref=ref,
        **kwargs
    )
