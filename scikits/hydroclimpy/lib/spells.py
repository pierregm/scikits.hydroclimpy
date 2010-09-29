"""
.. .. module:: scikits.hydroclimpy.lib.spells
.. currentmodule:: scikits.hydroclimpy.lib.spells

=======================
Wet/dry spells analysis
=======================


This module introduces utilities to analyze dry and wet spells in a
precipitation time series.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   clump_masked
   clump_unmasked


.. autosummary::
   :nosignatures:
   :toctree: generated/

   dry_spells
   wet_spells

"""

import itertools

import numpy as np
import numpy.ma as ma
from numpy.ma import masked, nomask

import scikits.hydroclimpy.enso as enso


def _ezclump(mask):
    """
    Finds the clumps (groups of data with the same values) for a 1D bool array.

    Returns a series of slices.
    """
    #def clump_masked(a):
    if mask.ndim > 1:
        mask = mask.ravel()
    idx = (mask[1:] - mask[:-1]).nonzero()
    idx = idx[0] + 1
    slices = [slice(left, right)
              for (left, right) in zip(itertools.chain([0], idx),
                                       itertools.chain(idx, [len(mask)]),)]
    return slices


def clump_unmasked(a):
    """
    Returns a list of slices corresponding to the unmasked clumps of a 1D array.

    Examples
    --------
    >>> a = ma.masked_array(np.arange(10))
    >>> a[[0, 1, 2, 6, 8, 9]] = ma.masked
    >>> clump_unmasked(a)
    [slice(3, 6, None), slice(7, 8, None)]

    .. versionadded:: 1.4.0
    """
    mask = getattr(a, '_mask', nomask)
    if mask is nomask:
        return [slice(0, a.size)]
    slices = _ezclump(mask)
    if a[0] is masked:
        result = slices[1::2]
    else:
        result = slices[::2]
    return result
np.ma.clump_unmasked = clump_unmasked

def clump_masked(a):
    """
    Returns a list of slices corresponding to the masked clumps of a 1D array.

    Examples
    --------
    >>> a = ma.masked_array(np.arange(10))
    >>> a[[0, 1, 2, 6, 8, 9]] = ma.masked
    >>> clump_masked(a)
    [slice(0, 3, None), slice(6, 7, None), slice(8, None, None)]

    .. versionadded:: 1.4.0
    """
    mask = ma.getmask(a)
    if mask is nomask:
        return []
    slices = _ezclump(mask)
    if len(slices):
        if a[0] is masked:
            slices = slices[::2]
        else:
            slices = slices[1::2]
    return slices
ma.clump_masked = clump_masked



def dry_spells(rainfall, threshold=0.254, start=True):
    """
    Compute the dry spells for a series of precipitations

    Parameters
    ----------
    rainfall : TimeSeries
        TimeSeries of precipitations.
    threshold : float, optional
        Minimum amount of precipitation defining a wet day.
    start : boolean, optional
        Whether the spells are associated with the first or last day of the spell.
    

    Returns
    -------
    dry_spells : TimeSeries
        A :class:`TimeSeries` giving the duration of the spell at either the first 
        or last date.
    """
    rdates = getattr(rainfall, 'dates', None)
    rdata = ma.masked_array(rainfall, subok=False)
    condition = ma.masked_greater(rdata, threshold)
    slices = ma.clump_unmasked(condition)
    # Get the durations and starting dates of each spell
    durations = [s.stop - s.start for s in slices]
    # Give up if we have no dates
    if rdates is None:
        return np.array(durations)
    # Get the dates, then...
    if start:
        dates = rdates[[s.start for s in slices]]
    else:
        dates = rdates[[s.stop - 1 for s in slices]]
    ensoi = getattr(rainfall, 'ensoindicator', None)
    spells = enso.climate_series(durations, dates=dates, ensoindicator=ensoi)
    return spells



def wet_spells(rainfall, threshold=0, start=True, mode='both'):
    """
    Compute the wet spells for a series of precipitations

    Parameters
    ----------
    rainfall : TimeSeries
        TimeSeries of precipitations.
    threshold : float, optional
        Minimum amount of precipitation defining a wet day.
    start : boolean, optional
        Whether the spells are associated with the first or last day of the spell.
    mode : {'durations', 'intensities', 'all'}
        Whether to return only durations, intensities or both.
    

    Returns
    -------
    wet_spells : TimeSeries
        A :class:`TimeSeries` giving the duration and/or intensity of a spell
        at either the first or last date.
    """
    rdates = getattr(rainfall, 'dates', None)
    rdata = ma.masked_array(rainfall, subok=False)
    condition = ma.masked_less(rdata, threshold)
    slices = ma.clump_unmasked(condition)
    # Get the durations and starting dates of each spell
    mode = (mode or 'both').lower()[0]
    if mode == 'd':
        # Durations
        result = np.array([s.stop - s.start for s in slices], dtype=int)
    elif mode == 'i':
        # Intensities
        result = np.array([rdata[s].sum() for s in slices])
    else:
        durations = [s.stop - s.start for s in slices]
        intensities = [rdata[s].sum() for s in slices]
        result = np.array(zip(durations, intensities),
                          dtype=[('durations', int), ('intensities', float)],)
    if rdates is None:
        return result
    if start:
        dates = rdates[[s.start for s in slices]]
    else:
        dates = rdates[[s.stop - 1 for s in slices]]
    ensoi = getattr(rainfall, 'ensoindicator', None)
    if mode in 'id':
        spells = enso.climate_series(result, dates=dates, ensoindicator=ensoi)
    else:
        spells = enso.climate_records(result, dates=dates, ensoindicator=ensoi)
    return spells



