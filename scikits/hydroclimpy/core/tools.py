"""
.. currentmodule:: scikits.hydroclimpy.core.tools

The :mod:`scikits.hydroclimpy.core.tools` module defines various classes and
function that can be used in a wide range of situations:



.. autosummary::
   :nosignatures:

   Bag
   Cluster

   flatten_sequence
   flatargs
   timed_question



:class:`Bag` object
===================

.. autoclass:: Bag



:class:`Cluster` object
=======================

.. autoclass:: Cluster



Convenience Functions
=====================

.. autofunction:: flatten_sequence
.. autofunction:: flatargs
.. autofunction:: timed_question


"""
import sys

import numpy as np
import numpy.ma as ma

#from scikits.timeseries._tools import deprecated, deprecated_for, docwrapper

__all__ = [#'deprecated', 'deprecated_for',
           'flatten_sequence',
           'Bag',
           'Cluster']


def flatten_sequence(iterable):
    """
    Flattens a compound of nested iterables, and returns a generator object
    of the flattened sequence

    Examples
    --------
    >>> x = flatten_sequence([0,[1,2,[3,4]],[5,6]])
    >>> [_ for _ in x]
    [0,1,2,3,4,5,6,]

    """
    itm = iter(iterable)
    for elm in itm:
        if hasattr(elm, '__iter__') and not isinstance(elm, basestring):
            for f in flatten_sequence(elm):
                yield f
        else:
            yield elm

def flatargs(*args):
    """
    Flattens the arguments.

    """
    if not hasattr(args, '__iter__'):
        return args
    else:
        return flatten_sequence(args)



class Bag(dict):
    """
    Generic extension to dictionaries, where elements can be accessed both by
    keys and attributes.

    The parameters required to create a :class:`Bag` can be:

    * A regular dictionary or another :class:`Bag`.
    * A list of tuples (key, value).
    * A list of optional arguments.

    As a subclass of :class:`dict`, a :class:`Bag` object inherits
    its attributes and methods.

    Examples
    --------
    >>> bag = Bag(a=1, b=2, c=3)
    {'a': 1, 'c': 3, 'b': 2}
    >>> bag = Bag(bag)
    {'a': 1, 'c': 3, 'b': 2}
    >>> bag = Bag([('a',1), ('b',2), ('c',3)])
    {'a': 1, 'c': 3, 'b': 2}

    See Also
    --------
    dict
        ... for the description of inherited methods.
    """
    #
    def __init__(self, obj=None, **kw):
        obj = dict(obj or {})
        obj.update(kw)
        dict.__init__(self, **obj)
        self.__dict__['_dict'] = obj
    #
    def __getattr__(self, item):
        _dict = self.__dict__.get('_dict', {})
        if item in _dict:
            return _dict[item]
        return dict.__getattr__(self, item)
    #
    def __setattr__(self, item, value):
        _dict = self.__dict__['_dict']
        # Don't overwrite the attributes of the dictionary
        if getattr(_dict, item, getattr(self, item, None)) is None:
            self[item] = _dict[item] = value
    #
    def set_defaults(self, defaults=None, **kwargs):
        """
    Sets the defaults of a :class:`Bag` instance.

    Parameters
    ----------
    default : dictionary, optional
        A dictionary of defaults
        """
        defaults = dict(defaults or {})
        defaults.update(kwargs)
        missing = dict([(k, v) for (k, v) in defaults.items() if k not in self])
        self.update(missing)



class Cluster(object):
    """
    Groups consecutive data from an array according to a clustering condition.

    A cluster is defined as a group of consecutive values differing by at most the
    increment value.

    Missing values are **not** handled: the input sequence must therefore be free 
    of missing values.

    Parameters
    ----------
    darray : ndarray
        Input data array to clusterize.
    increment : {float}, optional
        Increment between two consecutive values to group.
        By default, use a value of 1.
    operator : {function}, optional
        Comparison operator for the definition of clusters.
        By default, use :func:`numpy.less_equal`.


    Attributes
    ----------
    :attr:`inishape` : tuple
        Shape of the argument array (stored for resizing).
    :attr:`inisize` : integer
        Size of the argument array.
    :attr:`uniques` : sequence
        List of unique cluster values, as they appear in chronological order.
    :attr:`slices` : sequence
        List of the slices corresponding to each cluster of data.
    :attr:`starts` : ndarray
        Array of the indices at which the clusters start.
    :attr:`clustered` : list
        List of clustered data.

    Methods
    -------
    .. automethod:: markonsize
    .. automethod:: mark_greaterthan
    .. automethod:: grouped_slices
    .. automethod:: grouped_limits


    Examples
    --------
    >>> A = [0, 0, 1, 2, 2, 2, 3, 4, 3, 4, 4, 4]
    >>> klust = cluster(A,0)
    >>> [list(_) for _ in klust.clustered]
    [[0, 0], [1], [2, 2, 2], [3], [4], [3], [4, 4, 4]]
    >>> klust.uniques
    array([0, 1, 2, 3, 4, 3, 4])
    
    >>> x = [ 1.8, 1.3, 2.4, 1.2, 2.5, 3.9, 1. , 3.8, 4.2, 3.3, 
    ...       1.2, 0.2, 0.9, 2.7, 2.4, 2.8, 2.7, 4.7, 4.2, 0.4]
    >>> Cluster(x,1).starts
    array([ 0,  2,  3,  4,  5,  6,  7, 10, 11, 13, 17, 19])
    >>> Cluster(x,1.5).starts
    array([ 0,  6,  7, 10, 13, 17, 19])
    >>> Cluster(x,2.5).starts
    array([ 0,  6,  7, 19])
    >>> Cluster(x,2.5,greater).starts
    array([ 0,  1,  2,  3,  4,  5,  8,  9, 10, 
    ...    11, 12, 13, 14, 15, 16, 17, 18])
    >>> y = [ 0, -1, 0, 0, 0, 1, 1, -1, -1, -1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]
    >>> Cluster(y,1).starts
    array([ 0,  1,  2,  5,  7, 10, 12, 16, 18])

    """
    def __init__(self, darray, increment=1, operator=np.less_equal):
        """
        """
        if hasattr(darray, 'mask') and darray.mask.any():
            raise ma.MAError("Masked arrays should be filled prior clustering.")
        else:
            darray = np.asanyarray(darray)
        n = darray.size
        self.inishape = darray.shape
        self.inisize = darray.size
        darray = darray.ravel()
        clustercond = 1 - operator(np.absolute(np.diff(darray)), increment)
        sid = np.concatenate(([0, ],
                              np.arange(1, n).compress(clustercond),
                              [n, ]))
        slobj = np.asarray([slice(i, d)
                            for (i, d) in np.broadcast(sid[:-1], sid[1:])])
        #
        self.uniques = darray[sid[:-1]]
        self.clustered = [darray[k] for k in slobj]
        self.sizes = np.asarray(np.diff(sid))
        self.slices = slobj
        self.starts = sid[:-1]

    def markonsize(self, operator, sizethresh):
        """
    Creates a **mask** for the clusters that do not meet a size requirement.
    Thus, outputs ``False`` if the size requirement is met, ``True`` otherwise.
    
    Parameters
    ----------
    operator : function
        Comparison operator
    sizethresh : float
        Requirement for the sizes of the clusters

        """
        resmask = np.empty(self.inisize, dtype=bool)
        resmask[:] = True
#        for k in self.slices.compress(operator(self.sizes,sizethresh)):
        for k in self.slices[operator(self.sizes, sizethresh)]:
            resmask[k] = False
        return resmask.reshape(self.inishape)

    def mark_greaterthan(self, sizemin):
        """
    Shortcut for :meth:`markonsize(greater_equal,sizemin)`.
    Thus, the command outputs ``False`` for clusters larger than ``sizemin``, and
    ``True`` for clusters smaller than ``sizemin``.
    
    Parameters
    ----------
    sizemin : int
        Minimum size of the clusters.

    See Also
    --------
    :meth:`markonsize`
        Creates a **mask** for the clusters that do not meet a size requirement.
    """
        return self.markonsize(np.greater_equal, sizemin)

    def grouped_slices(self):
        """
    Returns a dictionary with the unique values of ``self`` as keys, and a list
    of slices for the corresponding values.

    See Also
    --------
    Cluster.grouped_limits
        that does the same thing
        """
        #
        uniques = self.uniques.view(np.ndarray)
        output = dict([(k, []) for k in np.unique1d(uniques)])
        for (k, v) in zip(self.uniques, self.slices):
            output[k].append(v)
        return output

    def grouped_limits(self):
        """
    Returns a dictionary with the unique values of ``self`` as keys, and a list
    of tuples (starting index, ending index) for the corresponding values.

    See Also
    --------
    Cluster.grouped_slices
        """
        output = dict([(k, []) for k in np.unique1d(self.uniques)])
        for (k, v) in zip(self.uniques, self.slices):
            output[k].append((v.start, v.stop))
        for k in output:
            output[k] = np.array(output[k])
        return output



###--- Timers ---

def alarm_handler(*args):
    """Raise a `TimeOut` exception."""
    raise Exception("TimeOut")

def timed_question(prompt, timeout=5, default_answer=''):
    """
    Waits a given delay for an answer to the prompted question.
    If no answer is given before the delay is over, outputs a predefined answer.

    Parameters
    ----------
    prompt : string
        Question asked at the prompt.
    timeout : {int}, optional
        Waiting time before the default answer is output(in s).
    default_answer : {string}, optional
        Default answer.

    """
    import signal
    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(timeout)
    sys.stdout.write(prompt)
    try:
        text = sys.stdin.readline()
    except KeyboardInterrupt:
        print "EXIT!"
        sys.exit()
    except:
        text = default_answer
    signal.alarm(0)
    return text



