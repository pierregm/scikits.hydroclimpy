import warnings

# Try importing mpl_addons : that should fail if matplotlib is not installed
import mpl_addons
from mpl_addons import *
import flowplots
from flowplots import *
mpl_activated = True

# Try to import maptools....
try:
    import maptools
    from maptools import *
except ImportError:
    msg = "The basemap module is not installed: mapping functions deactivated"
    warnings.warn(msg)

from matplotlib.pyplot import acorr, annotate, arrow, autumn, axes, axhline, \
                              axhspan, axis, axvline, axvspan, bar, barbs, \
                              barh, bone, box, boxplot, broken_barh, cla, \
                              clabel, clf, clim, close, cm, cohere, colorbar, \
                              colormaps, colors, connect, contour, contourf, \
                              cool, copper, csd, dedent, delaxes, disconnect, \
                              draw, draw_if_interactive, errorbar, figaspect, \
                              figimage, figlegend, figtext, figure, fill, \
                              findobj, flag, gca, gcf, gci, get, get_backend, \
                              get_cmap, get_current_fig_manager, \
                              get_plot_commands, get_scale_docs, \
                              get_scale_names, getp, ginput, gray, grid, \
                              hexbin, hist, hlines, hold, hot, hsv, imread, \
                              imshow, interactive, ioff, ion, is_numlike, \
                              is_string_like, ishold, isinteractive, jet, \
                              legend, loglog, matplotlib, matshow, mlab, \
                              new_figure_manager, normalize, over, pcolor, \
                              pcolormesh, pie, pink, plot, plot_date, \
                              plotfile, plotting, polar, prism, psd, \
                              pylab_setup, quiver, quiverkey, rc, \
                              rcParams, rcParamsDefault, rcdefaults, rgrids, \
                              savefig, scatter, sci, semilogx, semilogy, \
                              setp, show, silent_list, specgram, spectral, \
                              spring, spy, stem, step, subplot, subplot_tool, \
                              subplots_adjust, summer, suptitle, \
                              switch_backend, table, text, thetagrids, \
                              title, twinx, twiny, vlines, \
                              waitforbuttonpress, winter, xcorr, xlabel, \
                              xlim, xscale, xticks, ylabel, ylim, yscale, \
                              yticks
import scikits.timeseries.lib.plotlib
from scikits.timeseries.lib.plotlib import *
