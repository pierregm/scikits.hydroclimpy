
import numpy as np
import scikits.hydroclimpy as hydro
import scikits.hydroclimpy.enso as enso
import scikits.hydroclimpy.io.usgs as usgs
import scikits.hydroclimpy.plotlib as cpl


series = usgs.load_usgs_flows('02217770')

fig = cpl.figure()
fsp = fig.add_subplot(111)

cpl.plot_fdc(series, ax=fsp, lw=2, ls='-', c='k', zorder=10)
cpl.plot_fdc(series, ax=fsp, starting_month=4)
fsp.set_ylabel("Flows (cfs)", fontweight='bold')
fig.suptitle("Flow duration curve for the North Oconee River at Athens, GA",
             fontweight="bold", fontsize=12)
