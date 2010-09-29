
import numpy as np
import scikits.hydroclimpy as hydro
import scikits.hydroclimpy.enso as enso
import scikits.hydroclimpy.io.usgs as usgs
import scikits.hydroclimpy.plotlib as cpl

ONI = enso.load_oni()
series = usgs.load_usgs_flows('02217770').fill_missing_dates()
series = enso.set_ensoindicator(series, ONI)

fig = cpl.figure()
fsp = fig.add_subplot(111)
cpl.plot_fdc(series, plot_enso=True, marker='o', markersize=6, ax=fsp)
cpl.plot_fdc(series, plot_enso=True, starting_month=4, ax=fsp)
fsp.legend()
fsp.set_ylabel("Flows (cfs)", fontweight='bold')
fig.suptitle("Flow duration curve for the North Oconee River at Athens, GA",
             fontweight="bold", fontsize=12)

