
import numpy as np
import numpy.ma as ma
import scikits.hydroclimpy as hydro

import scikits.hydroclimpy.io.coaps as coaps
import scikits.hydroclimpy.io.usgs as usgs

import matplotlib.pyplot as pyplot
import scikits.hydroclimpy.plotlib as cpl

weatherdata = coaps.load_coaps_stationdata(90435).fill_missing_dates()
rainfall = weatherdata['rain']
flowdata = usgs.load_usgs_flows('02217770').fill_missing_dates()

rainfall = hydro.adjust_endpoints(rainfall, start_date=flowdata.dates[0])
flowdata = hydro.adjust_endpoints(flowdata, end_date=rainfall.dates[-1])

fig = cpl.hydrograph(rainfall, flowdata, figsize=(12, 6))

fig.hyeto.set_ylabel("Rainfall (mm)", fontweight='bold')
fig.hydro.set_ylabel("Flows (cfs)", fontweight='bold')
fig.suptitle("Hydrograph for the North Oconee River at Athens, GA",
             fontweight="bold", fontsize=12)

fig.hyeto.set_datelimits('2005-01-01', '2005-12-31')

