
import numpy as np
import numpy.ma as ma

import scikits.hydroclimpy as hydro
import scikits.hydroclimpy.io.coaps as coaps
import scikits.hydroclimpy.enso as enso


stationdict = coaps.ids_bystate('GA')
stationid = [v for (k, v) in stationdict.items() if 'Athens' in k.capitalize()]
print stationid

data = coaps.load_coaps_data(90435)
rainfall = data['rain']
print rainfall.freqstr
print rainfall.dates[[0,-1]]

mrainfall = rainfall.convert('M',func=ma.sum)
arainfall = mrainfall.convert('A')
assert(arainfall.shape == (64, 12))

monthly_means = arainfall.mean(axis=0).round(1)
print monthly_means
"""   [ 115.9  112.4  130.8   94.1   99.4  103.5  122.3   90.2   94.5   77.8
         92.7   98.6]"""


ONI = enso.load_oni()
print ONI.dates[[0,-1]]
"""   [Jan-1950 Dec-2008]"""

mrainfall = enso.set_ensoindicator(mrainfall, ONI)
print type(mrainfall)
assert(isinstance(mrainfall, enso.ClimateSeries))
print mrainfall.ensoindicator.dates[[0,-1]]
"""   [Jan-1944 Dec-2007]"""

mrainfall_2K = mrainfall[(mrainfall.year == 2000)]
print mrainfall_2K.cold
"""[113.792 50.546 86.614 43.18 54.864 50.292 -- -- -- 5.842 106.68 87.884]"""
print mrainfall_2K.neutral
"""[-- -- -- -- -- -- 85.344 93.218 122.174 -- -- --]"""
print mrainfall_2K.warm
"""[-- -- -- -- -- -- -- -- -- -- -- --]"""


mdtype = [('cold', float), ('neutral', float), ('warm', float),
          ('global', float)]
monthly_means_enso = ma.empty(12, dtype=mdtype)
monthly_means_enso['global'] = monthly_means

for phase in ('cold', 'neutral', 'warm'):
    mcurrent = getattr(mrainfall, phase)
    acurrent = mcurrent.convert('A')
    monthly_means_enso[phase] = acurrent.mean(axis=0).round(1)


