"""
Computing potential evapotranspiration.

Please consult the documentation for a more exhaustive presentation.
"""

__all__ = ['SolarInformation', 'PotentialEvapoTranspiration', 'PETModel',
           'aerodynamic_resistance']

import datetime
import warnings

import numpy as np
import numpy.ma as ma

import scikits.hydroclimpy as hydro
from scikits.hydroclimpy import force_reference, \
                                Date, DateArray, TimeSeries, const as _c, \
                                date_array, time_series

import scikits.timeseries.const as _c

from scikits.timeseries.lib.interpolate import backward_fill

_doc_parameters = dict(
z="""z : float, optional
        Local elevation above the sea level [m].
""",
albedo="""albedo : float, optional
        Albedo.
        By default, a value of ``0.23`` is assumed for grass/alfalfa references.
""",
a_s="""a_s : float, optional
        Regression constant expressing the fraction of extraterrestrial solar
        radiation reaching the ground on overcast days.
        By default, ``a_s = 0.25``.
""",
b_s="""b_s : float, optional
       Regression constant related to the amount of extraterrestrial solar
       radiation reaching the ground on clear day.
       By default, ``b_s = 0.5``.
""",
RHmin="""RHmin : TimeSeries, optional
        Minimum relative humidity ``[0-1]``.
""",
RHmax="""RHmax : TimeSeries, optional
        Maximum relative humidity ``[0-1]``.
""",
RHmean="""RHmean : TimeSeries, optional
        Mean relative humidity ``[0-1]``.
""",
crop_ref="""crop_ref : {'short','tall'}, optional
        String giving the kind of reference crop. 
        Use ``'short'`` for clipped grass (default), ``'tall'`` for alfalfa.
""",
u2="""u2 : float, optional
        Mean wind speed at an elevation of 2m [m.s\ :sup:`-1`]
""",
Kt="""Kt : {None, float}, optional
        Regression constant used to approximate solar radiations from temperatures. 
        The default (:math:`K_t=0.171` corresponds to the standard value suggested
        by Hargreaves and Samani [Hargreaves_Samani_1985]_.
        If None, :math:`K_t` is evaluated from the temperatures range :math:`\Delta T`
        as :math:`K_t = 0.00185 {\Delta T}^2 - 0.0433 {\Delta T} + 0.4023`.
        [REF_TO_FIND_19xx]_
""",
Rs="""Rs : TimeSeries, optional
        Incoming solar radiations [MJ.m\ :sup:`-2`.d\ :sup:`-1`].
        If None, ``Rs`` is estimated from extraterrestrial solar radiations 
        (``Ra``) and minimum and maximum temperatures as 
        :math:`K_{t} \sqrt{\Delta\;T} Ra`.
""",
result="""PET : TimeSeries
        Potential Evapo-Transpiration [mm.d\ :sup:`-1`].
""",
reference="""reference : string, optional
        Name of the PET model used as reference.
        If ``None``, use the model defined in :attr:`reference`.
""",
models="""models : {string, sequence}, optional
        A list of strings corresponding to the models to investigate.
        Each entry should be part of the :attr:`models` list.
        If one or several antries do not satisfy this requirement, they are 
        simply discarded.
        If ``None``, the list defaults to the :attr:`models` attribute.
""",
freq="""freq : {string, int}, optional
        A valid frequency specifier.
        If ``None``, the parameter is set to the :attr:`freq` attribute of
        the PET estimates.
""",
)



def _checktimeseries(*series):
    """
    Checks the compatibility of two series.

    """
    for _ in series:
        if not isinstance(_, TimeSeries):
            raise TypeError(
                  "Input temperatures should be valid TimeSeries objects!")
    series = hydro.align_series(*series)
    freq = series[0].freq
    if freq < _c.FR_MTH:
        raise ValueError("The input frequency should be at most monthly!")
    elif freq >= _c.FR_HR:
        raise ValueError("The input frequency should be at least daily!")
    if len(series) == 1:
        series = series[0]
    return series



def frac_day_of_year(date, offset=False):
    """
    Returns the fractional day of year, where 0 or 1./365 represents Jan., 1st
    and 1. represents Dec., 31st.

    Parameters
    ----------
    date : var
        A :class:`~scikits.timeseries.Date`, :class:`~scikits.timeseries.DateArray` or
        :class:`datetime.datetime` object.
    offset : {False, True}, optional
        Whether the first day of year corresponds to 1./365 (False) or 0 (True)

    """
    if isinstance(date, (Date, datetime.datetime)):
        shp = (1,)
    elif isinstance(date, DateArray):
        shp = len(date)
    else:
        raise TypeError("Unrecognizable input type: got '%', expected a "\
                        "Date, DateArray or datetime.datetime object" % \
                        type(date))
    years = date.year
    yrdays = 365 * np.ones(shp, dtype=float)
    yrdays[(years % 4 == 0) | (years % 400 == 0)] = 366
    freq = getattr(date, 'freq', _c.FR_DAY)
    if freq >= _c.FR_WK:
        doy = date.day_of_year
    elif freq <= _c.FR_QTR:
        raise ValueError("The frequency should be at most monthly!")
    else:
        doy = np.array(date.month / 12. * 365 - 15, dtype=float)
    if offset:
        doy -= 1
    return doy / yrdays




class SolarInformation(object):
    """
    Defines some astronomical information for a series of dates at a given latitude.

    The series of dates must be a valid :class:`~scikits.timeseries.DateArray`
    object.
    If its frequency is weekly or less (monthly, quarterly), the series is 
    internally converted to a daily frequency, but its original frequency is
    saved in the :attr:`freq` attribute. 
    A :exc:`NotImplementedError` is raised if the initial frequency is hourly
    or higher.

    All of the attributes are stored internally with a daily frequency.
    Public versions are available: in that case, the attribute is reconverted
    to the initial frequency :attr:`freq`.

    Parameters
    ----------
    dates : DateArray
        Series of dates on which to compute the astronomical information.
    latitude : float
        Latitude of the site (in degrees).


    """
    #
    def __init__(self, dates, latitude):
        if not isinstance(dates, DateArray):
            raise TypeError("The input should be a valid DateArray object "\
                            "(got '%s' instead)" % type(dates))
        self.freq = dates.freq
        if self.freq < _c.FR_DAY:
            self._dates = date_array(start_date=dates[0].asfreq('D', 'S'),
                                     end_date=dates[-1].asfreq('D', 'E'),
                                     freq='D')
            self._convertfunc = ma.mean
        elif self.freq > _c.FR_DAY:
            err_msg = "Solar information is not available for series with a "\
                      "frequency higher than daily."
            raise NotImplementedError(err_msg)
        else:
            self._dates = dates
            self._convertfunc = None
        self._latitude = latitude * np.pi / 180.
        #
        self._etrad = None
        self._csrad = None
        # Earth-Sun distance
        esa = 2 * np.pi * frac_day_of_year(self._dates, offset=False)
        self._earth_sun_distance = 1. + 0.033 * np.cos(esa)
        # Solar declination
        fyr = 2 * np.pi * frac_day_of_year(self._dates, offset=True)
        solar_declination = (0.006918 - \
                             0.399912 * np.cos(fyr) + 0.070257 * np.sin(fyr) - \
                             0.006758 * np.cos(2 * fyr) + 0.000907 * np.sin(2 * fyr) - \
                             0.002697 * np.cos(3 * fyr) + 0.001480 * np.sin(3 * fyr))
        self._solar_declination = solar_declination
        # Solar angle
        cos_ws = -np.tan(solar_declination) * np.tan(self._latitude)
        self._sunset_angle = np.arccos(np.clip(cos_ws, -1., 1.))
        # Daylight hours
        self._daylighthours = self._sunset_angle * 24 / np.pi
        # Extraterrestrial solar radiation [MJ/m2/d]
        solconst = 4.921 * 24
        sunangle = self._sunset_angle
        latitude = self._latitude
        etsolrad = (sunangle * np.sin(latitude) * np.sin(solar_declination) + \
                    np.cos(latitude) * np.cos(solar_declination) * np.sin(sunangle))
        self._etrad = (1. / np.pi * solconst * self._earth_sun_distance) * etsolrad
        self._csrad = 0.75 * self._etrad


    def _get_latitude(self):
        """
    Returns the latitude of the site (in degrees).
    The latitude is stored in radians internally.
        """
        return self._latitude * 180. / np.pi
    latitude = property(fget=_get_latitude, doc=_get_latitude.__doc__)


    def _get_solar_declination(self):
        """
    Returns the solar declination :math:`\delta` as [Spencer_1971]:
    
    .. math::

       \\frac{\pi}{180^\circ} \delta &= 0.006918 \\\\
                                     &- 0.399912 \cos(\phi)  + 0.070257 \sin(\phi) \\\\
                                     &- 0.006758 \cos(2\phi) + 0.000907 \sin(2\phi) \\\\
                                     &- 0.002697 \cos(3\phi) + 0.00148 \sin(3\phi)

    where :math:`\phi= 2 \pi (n-1)/N` is the fractional year in radians,
    :math:`n` the day of the year (1 for Jan, 1st) and :math:`N` the number of days 
    in the year.
        """
        out = time_series(self._solar_declination, dates=self._dates)
        if self.freq < _c.FR_DAY:
            return out.convert(self.freq, self._convertfunc)
        else:
            return out
    solar_declination = property(fget=_get_solar_declination,
                                 doc=_get_solar_declination.__doc__)


    def _get_earth_sun_distance(self):
        """
    Returns the average inverse relative distance Earth-sun :math:`\delta_r`
    calculated as:
    
    .. math::
       \delta_r = 1 + 0.033 \cos\left(2\pi f\\right)

    where :math:`f` is the fractional day of the year (1./365 for Jan., 1st and
    1. for Dec.,31st).
    
        """
        esd = time_series(self._earth_sun_distance, dates=self._dates)
        if self.freq < _c.FR_DAY:
            return esd.convert(self.freq, self._convertfunc)
        else:
            return esd
    earth_sun_distance = property(fget=_get_earth_sun_distance,
                                  doc=_get_earth_sun_distance.__doc__)


    def _frac_day_of_year(self, offset=False):
        """
    Returns the fractional days of the year, as a ndarray.
    The fractional day of the year is a float between 0 (or 1./365) for Jan., 1st
    and 1. for Dec., 31st.

    Parameters
    ----------
    offset : {False, True}
        Whether the first day of year corresponds to 1./365 (False) or 0 (True)
        """
        return frac_day_of_year(self._dates, offset)
    #
    def frac_day_of_year(self, offset=False):
        """
    Returns the fractional days of the year, as a :class:`~scikits.timeseries.TimeSeries`.
    The fractional day of the year is a float between 0 (or 1./365) for Jan., 1st
    and 1. for Dec., 31st.

    Parameters
    ----------
    offset : {False, True}
        Whether the first day of year corresponds to 1./365 (False) or 0 (True)
        """
        fracday = time_series(self._frac_day_of_year(offset), dates=self._dates)
        if self.freq < _c.FR_DAY:
            return fracday.convert(self.freq, func=self._convertfunc)
        return fracday


    def _get_sunset_angle(self):
        """
    Returns the sunset hour angle, as a :class:`~scikits.timeseries.TimeSeries`.
    The sunset hour angle :math:`\omega_s` is calculated from the solar declination
    :math:`\delta` and the latitude :math:`\Phi` as
    
    .. math::
       \cos(\omega_s) = -\\tan(\Phi) \cdot \\tan(\delta)

        """
        out = time_series(self._sunset_angle, dates=self._dates)
        if self.freq < _c.FR_DAY:
            return out.convert(self.freq, func=self._convertfunc)
        return out
    sunset_angle = property(fget=_get_sunset_angle,
                            doc=_get_sunset_angle.__doc__)


    def _get_daylength(self):
        """
    Returns the number of hours of daylights, as a :class:`~scikits.timeseries.TimeSeries`.
    The day length :math:`D` depends on the sunset hour angle :math:`\omega_s`
    through the relation :math:`D = 24 \cdot \omega_s / \pi`.

        """
        out = time_series(self._daylighthours, dates=self._dates)
        if self.freq < _c.FR_DAY:
            return out.convert(self.freq, func=self._convertfunc)
        return out
    daylength = property(fget=_get_daylength, doc=_get_daylength.__doc__)
    daylighthours = property(fget=_get_daylength,
                             doc="Alias to :attr:`daylength`.")


    def _get_etrad(self):
        """
    The extraterrestrial solar radiation :math:`R_a` is the solar radiation
    received at the top of the atmosphere on an horizontal surface.
    It depends on the solar constant :math:`\sigma`, the solar radiation hitting
    a surface perpendicular to the sun's rays at the top of the atmosphere,
    and the angle between the sun rays and the normal to the atmosphere.
    
    An expression of :math:`R_a`  [MJ.m\ :sup:`-2`.d\ :sup:`-1`] is given by
    Duffie & Beckman (1981) [Duffie_Beckman_1981]_:
    
    .. math::
       R_a = \\frac{\sigma \delta_r}{\pi} \cdot
             \left[\omega_s \sin(\Phi) \sin(\delta_s) + \cos(\Phi) \cos(\delta_s) \sin(\omega_s)\\right]

    where :math:`\omega_s` is the sunset hour angle (:attr:`sunset_angle`),
    :math:`\delta_r` the inverse relative Earth-Sun distance(:attr:`earth_sun_distance`),
    :math:`\delta_s` the solar declination (:attr:`solar_declination`),
    :math:`\Phi` the site latitude (:attr:`latitude`).


        """
        out = time_series(self._etrad, dates=self._dates)
        if self.freq < _c.FR_DAY:
            return out.convert(self.freq, func=self._convertfunc)
        return out
    extraterrestrial_solar_radiations = property(fget=_get_etrad,
                                                 doc=_get_etrad.__doc__)


    def clearsky_solar_radiations(self, z=0):
        """
    Returns the clear sky solar radiation :math:`R_{so}`
    [MJ.m\ :sup:`-2`.d\ :sup:`-1`], as:
    
    .. math::
       R_{so} = (0.75 + 2 10^{-5} z) R_a

    Parameters
    ----------
    z : {float}, optional
        Elevation from sea  level [m]

        """
        if z:
            csrad = self._csrad + 2e-5 * z * self._etrad
        else:
            csrad = self._csrad
        out = time_series(csrad, dates=self._dates)
        if self.freq <= _c.FR_DAY:
            return out.convert(self.freq, func=self._convertfunc)
        return out


    def solar_radiations_from_sunshinehours(self, sunshinehours,
                                            a_s=0.25, b_s=0.50):
        """
    Computes the incoming solar (shortwave) radiations :math:`Rs` from actual
    sunshine hours as :
    
    .. math::
       R_{s} = \left(a_s + b_s \\frac{D_{act}}{D}\\right) R_a

    where :math:`R_a` is the extraterrestrial radiation
    (:attr:`extraterrestrial_solar_radiations`); 
    :math:`D_{act}/D` is the relative sunshine duration as the ratio of the 
    actual duration of sunshine (:math:`D_{act}`) on the maximum possible
    daylight hours (:math:`D`);
    and :math:`a_s` and :math:`b_s` two regression constants.

    In absence of clouds, :math:`R_s = (a_s+b_s) R_a \\approx 0.75 R_a`.

    Parameters
    ----------
    sunshinehours : TimeSeries
        Average actual hours of sunshine [hr].
    %(a_s)s
    %(b_s)s

        """
        if not isinstance(sunshinehours, TimeSeries):
            raise TypeError("The input should be a valid TimeSeries!")
        freq = sunshinehours.freq
        if freq != _c.fr_DAY:
            sunshinehours = backward_fill(sunshinehours.convert('D', 'E'))
        effsunshine = sunshinehours / self._daylighthours
        sol_rad = (a_s + b_s * effsunshine) * self._etrad
        if freq != _c.FR_DAY:
            return sol_rad.convert(freq, func=ma.mean)
        return sol_rad


    def solar_radiations_from_temperatures(self, tmin, tmax, Kt=0.170):
        u"""
    Computes the incoming solar radiations [MJ.m\ :sup:`-2`.d\ :sup:`-1`] from
    temperatures, as:
    
    .. math::
       R_s = K_t \cdot \sqrt{T_{max}-T_{min}} \cdot R_a

    where :math:`K_t` is an empirical parameter.
    By default, :math:`K_t=0.171`. Otherwise, :math:`K_t` is estimated from the
    temperature range :math:`\Delta T=T_{max}-T_{min}` as

    .. math::
       K_t = 0.00185 (\Delta T)^2 - 0.0433 \Delta T + 0.4023

    Parameters
    ----------
    tmin, tmax : TimeSeries
        Minimum and maximum temperatures [\u00B0C].
    %(Kt)s

        """
        etrad = self.extraterrestrial_solar_radiations
        (tmin, tmax) = hydro.align_with(etrad, *_checktimeseries(tmin, tmax))
        trange = (tmax - tmin)
        if Kt is None:
            Kt = 0.00185 * trange ** 2 - 0.0433 * trange + 0.4023
        if trange.ndim > 1:
            etrad = etrad[:, None]._series
        sol_rad = Kt * ma.sqrt(trange) * etrad
        return sol_rad
    solar_radiations_from_temperatures.__doc__ = \
    (solar_radiations_from_temperatures.__doc__ or '') % _doc_parameters


def aerodynamic_resistance(uz, zw=2, zh=2, h=0.12):
    """
    Computes the aerodynamic resistance for neutral stability conditions.

    Parameters
    ----------
    uz : float
        Windspeed at height zw [m.s^{-1}]
    zw : {2., float}, optional
        Height of wind measurement [m]. A default of 2m is assumed.
    zh : {2., float} optional
        Height of humidity and air temperature measurements [m]
    h : {0.12, float} optional
        Mean height of vegetation [m]
    """
    # d : zero plane displacement height [m]
    d = 2. / 3. * h
    # zom : roughness length for momentum transfer [m]
    zom = 0.123 * h
    # zoh : roughness length for heat and vaopr transfer [m]
    zoh = 0.0123 * h
    # k: von Karman's constant
    k = 0.41
    ra = np.log((zw - d) / zom) * np.log((zh - d) / zoh) / k ** 2 / uz
    return ra


#### --------------------------------------------------------------------------


class PotentialEvapoTranspiration(SolarInformation):
    """
    Computes the Potential EvapoTranspiration (PET) through different models.
    
    The class instance stores the information (temperatures series, solar
    radiations, altitude...) required to calculate the PET for one given site
    over several dates with different models.

    The PET estimated with a given model is the output of the instance method
    corresponding to the model.
    Thus, it is very easy to compare PETs estimated with different models.

    Parameters
    ----------
    tmin : TimeSeries
        Minimum air temperatures [\ :math:`\\textrm{\\textcelsius}`].
    tmax : TimeSeries
        Maximum air temperatures [\ :math:`\\textrm{\\textcelsius}`].
    latitude : float
        Latitude of the site (in degrees).
    z : float
        Altitude of the site [m].
    %(Kt)s
    """
    __doc__ = (__doc__ or '') % _doc_parameters
    #
    def __init__(self, tmin, tmax, latitude,
                 z=0, Kt=0.170):
        (tmin, tmax) = _checktimeseries(tmin, tmax)
        (self.tmin, self.tmax) = (tmin, tmax)
        self.tavg = (self.tmin + self.tmax) / 2.
        SolarInformation.__init__(self, tmin._dates, latitude)
        self._z = z
        self._Kt = None
        self.Kt = Kt
        #
        self._lambda = None
        self.lambdah = 2.451
        self.atm_pressure = 101.3
        self._es = None
        self._Delta = None


    def _get_altitude(self):
        """Altitude of the site [m]."""
        return self._z
    def _set_altitude(self, z):
        "(sets the attribute)"
        self._z = z
    altitude = property(fget=_get_altitude, fset=_set_altitude,
                        doc=_get_altitude.__doc__)

    def _get_Kt(self):
        """
    Empirical coefficient for the estimation of solar radiations from
    temperatures.
    A value of 0.171 is used by default.

    If set to ``None``, :math:`K_t` is computed from :math:`\Delta T`, 
    the temperature range, as
    :math:`K_t = 0.00185 {\Delta T}^2 - 0.0433 {\Delta T} + 0.4023`.  
    """
        return self._Kt
    def _set_Kt(self, value):
        "(sets Kt)"
        if value is None:
            trange = self.tmax - self.tmin
            value = 0.00185 * trange ** 2 - 0.0433 * trange + 0.4023
        self._Kt = value
    Kt = property(fget=_get_Kt, fset=_set_Kt, doc=_get_Kt.__doc__)

    def _get_lambda(self):
        """
    Latent heat of vaporization [MJ.kg\ :sup:`-1`].
    A default of 2.451 is assumed. If set to ``None``, ``lambdah`` is
    computed from the average air temperature :math:`T_{avg}` as
    :math:`\lambda_H = 2.501 - 0.0002361  T_{avg}`
        """
        return self._lambda
    def _set_lambda(self, value):
        "(sets the attribute)"
        if value is None:
            tair = (self.tmax + self.tmin) / 2.
            self._lambda = 2.501 - (2.361e-3) * tair
        else:
            self._lambda = 2.451
        return
    lambdah = property(fget=_get_lambda, fset=_set_lambda,
                       doc=_get_lambda.__doc__)


    def _get_Delta(self):
        """
    Slope of the vapor pressure-temperature relationship, as

    .. math::
       \Delta = 2502.992210 \cdot \\exp{\left[17.2693882\cfrac{T_{avg}}{T_{avg}+237.3}\\right]} \cdot (T_{avg}+237.3)^{-2}
    
    where :math:`T_{avg}` is average air temperature 
    [:math:`\\textrm{\\textcelsius}`].

        """
        if self._Delta is None:
            _t = (self.tavg + 237.3)
            self._Delta = 2502.992210 * ma.exp(17.2693882 * self.tavg / _t) / _t ** 2
        return self._Delta
    Delta = property(fget=_get_Delta, doc=_get_Delta.__doc__)

    #.................................................
    def _get_satvap(self, tair):
        return 0.61078 * ma.exp(17.2693882 * tair / (tair + 237.3))
    #
    def _get_saturated_vapor_pressure(self):
        """
    Returns the vapor pressure at saturation [kPa].
    The saturated vapor pressure for the time step is calculated as the average
    of the saturated vapor pressures corresponding to the minimum and maximum
    temperatures for the time step.
    
    The relationship between saturated vapor pressure :math:`e_s^*` and
    temperature :math:`T` is:
    
    .. math::
       e_s^*(T) = 0.61078 \cdot \\exp\left[\\frac{17.2693882 T}{T+237.3} \\right]
    
        """
        if self._es is None:
            es = self._get_satvap
            self._es = (es(self.tmin) + es(self.tmax)) / 2.
        return self._es
    saturated_vapor_pressure = property(fget=_get_saturated_vapor_pressure,
                                    doc=_get_saturated_vapor_pressure.__doc__)
    #
    def vapor_pressure_from_humidity(self, RHmin=None, RHmax=None, RHmean=None):
        """
    Returns the vapor pressure :math:`e_a` [kPa] calculated from relative
    humidity.
    
    The relative humidity is the ratio of the vapor pressure :math:`e_a` on the
    vapor pressure at saturation :math:`e_s^*` at the same temperature.
    
    If maximum and minimum relative humidities (:math:`RH_{max}` and 
    :math:`RH_{min}`) are available, the vapor pressure is calculated as
    
    .. math::
       e_a = \left[ e_s^*(T_{min}) \cdot RH_{max} + e_s^*(T_{max}) \cdot RH_{min}\\right] / 2
    
    If :math:`RH_{min}` is unavailable, the previous equation simplifies in
    
    .. math::
       e_a = e_s^*(T_{min}) \cdot RH_{max}
    
    If only the average relative humidity :math:`RH_{mean}` is available, 
    :math:`e_a` is estimated with:
    
    .. math::
       e_a = RH_{mean} \cdot \left[ e_s^*(T_{min}) + e_s^*(T_{max})\\right] / 2
    
    
    Parameters
    ----------
    RHmin : None, TimeSeries
       Minimum relative humidity for the time step.
    RHmax : None, TimeSeries
       Maximum relative humidity for the time step.
    RHmean : None, TimeSeries
       Mean relative humidity for the time step.
        """
        es = self._get_satvap
        if RHmin is None:
            if RHmax is None:
                if RHmean is None:
                    raise ValueError("No humidity data !")
                else:
                    return RHmean * self.saturated_vapor_pressure
            return RHmax * es(self.tmin)
        elif RHmax is None:
            return RHmin * es(self.tmax)
        ea = (RHmax * es(self.tmin) + RHmin * es(self.tmax)) / 2.
        return ea
    # 
    def vapor_pressure_from_dewpoint(self, tdew=None):
        """
    Returns the vapor pressure [kPa] calculated from the dew point temperature,
    :math:`T_{dew}`:
    
    .. math::
       e_a = e_s^*(T_{dew})
    
    If the dew point temperature is None, it is approximated by the minimum air
    temperature :math:`T_{min}`.
    
    Parameters
    ----------
    tdew : None, TimeSeries
       Dew point temperature [:math:`\\textrm{\\textcelsius}`].
       If None, use the minimum temperature.
       """
        if tdew is None:
            tdew = self.tmin
        else:
            tdew = force_reference(tdew, self.tmin)
        return self._get_satvap(tdew)


    #.................................................
    def set_atmospheric_pressure(self, P, zo=0, tmean=0):
        """
    Sets the atmospheric pressure.
    
    If ``P`` is None, it is calculated from :math:`P_o` [101.3 kPa],
    the average atmospheric pressure at sea level.

    Parameters
    ----------
    P : TimeSeries
        Series of atmospheric pressures [kPa].
    z0 : {0., float} optional
        Elevation at reference level [m].
    tmean : {20., float} optional
        Mean temperature at the reference level
        [:math:`\\textrm\{\\textcelsius\}`].

    References
    ----------
    Burman et al., 1987
        """
        if P is None:
            # a1: Constant lapse moist air [K.m^{-1}] : 0.0065
            # R: Specific gas constant [J.kg^{-1}.K^{-1}]: 286.9
            # g: Gravitational acceleration [m.s^{-2}] : 9.087
            # Po: Atmospheric pressure at sea level [kPa]
            (a1, R, g, Po) = (0.0065, 286.9, 9.087, 101.3)
            #....
            tmean = self.tavg + 273.16
            z = self.altitude
            self.atm_pressure = Po * (1. - (a1 * (z - zo)) / tmean) ** (g / (a1 * R))
        else:
            self.atm_pressure = force_reference(P, self.tavg)
        return


    def _get_gamma(self):
        """
    Psychrometric constant 
    [kPa/\ :math:`\\textrm{\\textcelsius}`], as given by 
    :math:`\gamma = \\frac{c_p}{\\epsilon} \\frac{P}{\lambda_h}`
    where :math:`c_p = 1.013 10^{-3}` 
    [MJ.kg\ :sup:`-1`.\ :math:`\\textrm{\\textcelsius}`\ :sup:`-1`] is the 
    specific heat of water at constant pressure and :math:`\epsilon = 0.622`
    is the ratio molecular weight of water vapor on dry air.
        """
        # cp : Specific heat of moist air [MJ.kg^{-1}.oC^{-1}]: 1.013e-3
        # epsilon : Ratio molecular weight water vapor/dry air: 0.622
        return (1.013e-3 / 0.622) * self.atm_pressure / self.lambdah
    gamma = property(fget=_get_gamma,
                     doc=_get_gamma.__doc__)


    def apx_solar_radiations(self, Kt=None):
        """
    Returns the solar radiations [MJ.m\ :sup:`-2`.d\ :sup:`-1`] approximated
    from temperatures, as:
    
    .. math::
       R_s = K_t \cdot \sqrt{T_{max}-T_{min}} \cdot R_a

    where :math:`K_t` is an empirical parameter.
    By default, :math:`K_t=0.171`. Otherwise, :math:`K_t` is estimated from the
    temperature range :math:`\Delta T=T_{max}-T_{min}` as
    :math:`K_t = 0.00185 (\Delta T)^2 - 0.0433 \Delta T + 0.4023`.
    
    Parameters
    ----------
    %(Kt)s
        """
        if Kt is None:
            Kt = self._Kt
        return self.solar_radiations_from_temperatures(self.tmin, self.tmax, Kt)
    apx_solar_radiations.__doc__ = (apx_solar_radiations.__doc__ or '') % \
                                   _doc_parameters

    def net_lw_radiations(self, ea, solar_radiations):
        """
    Computes the net longwave radiation [MJ.m\ :sup:`-2`.d\ :sup:`-1`] as a
    function of temperatures:
    
    .. math::
       R_{nl} = \sigma \cdot \left[\\cfrac{T_{min}^4 + T_{max}^4}{2}\\right] \cdot
                (0.34 - 0.139\sqrt{(e_a)} \cdot
                \left[ 1.35 \cfrac{R_s}{R_{so}} - 0.35\\right]

    where :math:`\sigma = 4.903 10^{-9}` [MJ.K\ :sup:`-4`.m\ :sup:`-2`.d\ :sup:`-1`]
    is the Stefan-Bolztmann constant,
    :math:`e_a` the actual vapor pressure,
    :math:`R_s` the actual incoming solar radiation,
    :math:`R_{so}` the incoming clear-sky solar radiation.
    The temperatures are expressed in [K]

    Parameters
    ----------
    ea : float
        Actual vapor pressure [kPa]
    solar_radiations : TimeSeries
        Mean incoming solar radiations :math:`R_s` [MJ.m\ :sup:`-2`.d\ :sup:`-1`]

        """
        (tmin, tmax) = (self.tmin, self.tmax)
        # sigma : Stefan-Boltzmann constant 4.903e-9 MJ.K^{-4}.m^{-2}.d^{-1}
        # Transforms the temperatures to absolute
        tKmin = tmin + 273.16
        tKmax = tmax + 273.16
        tKfactor = (tKmin ** 4 + tKmax ** 4) / 2.
        #
        Rs = solar_radiations
        Rso = self.clearsky_solar_radiations(self._z)
        # Make sure that the relative shortwave radiation ratio is in [0.2, 1.]
        rel_sw = np.clip(Rs / Rso, 0.2, 1.0)
        # Cloudiness factor
        rsw_factor = (1.35 * rel_sw - 0.35)
        return 4.903e-9 * tKfactor * (0.34 - 0.139 * np.sqrt(ea)) * rsw_factor


    def _get_soil_heat_flux_density(self):
        """
    Computes the soil heat flux density.
    The soil heat flux density is set to zero when the time step is less
    a month.
    Otherwise, it is computed as
    :math:`G_{i} = 0.07 \left( t_{air,i+1} - t_{air,i-1} \\right)`
        """
        tair = self.tavg
        G = hydro.empty_like(tair)
        G.flat = 0
        if tair.freq <= _c.FR_MTH and tair.size > 3:
            G[1:-1] = 0.07 * (tair[2:] - tair[:-2])
        return G
    soil_heat_flux_density = property(fget=_get_soil_heat_flux_density,
                                      doc=_get_soil_heat_flux_density.__doc__)


    def PenmanMonteith(self,
                       Rs=None, sunshinehours=None, a_s=0.25, b_s=0.5, Kt=None,
                       tdew=None, RHmin=None, RHmax=None, RHmean=None,
                       crop_ref='short', u2=2, albedo=0.23):
        """
    Returns the potential evapotranspiration estimated by the 
    Penman-Monteith equation, using the approximation of the ASCE.

    The Penman-Monteith equation is

    .. math::
       \lambda_H ET = \cfrac{\Delta (R_n-G) + \\rho_a c_p (e_s-e_a)/r_a}{
                             \Delta + \gamma \left(1 + \\frac{r_s}{r_a}\\right)}
    
    where 
    :math:`\lambda` is the
    :math:`R_n = R_{ns}-R_{nl}` the net radiation as the difference between the
    incoming shortwave radiation :math:`R_{ns}` and incoming longwave radiation
    R_{nl};
    :math:`G` the soil flux density;
    :math:`\\rho_a` the mean air density at constant pressure;
    :math:`c_p` the specif heat of air;
    :math:`e_s-e_a` the vapor pressure deficit of air;
    :math:`r_s` and :math:`r_a` the bulk surface and aerodynamic resistances.

    The net incoming shortwave radiation :math:`R_{ns}` is calculated from the 
    incoming solar radiations :math:`R_s` as :math:`R_{ns} = (1-\\alpha) R_s`
    where :math:`\\alpha` is the albedo.
    
    The aerodynamic resistance is calculated assuming a standardized measurement
    height for wind and humidity of 2m: 

    .. math::
     r_a = \cfrac{\ln{\left[\cfrac{2-2/3h}{0.123h}\\right]} \ln{\left[\cfrac{2-2/3h}{0.0123h}\\right]}}{(0.41)^2 u_2}.

    where :math:`h` is the crop height.
    For small crops (clipped grass), a value :math:`h=0.12` is assumed, with a
    bulk surface resistance :math:`r_s=70` [m.s\ :sup:`-1`].
    For tall crops (alfalfa), a value :math:`h=0.5` is assumed, with a
    bulk surface resistance :math:`r_s=45` [m.s\ :sup:`-1`].

    

    Parameters
    ----------
    %(Rs)s
    sunshineshours : TimeSeries, optional
        Series of actual sunshine hours.
    a_s : float, optional
        Regression constant.
    b_s : float, optional
        Regression constant.
    %(Kt)s
    %(crop_ref)s 
    %(z)s
    %(u2)s
    %(albedo)s
    %(RHmin)s
    %(RHmax)s
    %(RHmean)s

        """
        tair = self.tavg
        (Delta, gamma, lambdah) = (self.Delta, self.gamma, self.lambdah)
        # es, ea : mean saturation vapor pressure, actual vapor pressure
        es = self.saturated_vapor_pressure
        if (RHmin, RHmax, RHmean) != (None, None, None):
            ea = self.vapor_pressure_from_humidity(RHmin, RHmax, RHmean)
        else:
            ea = self.vapor_pressure_from_dewpoint(tdew)
        # Rns/Rnl : Net shortwave/longwave radiations [MJ.kg^{-1}.d^{-1}]
        if Rs is None:
            if sunshinehours is None:
                Rs = self.apx_solar_radiations(Kt=Kt)
            else:
                Rs = self.solar_radiations_from_sunshine(sunshinehours,
                                                         a_s=a_s, b_s=b_s)
        # Get the radiation term ($\Delta [R_{ns}-R_{nl}-G] / \lambda$)
        Rns = (1. - albedo) * Rs
        Rnl = self.net_lw_radiations(ea, Rs)
        G = self.soil_heat_flux_density
        radiation_term = Delta * (Rns - Rnl - G) / lambdah
        # Get the aerodynamic resistance [s/m]
        if crop_ref == 'tall':
            (h, rs) = (0.5, 45)
        elif crop_ref == 'short':
            (h, rs) = (0.12, 70)
        (zom, zoh) = (0.123 * h, 0.0123 * h)
        ra = np.log((2 - 2 / 3. * h) / zom) * np.log((2 - 2 / 3. * h) / zoh) / (0.41) ** 2 / u2
        # Get the aerodynamic term#
        P = self.atm_pressure
        rhoa = P / (0.287 * (tair + 273.16)) * (1. - 0.378 * ea / P)
        cp = 1.013e-3
        aerodynamic_term = (rhoa * cp) / lambdah * (es - ea) * 86400. / ra
        #
        num = radiation_term + aerodynamic_term
        denom = Delta + gamma * (1. + rs / ra)
        return num / denom
    PenmanMonteith.__doc__ = (PenmanMonteith.__doc__ or '') % _doc_parameters

    def PenmanMonteithASCE(self,
                           Rs=None, sunshinehours=None, a_s=0.25, b_s=0.5, Kt=None,
                           tdew=None, RHmin=None, RHmax=None, RHmean=None,
                           crop_ref='short', u2=2, albedo=0.23):
        """
    Returns the potential evapotranspiration estimated by the 
    Penman-Monteith equation, using the approximation of the ASCE [Allen_etal_1998]._

    .. math::
       \lambda_H ET_o = \cfrac{\Delta (R_n-G) + \gamma \cfrac{C_n}{T_{avg}+273} u_2 (e_s^* - e_a)}{
                           \Delta + \gamma \left(1 + C_d u_2\\right)}
    
    where 
    :math:`u_2` is the wind speed measured at an altitude of 2m [m.s\ :sup:`-1`],
    :math:`C_p` and :math:`C_d` two constants depending on the reference crop
    and the time interval: for daily to monthly time steps, 
    :math:`C_n=900` and :math:`C_d=0.34` for clipped grass, 
    and :math:`C_n=1600` and :math:`C_d=0.38` for alfalfa.
    
    The values of these coefficients were determined by assuming a bulk
    resistance of 70 [s.m\ :sup:`-1`] for clipped grass of an height of 0.12 m, 
    and of 45 [s.m\ :sup:`-1`] for alfalfa at a height of 0.50~m.
    A constant albedo :math:`\\alpha \\approx 0.23` is also assumed.



    Parameters
    ----------
    %(Rs)s
    sunshineshours : TimeSeries, optional
        Series of actual sunshine hours.
    a_s : float
        Regression constant.
    b_s : float
        Regression constant.
    %(Kt)s
    %(crop_ref)s
    %(z)s
    %(u2)s
    %(albedo)s
    %(RHmin)s
    %(RHmax)s
    %(RHmean)s

        """
        tair = self.tavg
        # Crop constants
        # For clipped grass, Cn=897.211450 and Cd=0.33708287
        constants = {'short':{'Cn':900, 'Cd':0.34},
                     'tall' :{'Cn':1600, 'Cd':0.38}}
        Cn = constants[crop_ref]['Cn']
        Cd = constants[crop_ref]['Cd']
        #..........
        (Delta, gamma, lambdah) = (self.Delta, self.gamma, self.lambdah)
        # es, ea : mean saturation vapor pressure, actual vapor pressure
        es = self.saturated_vapor_pressure
        if (RHmin, RHmax, RHmean) != (None, None, None):
            ea = self.vapor_pressure_from_humidity(RHmin, RHmax, RHmean)
        else:
            ea = self.vapor_pressure_from_dewpoint(tdew)
        # Rns/Rnl : Net shortwave/longwave radiations [MJ.kg^{-1}.d^{-1}]
        if Rs is None:
            if sunshinehours is None:
                Rs = self.apx_solar_radiations(Kt=Kt)
            else:
                Rs = self.solar_radiations_from_sunshine(sunshinehours,
                                                         a_s=a_s, b_s=b_s)
        Rns = (1. - albedo) * Rs
        Rnl = self.net_lw_radiations(ea, Rs)
        G = self.soil_heat_flux_density
        radiation_term = Delta * (Rns - Rnl - G) / lambdah
        aerodynamic_term = gamma * Cn / (tair + 273) * u2 * (es - ea)
        #
        num = radiation_term + aerodynamic_term
        denom = Delta + gamma * (1. + Cd * u2)
        return num / denom
    PenmanMonteithASCE.__doc__ = (PenmanMonteithASCE.__doc__ or '') % _doc_parameters


    def PriestleyTaylor(self,
                        Rs=None, albedo=0.23,
                        tdew=None, RHmin=None, RHmax=None, RHmean=None):
        """
    Computes the PET with the Priestley-Talyor method [Priestley_Taylor_1972]_.

    .. math::
       \lambda_H ET_o = a \cfrac{\Delta}{\Delta+\gamma}
                         [(1-\\alpha)R_s - R_{nl} - G]

    where :math:`a` is an empirical coefficient set to 1.26 in humid regions.


    Parameters
    ----------
    %(Rs)s
    %(albedo)s
    %(RHmin)s
    %(RHmax)s
    %(RHmean)s

        """
        (Delta, gamma, lambdah) = (self.Delta, self.gamma, self.lambdah)
        # ea : mean saturation vapor pressure, actual vapor pressure
        if (RHmin, RHmax, RHmean) != (None, None, None):
            ea = self.vapor_pressure_from_humidity(RHmin, RHmax, RHmean)
        else:
            ea = self.vapor_pressure_from_dewpoint(tdew)
        # Solar radiation
        if Rs is None:
            Rs = self.apx_solar_radiations(Kt=self.Kt)
        Rns = (1. - albedo) * Rs
        Rnl = self.net_lw_radiations(ea, Rs)
        Rn = Rns - Rnl
        G = self.soil_heat_flux_density
        return 1.26 * Delta / (Delta + gamma) * (Rn - G) / lambdah
    PriestleyTaylor.__doc__ = (PriestleyTaylor.__doc__ or '') % _doc_parameters


    def Makkink(self, Rs=None):
        """
    Computes the PET with the Makkink (1957) method [Makkink_1957]_:

    .. math::
       PET = 0.61 \cfrac{\Delta}{\Delta+\gamma} \cfrac{R_s}{\lambda_H} - 0.12

    Parameters
    ----------
    %(Rs)s

        """
        (Delta, gamma, lambdah) = (self.Delta, self.gamma, self.lambdah)
        # Rs : solar radiation
        if Rs is None:
            Rs = self.apx_solar_radiations(Kt=self.Kt)
        else:
            Rs = force_reference(Rs, self.tavg)
        PET = 0.61 * Delta / (Delta + gamma) * Rs / lambdah - 0.12
        ma.putmask(PET, (PET._series.filled(0) < 0), 0)
        return PET
    Makkink.__doc__ = (Makkink.__doc__ or '') % _doc_parameters


    def Hansen(self, Rs=None):
        """
    Computes the PET with the Hansen (1957) method [Hansen_1957]_.

    .. math::
       PET = 0.70 \cfrac{\Delta}{\Delta+\gamma} \cfrac{R_s}{\lambda_H}

    Parameters
    ----------
    %(Rs)s

        """
        (Delta, gamma, lambdah) = (self.Delta, self.gamma, self.lambdah)
        # Rs : solar radiation
        if Rs is None:
            Rs = self.apx_solar_radiations(Kt=self.Kt)
        else:
            Rs = force_reference(Rs, self.tavg)
        return 0.70 * Delta / (Delta + gamma) * Rs / lambdah
    Hansen.__doc__ = (Hansen.__doc__ or '') % _doc_parameters



    def Turc(self, Rs=None, RHmean=0.5):
        """
    Computes the PET with the Turc (1961) method [Turc_1961]_:

    .. math::
       PET = a \; 0.013 \\frac{T_{avg}}{T_{avg}+15} \\frac{23.8856 R_s + 50}{\lambda}
    
    where :math:`a` is a coefficient dependent on the relative humidity :math:`RH`:
    
    .. math::

       \\begin{aligned}
       a &= 1.                &\\qquad RH \geq 0.5 \\\\
         &= 1. + (0.5-RH)/0.7 &\\qquad RH < 0.5
       \end{aligned} 

    This model should not be applied when the mean air temperature is less than 
    :math:`-10\\textcelsius`: the results will be masked when :attr:tavg < -10.

    Parameters
    ----------
    %(Rs)s
    %(RHmean)s

        """
        # Mask the temperatures less than -10oC, as the formula is not valid
#        tair = ma.masked_less(self.tavg, -10.)
        # In fact, let's mask negative temperatures...
        tair = ma.masked_less(self.tavg, 0.)
        # Rs : solar radiation
        if Rs is None:
            Rs = self.apx_solar_radiations(Kt=self.Kt)
        else:
            Rs = force_reference(Rs, tair)
        #
        at = ma.where(RHmean >= 0.5, 1., 1. + (0.5 - RHmean) / 0.70)
        return at * 0.013 * tair / (tair + 15) * (23.8856 * Rs + 50) / self.lambdah
    Turc.__doc__ = (Turc.__doc__ or '') % _doc_parameters


    def Hargreaves(self, Kt=None):
        """
    Computes the PET [mm.d\ :sup:`-1`] with the 
    method of Hargreaves and Samani [Hargreaves_Samani_1985]_.

    .. math::
       PET = 0.0135 * (R_s/\lambda) * (T_{avg} + 17.8)
    
    where the solar radiations :math:`R_s` are estimated from the extraterrestrial
    solar radiations as :math:`R_s = K_t \sqrt{(T_{max}-T_{min})} R_a`.

    Parameters
    ----------
    %(Kt)s

    Returns
    -------
    %(result)s

        """
        if Kt is None:
            Kt = self._Kt
        tair = self.tavg
        lambdah = self.lambdah
        Rs = self.apx_solar_radiations(Kt=Kt)
        PET_ = 0.0135 * Rs / lambdah * (tair + 17.8)
        return PET_
    Hargreaves.__doc__ = (Hargreaves.__doc__ or '') % _doc_parameters


    def Droogers(self, rainfall):
        """
    Computes the PET with the method of Hargreaves and Samani
    [Hargreaves_Samani_1985]_, modified by Droogers and Allen
    [Droogers_Allen_2002]_.
    
    .. math::
       ET_o = 0.0013 (R_a/\lambda_H) (T_{avg}+17) (\Delta T - 0.0123 P_{mm})^{0.76}
    
    where :math:`P_{mm}` is the cumulative rainfall [mm].
    
    Parameters
    ----------
    rainfall : TimeSeries
        Monthly cumulative rainfall [mm],
        as a :class:`~scikits.timeseries.TimeSeries` object.
    
    Notes
    -----
    
    This relation was derived for monthly data, assuming that total precipitation
    is to some extent representative of relative humidity.
    A warning will be issued if the frequency of the
    :class:`PotentialEvapoTranspiration` is not monthly.
    
    [Droogers_Allen_2002]_ compared their estimate of :math:`ET_o` 
    with the estimates calculated with the ASCE version of Penman-Monteith  and
    the Hargreaves model [Hargreaves_Samani_1985]_ to assess the effiency of each
    model in the presence of inaccurate climatic data.
    They concluded that their model was better at approximating the Penman-Monteith
    ASCE equation, and that it was more accurate than this latter when random
    errors were introduced in the dataset.
    
        """
        if not (np.isscalar(rainfall) or isinstance(rainfall, TimeSeries)):
            err_msg = "Precipitation input should be a valid TimeSeries object!"
            raise TypeError(err_msg)
        if self.freq != _c.FR_MTH:
            warnings.warn("Droogers relation valid for monthly data only!",
                          UserWarning)
        tair = self.tavg
        trange = self.tmax - self.tmin
        rainfall = force_reference(rainfall, self.tmin)
        etrad = self.extraterrestrial_solar_radiations / self.lambdah
        trange_corrected = np.clip((trange - 0.0123 * rainfall), 0, 1e+20)
        PET_ = 0.0013 * etrad * (tair + 17) * trange_corrected ** 0.76
        return PET_
    Droogers.__doc__ = (Droogers.__doc__ or '') % _doc_parameters


#    #--- Temperature based -- --------------------------------------------

    def Hamon(self):
        """
    Computes the PET [mm.d\ :sup:`-1`] with the method of Hamon (1961)
    [Hamon_1961]_.

    .. math::
       PET = a \left(\\frac{D}{12}\\right)^2 \\frac{e_{s}^{*}(T_{avg})}{R (T_{avg} + 273.2)}
      
    where :math:`D` is the number of daylight hours [hr],
    :math:`T_{avg}` the average air temperature :math:`[\\textcelsius]`
    :math:`e_s^{*}` the saturated vapor pressure at temperature :math:`T_{avg}` [kPa]
    :math:`R=0.461495` the ideal gas constant [kJ.kg\ :sup:`-1`.K\ :sup:`-1`]
    and :math:`a = 0.14` a conversion constant

        """
        # Masks the dates for which the temperature is negative
        # We need to make a copy to avoid propagation of modifs to the mask
        tair = ma.masked_less(self.tavg, 0)
        # Get the daylengths
        daylengths = self.daylength / 12.
        # Saturated water pressure (kPa)
        swp = self._get_satvap(tair)
        # Ideal gas constant for water vapor [kJ/kg/K]
        Rw = 0.461495
        # Absolute humidity [kg/m3] : 
        absh = 1000. * swp / (Rw * (tair + 273.15))
        # Constant : 0.0055 [in/d] / 12**2 -> mm/d
        chih = 0.0055 * 25.4
        # PET Hamon (mm/d)
        return chih * daylengths ** 2 * absh
    Hamon.__doc__ = (Hamon.__doc__ or '') % _doc_parameters


    def Kharrufa(self):
        """
    Computes the PET with Kharrufa's approximation [Kharrufa_1985]_.

    .. math::
       PET = 0.34 * D_{\%} * T^{1.34}_{avg}
    
    where :math:`D_{\%}` is the percentage of total daylight hours per period
    with reference to the total annual daylight hours, 
    and :math:`T_{avg}` the mean temperature [:math:`\textcelsius`].

        """
        tair = self.tavg
        freq = self.freq
        # get the fractional daylight as ratio daylight(d)/daylight(y)
        daylight = hydro.time_series(self._daylighthours, dates=self._dates)
        anndaylight = daylight.convert('A', ma.sum)
        anndaylight = backward_fill(anndaylight.convert('D', position='END'))
        if len(anndaylight) > len(daylight):
            anndaylight = hydro.adjust_endpoints(anndaylight,
                                                 start_date=daylight.start_date,
                                                 end_date=daylight.end_date)
        fracdaylight = daylight / anndaylight
        if freq <= _c.FR_DAY:
            fracdaylight = fracdaylight.convert(freq, ma.sum)
        return 0.34 * fracdaylight * ma.where(tair > 0, tair ** 1.3, 0)


    def Thornthwaite(self):
        """
    Computes the potential evapotranspiration [mm.d\ :sup:`-1`] with the 
    method of Thornthwaite.
    
    .. math::
       PET = 16 \left(\\frac{10 T_{avg}}{I}\right)^{a} \\frac{D}{12} \\frac{N}{30}
    
    where :math:`D` is the number of daylight hours [hr],
    :math:`T_{avg}` the average air temperature [:math:`\textcelsius`]
    
    .. math::
        I = \sum_{i=1}^{12}{(\bar{T}_{avg, i}/5)^{1.514}}
    
    where :math:`\bar{T}_{avg,i}` is the mean air temperature for month `i`,
    and
    
    .. math::
       a = 0.49239 + 0.0179 I - 77.1e^{-6} I^{2} + 67.5e^{-8} I^{3}

        """
        tavg = self.tavg
        freq = self.freq
        if freq > _c.FR_MTH:
            tair = tavg.convert('M', func=ma.mean)
            daylighthours = self.daylength.convert('M', func=ma.mean)
        else:
            tair = tavg
            daylighthours = self.daylength
        # Get the day-per-month correction using a 30d month as a reference
        dpmcorr = time_series(tair.days / 30., dates=tair.dates)
        # Mask negative temperatures
        tair = ma.masked_less(tair, 0)
        # Get the average monthly temperatures
        tairm = tair.convert('A').mean(0)
        # Get the temperature indices
        I = ((tairm / 5.) ** 1.514).sum()
        a = 0.49239 + 0.0179 * I - 77.1e-6 * I ** 2 + 67.5e-8 * I ** 3
        # Compute the PET
        PET = (10 * tair / I) ** a * 16. * daylighthours / 12. * dpmcorr
        # Take the high temperatures into account
        ma.putmask(PET, tair >= 26.5, -415.85 + 32.25 * tair - 0.43 * tair ** 2)
        PET /= tair.days
        if freq >= _c.FR_WK:
            _PET = backward_fill(PET.convert(freq, position='END'))
            PET = hydro.align_with(tavg, _PET)
        return PET
    Thornthwaite.__doc__ = (Thornthwaite.__doc__ or '') % _doc_parameters

PETModel = PotentialEvapoTranspiration




class PETEvaluation(object):
    """
    A holdall class for the comparison of PET estimates.

    The PET models listed in the :keyword:`models` are all evaluated from
    the series of input temperatures :keyword:`tmin` and :keyword:`tmax`, for
    the given :keyword:`latitude`.

    Parameters
    ----------
    tmin : TimeSeries
        Series of minimum air temperatures [\ :math:`\\textrm{\\textcelsius}`].
    tmax : TimeSeries
        Series of maximum air temperatures [\ :math:`\\textrm{\\textcelsius}`].
    latitude : float
        Latitude of the site (in degrees).
    rain : TimeSeries, optional
        Series of precipitation [mm].
    %(Kt)s
    %(models)s
    %(reference)s

    """
    __doc__ = (__doc__ or '') % _doc_parameters
    #
    _defaultmodels = ['Droogers',
                      'Hamon',
                      'Hansen',
                      'Hargreaves',
                      'Kharrufa',
                      'Makkink',
                      'PenmanMonteith',
                      'PenmanMonteithASCE',
                      'PriestleyTaylor',
                      'Thornthwaite',
                      'Turc',
                      ]


    def __init__(self, tmin, tmax, latitude, rain=None, Kt=0.17,
                 models=None, reference=None):
        # Check the reference
        if reference is None:
            self.reference = 'PenmanMonteith'
        elif reference in self._defaultmodels:
            self.reference = reference
        else:
            errmsg = "Invalid reference '%s'.\nPlease choose in %s)"
            raise ValueError(errmsg % (reference, self._defaultmodels))
        # Check the models and define the default output
        (self.models, self.ndtype) = self._checkmodels(models,
                                                       self._defaultmodels)
        # Initialize the PET class
        _PET = PETModel(tmin, tmax, latitude, Kt=Kt)
        # Construct a dictionary {model : PET estimate}
        PET = dict([(m, getattr(_PET, m).__call__()) for m in self.models[1:]])
        if 'Droogers' in self.models:
            if rain is None:
                PET['Droogers'] = ma.array(np.empty_like(PET[self.models[1]]),
                                           mask=True)
            else:
                PET['Droogers'] = _PET.Droogers(rain)
        # Construct a temporary ClimateRecords array
        pet = ma.empty(len(tmin), dtype=self.ndtype).view(hydro.TimeSeries)
        pet._dates = tmin.dates
        for (m, p) in sorted(PET.iteritems()):
            setattr(self, m, p)
            pet[m] = p
        # Sets some basic attributes
        self._PET = _PET
        self.estimates = pet

    def _checkreference(self, reference):
        """Checks the reference PET."""
        models = self.models
        if reference is None:
            reference = self.reference
        if isinstance(reference, basestring):
            if not (reference in models):
                raise ValueError("The reference should be in %s ! "\
                                 "(got ['%s'] instead)" % (models, reference))
            reference = self.estimates[reference]
        return reference

    def _checkmodels(self, models, default):
        "Check the validity of models"
        if models is None:
            return (default, np.dtype([(m, float) for m in default]))
        if isinstance(models, basestring):
            models = [models, ]
        checked = [mod for mod in models if mod in default]
        if not len(checked):
            if not len(self.models):
                rmsg = "No valid models were specified.\n Please choose in %s."
                raise ValueError(rmsg % self._defaultmodels)
        ndtype = [(m, float) for m in checked]
        return (checked, ndtype)

    def _select_estimates(self, models):
        """
        Private function: selects the estimates from a list of models.
        """
        #
        estimates = self.estimates
        (models, ndtype) = self._checkmodels(models, self.models)
        if list(models) == self.models:
            return estimates
        pet = ma.empty(len(estimates), ndtype).view(hydro.TimeSeries)
        for m in models:
            pet[m] = estimates[m]
        pet._dates = estimates._dates
        return pet


    def RMSE(self, reference=None, models=None):
        """
    Returns the Root Mean Square Error for the PET estimates.

    Parameters
    ----------
    %(reference)s
    %(models)s

    Returns
    -------
    RMSE : MaskedArray
        A structured :class:`~numpy.ma.MaskedArray`, where each field corresponds
        to the RMSE between a model and the reference
        """
        # Check the reference
        reference = self._checkreference(reference)
        # Select the estimates from the list of models
        pet = self._select_estimates(models)
        mpet = pet.view((float, len(pet.dtype)))
        #
        mpet_rmse = mpet - reference[:, None]
        mpet_rmse **= 2
        RMSE = ma.sqrt(mpet_rmse).mean(0).view(pet.dtype)
        return RMSE
    RMSE.__doc__ = (RMSE.__doc__ or '') % _doc_parameters


    def RMSE_per_period(self, freq=None, reference=None, models=None):
        """
    Returns the Root Mean Square Error for the PET estimates, per period.

    Parameters
    ----------
    %(reference)s
    %(models)s
    %(freq)s

    Returns
    -------
    result : MaskedArray
        A (n,) structured :class:`~numpy.ma.MaskedArray` object, with
        ``n`` the number of periods per year (``n=12`` for monthly, ``n=4`` for
        quarterly...), and where each field correspond to the RMSE between
        the given model and the reference.
        """
        pet = self.estimates
        # Get the reference
        reference = self._checkreference(reference)
        # Select the estimates from the list of models
        pet = self._select_estimates(models)
        #
        if freq is None:
            freq = self._PET.freq
        #
        tmpdict = {}
        for m in models:
            _tmp = (pet[m] - reference).convert(freq, func=ma.mean).convert('A')
            _tmp **= 2
            tmpdict[m] = ma.sqrt(_tmp.mean(0))
        result = ma.empty((len(tmpdict[models[0]]),), dtype=pet.dtype)
        for m in models:
            result[m] = tmpdict[m]
        return result
    RMSE_per_period.__doc__ = (RMSE_per_period.__doc__ or '') % _doc_parameters


    def monthly_average(self):
        """
    Returns the monthly averages for each model
        """
        pet = self.estimates
        mavg = ma.empty(12, dtype=self.ndtype)
        for m in self.models:
            mavg[m] = pet[m].convert('A').mean(0)
        return mavg


    def kendalltau(self, freq=None, reference=None, models=None):
        """
    Computes the Kendall tau between each model and the reference,
    for each period.
 
    Parameters
    ----------
    %(reference)s
    %(models)s
    %(freq)s

    Returns
    -------
    result : dictionary
        A dictionary (model <> array of results).
        For each key, the corresponding element is a (n+1,) array of two fields,
        (`tau`, `p-value`), of the Kendall :math:`tau` calculated between the
        given model (key) and the reference and its corresponding `p` value.
        The first n entries correspond to the period (``n=12`` for monthly
        results, ``n=4`` for quaterly results); the last entry corresponds to
        the Spearman rho calculated over the whole range of dates.
        """
        import scipy.stats.mstats as mstats
        reference = self._checkreference(reference)
        pet = self._select_estimates(models)
        #
        if freq is None:
            freq = self._PET.freq
        #
        pet_r = reference.convert(freq, func=ma.mean).convert('A')._series
        result = {}
        for m in pet.dtype.names:
            pet_m = pet[m].convert(freq, func=ma.mean).convert('A')._series
            res = []
            for k in np.arange(pet_m.shape[-1]):
                res.append(mstats.kendalltau(pet_r[:, k], pet_m[:, k],))
            res.append(mstats.kendalltau(pet_r, pet_m))
            result[m] = np.array(res, dtype=[('tau', float), ('pvalue', float)])
        return result
    kendalltau.__doc__ = (kendalltau.__doc__ or '') % _doc_parameters
    #
    def spearmanr(self, freq=None, reference=None, models=None):
        """
    Computes the Spearman rho per period.

    Parameters
    ----------
    %(reference)s
    %(models)s
    %(freq)s

    Returns
    -------
    result : dictionary
        A dictionary (model <> array of results).
        For each key, the corresponding element is a (n+1,) array of two fields,
        (`rho`, `p-value`), of the Spearman :math:`rho` calculated between the
        given model (key) and the reference and its corresponding `p` value.
        The first n entries correspond to the period (``n=12`` for monthly
        results, ``n=4`` for quaterly results); the last entry corresponds to
        the Spearman rho calculated over the whole range of dates.
        """
        import scipy.stats.mstats as mstats
        reference = self._checkreference(reference)
        pet = self._select_estimates(models)
        #
        if freq is None:
            freq = self._PET.freq
        #
        pet_r = reference.convert(freq, func=ma.mean).convert('A')._series
        result = {}
        for m in pet.dtype.names:
            pet_m = pet[m].convert(freq, func=ma.mean).convert('A')._series
            res = []
            for k in np.arange(pet_m.shape[-1]):
                res.append(mstats.spearmanr(pet_r[:, k], pet_m[:, k],))
            res.append(mstats.spearmanr(pet_r, pet_m))
            result[m] = np.array(res, dtype=[('tau', float), ('pvalue', float)])
        return result
    spearmanr.__doc__ = (spearmanr.__doc__ or '') % _doc_parameters

