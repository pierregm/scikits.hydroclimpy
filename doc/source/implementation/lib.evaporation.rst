
.. currentmodule:: scikits.hydroclimpy
.. module:: scikits.hydroclimpy.lib.evaporation

===================================================================
Calculating evapo-transpiration: the :mod:`~lib.evaporation` module
===================================================================

Potential evapotranspiration is an important component of the hydrological
budget.
It corresponds to the maximum amount of water that can be lost from soils 
and vegetation on a given time interval.
Reference ET, often noted :math:`ET_o` is the maximum rate at which water can be
lost from a specific crop, called 'reference crop'.
Two references crop are widely used: clipped grass [Doorenbos_1977]_ 
and full-cover alfalfa [Wright_1982]_.



The :mod:`scikits.hydroclimpy.lib.evaporation` introduces several classes to
calculate potential evapotranspiration.


====================================  ===================================================
:class:`SolarInformation`             Stores information about the 
                                      extraterrestrial solar radiations;
:class:`PotentialEvapoTranspiration`  Stores information required for the 
                                      computation of reference evapotranspiration
:class:`PETEvaluation`                Compares estimates of the PET for a given location.
====================================  ===================================================


.. toctree::
   :maxdepth: 1

   _evaporation.solar_radiations
   _evaporation.potential
   _evaporation.evaluation



Symbols
=======

.. tabularcolumns:: |c|l|c|

+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| Symbol            | Description                                                                        | Units                                                       |
+===================+====================================================================================+=============================================================+
| :math:`\lambda_H` | Latent heat of vaporization.                                                       | [MJ.kg\ :sup:`-1`]                                          |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`\rho_a`    | Mean air density at constant pressure.                                             | [kg.m\ :sup:`-3`]                                           |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`c_p`       | Specific heat of water at constant pressure (:math:`c_p = 1.013 10^{-3}`).         | [MJ.kg\ :sup:`-1`. :math:`\textrm{\textcelsius}`\ :sup:`-1`]|
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`\gamma`    | Psychrometric constant.                                                            | [kPa. :math:`\textrm{\textcelsius}`\ :sup:`-1`]             |
|                   |                                                                                    |                                                             |
|                   | .. math::                                                                          |                                                             |
|                   |    \gamma = \frac{c_p}{\epsilon} \frac{P}{\lambda_h}                               |                                                             |
|                   |                                                                                    |                                                             |
|                   | with :math:`P` the atmospheric pressure [kPa]                                      |                                                             |
|                   | and :math:`\epsilon = 0.622` the ratio molecular weight of water vapor on dry air. |                                                             |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`e_s^{*}`   | Saturated vapor pressure.                                                          | [kPa]                                                       |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`e_a`       | actual vapor pressure at temperature.                                              | [kPa]                                                       |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`\Delta`    | Slope of the vapor pressure - temperature relationship                             | [:math:`\textrm{\textcelsius}`]                             |
|                   |                                                                                    |                                                             |
|                   | .. math::                                                                          |                                                             |
|                   |    \Delta = 2502.992 \cdot \exp{\left[17.269\;T/T_{K}\right]}/T_{K}^2              |                                                             |
|                   |                                                                                    |                                                             |
|                   | with :math:`T_{K} = T + 237.3`                                                     |                                                             |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`r_s`       | (Bulk) surface resistance, that represents the resistance of vapor flow through    |                                                             |
|                   | the evaporating surface.                                                           |                                                             |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`r_a`       | Aerodynamic resistance, that controls the transfer of heat and water vapor from    |                                                             |
|                   | the evaporating surface to the atmosphere.                                         |                                                             |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`R_{ns}`    | Net incoming shortwave (solar) radiations:                                         | [MJ.m\ :sup:`-2`.d\ :sup:`-1`]                              |
|                   |                                                                                    |                                                             |
|                   | .. math::                                                                          |                                                             |
|                   |   R_{ns} = (1-\alpha) R_{s}                                                        |                                                             |
|                   |                                                                                    |                                                             |
|                   | where :math:`\alpha` is the albedo and :math:`R_s` the incoming solar radiations.  |                                                             |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`R_{nl}`    | Net incoming longwave (solar) radiations.                                          | [MJ.m\ :sup:`-2`.d\ :sup:`-1`]                              |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+
| :math:`G`         | Soil heat flux density.                                                            | [MJ.m\ :sup:`-2`.d\ :sup:`-1`]                              |
+-------------------+------------------------------------------------------------------------------------+-------------------------------------------------------------+


Other symbols that will appear in the following are reported in this table

.. tabularcolumns:: |c|l|c|

+-------------------+-------------------------------------------------------+---------------------------------+
| Symbol            | Description                                           | Units                           |
+===================+=======================================================+=================================+
| :math:`T_{min}`   | Minimum temperature for the time step.                | [:math:`\textrm{\textcelsius}`] |
+-------------------+-------------------------------------------------------+---------------------------------+
| :math:`T_{max}`   | Maximum temperature for the time step.                | [:math:`\textrm{\textcelsius}`] |
+-------------------+-------------------------------------------------------+---------------------------------+
| :math:`T_{avg}`   | Average temperature for the time step                 | [:math:`\textrm{\textcelsius}`] |
|                   | :math:`T_{avg} = (T_{max}+T_{min})/2`.                |                                 |
+-------------------+-------------------------------------------------------+---------------------------------+
| :math:`\Delta T`  | Temperature range for the time step                   | [:math:`\textrm{\textcelsius}`] |
|                   | :math:`\Delta T = (T_{max} - T_{min})`.               |                                 |
+-------------------+-------------------------------------------------------+---------------------------------+
| :math:`D`         | Number of daylight hours.                             |                                 |
+-------------------+-------------------------------------------------------+---------------------------------+
| :math:`R_a`       | Extra-terrestrial solar radiations.                   |                                 |
+-------------------+-------------------------------------------------------+---------------------------------+
| :math:`R_s`       | Incoming solar radiations.                            | [MJ.m\ :sup:`-2`.d\ :sup:`-1`]  |
|                   | By default, :math:`Rs = K_t \sqrt{\Delta T} R_a`.     |                                 |
+-------------------+-------------------------------------------------------+---------------------------------+



References
==========

.. [Allen_etal_1998]
   Allen, R.G., Pereira, L.S., Raes, D. and Smith, M. 1998.
   Crop evapotranspiration - Guidelines for computing crop water requirements
   requirements. *Irrigation and drainage paper 56*.
   United Nations Food and Agriculture Organization, Rome, Italy.
   url=http://www.fao.org/docrep/X0490E/X0490E00.htm ",
.. [Doorenbos_1977]
   Doorenbos, J. and W. O. Pruitt. 1977. Crop water requirements.
   *Irrigation and Drainage Paper 24 2nd ed.*, United Nations Food and Agriculture
   Organization, Rome, Italy.
.. [Droogers_Allen_2002]
   Droogers, P. and Allen, R.G. 2002. Estimating reference evapotranspiration
   under inaccurate data conditions. *Irrigation and Drainage Systems*, 16(1):33-45.
.. [Duffie_Beckman_1981]
   Duffie, J.A., and Beckman, W.A. 1981. *Solar engineering of thermal 
   processes*, 2nd ed. John Wiley & Sons, Inc., New York, N.Y.
.. [Hamon_1961]
   Hamon, W.R. 1961. Estimating the potential evapotranspiration. 
   *ASCE J. Hydraul. Div.*, 87:107-120.
.. [Hansen_1957]
   Hansen, X. 19xx.
.. [Hargreaves_1994]
   Hargreaves, G.H. 1994. Defining and using reference evapotranspiration.
   *Journal of Irrigation and Drainage Engineering*, 120(6):1132-1139. 
.. [Hargreaves_Samani_1985]
   Hargreaves, G.H. and Samani, Z.A. 1985. Reference Crop Evapotranspiration
   from Temperature. *Applied Engineering in Agriculture* 1(2):96-99.
.. [Kharrufa_1985]
   Kharrufa,  N. S. 1985. Simplified equation for evapotranspiration in arid regions.
   *Beitr\:age zur Hydrologie*, Sonderheft 5.1:39-47.
.. [Makkink_1957]
   Makkink, G.F. 1957. Testing the Penman formula by means of lysimeter,
   *J. Inst. of Water Eng.*, 11:277-288.
.. [Penman_Monteith_19xx]
   Penman, xx and Monteith, xx.
.. [Priestley_Taylor_1972]
   Priestley, C. H. B. and R. J. Taylor. 1972. On the assessment of the surface
   heat flux and evaporation using large-scale parameters.
   *Monthly Weather Review*, 100, 81-92.
.. [Spencer_1971]
   Spencer, J.W. 1971. Fourier series representation of the position of the Sun. 
   *Search*, 2(5), 172.
.. [Turc_1961]
   Turc, L. 1961. Estimation of Irrigation Water Requirements, Potential
   Evapotranspiration: A simple climatic formula evovled up to date,
   *Ann. Agronomy*, 12:13-49.
.. [Wright_1982]
    Wright, J.L. 1982. New evapotranspiration crop coefficients.
    *Journal of Irrigation and Drainage Division* 108(1):57-74.
.. [REF_TO_FIND_19xx]
   Reference to find.
