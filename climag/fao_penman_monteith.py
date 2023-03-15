"""fao_penman_monteith.py

The FAO Penman-Monteith equation to derive daily evapotranspiration using the
Met Éireann Reanalysis (MÉRA) dataset.

Allen, R. G., Pereira, L. S., Raes, D. and Smith, M. (1998). 'Reference
Evapotranspiration (ETo)', in Crop evapotranspiration: Guidelines for
computing crop water requirements, FAO Irrigation and Drainage Paper, Rome,
Italy, FAO - Food and Agriculture Organization of the United Nations, pp.
15-86. [Online]. Available at https://www.fao.org/3/X0490E/X0490E00.htm
(Accessed 16 August 2022).

Richards, M. (2022). 'PyETo', [Online]. Available at
https://github.com/woodcrafty/PyETo (Accessed 8 November 2022).

PyETo documentation:
https://pyeto.readthedocs.io/en/latest/fao56_penman_monteith.html

MÉRA variables:
- t_max [K]
- t_min [K]
- r_ns [J m⁻²]
- r_nl [J m⁻²]
- p_atm [Pa]
- u_10 [m s⁻¹]
- rh_mean

Overall equation:

```py
import math

t_mean = (t_max + t_min) / 2 - 273.15

e_s = (
    0.6108 * math.exp((17.27 * t_max - 273.15) / (t_max - 273.15 + 237.3)) +
    0.6108 * math.exp((17.27 * t_min - 273.15) / (t_min - 273.15 + 237.3))
) / 2

delta = (4098 * e_s) / math.pow((t_mean + 273.3), 2)

r_n = (r_ns - r_nl) / 1e6

gamma = 0.665 * p_atm / 1e6

u_2 = u_10 * (4.87 / math.log((67.8 * 10) - 5.42))

e_deficit = e_s - rh_mean / 100 * e_s

ETo = (
    (
        (0.408 * delta * r_n) +
        (gamma * (900 / (t_mean + 273)) * u_2 * e_deficit)
    ) /
    (delta + (gamma * (1 + 0.34 * u_2)))
)
```
"""

import math


def mean_air_temperature(t_min: float, t_max: float) -> float:
    """
    Due to the non-linearity of humidity data required in the FAO
    Penman-Monteith equation, the vapour pressure for a certain period should
    be computed as the mean between the vapour pressure at the daily maximum
    and minimum air temperatures of that period.

    The daily maximum air temperature and daily minimum air temperature are,
    respectively, the maximum and minimum air temperature observed during the
    24-hour period, beginning at midnight.

    The mean daily air temperature is only employed in the FAO Penman-Monteith
    equation to calculate the slope of the saturation vapour pressure curves
    and the impact of mean air density as the effect of temperature variations
    on the value of the climatic parameter is small in these cases.

    For standardisation, the mean temperature for 24-hour periods is defined
    as the mean of the daily maximum and minimum temperatures rather than as
    the average of hourly temperature measurements.

    Equation (9) in Allen et al. (1998).

    T_mean = (T_max + T_min) / 2

    - T_mean: mean air temperature at 2 m height [°C]
    - T_max: maximum air temperature at 2 m height [°C]
    - T_min: minimum air temperature at 2 m height [°C]

    Parameters
    ----------
    t_min : Minimum air temperature at 2 m height [°C]
    t_max : Maximum air temperature at 2 m height [°C]

    Returns
    -------
    - Mean air temperature at 2 m height [°C]
    """

    return (t_max + t_min) / 2


def saturation_vapour_pressure_temp(t_air: float) -> float:
    """
    Equation (11) in Allen et al. (1998).

    e^o(T) = 0.6108exp(17.27T / (T + 237.3))

    - e^o(T): saturation vapour pressure at the air temperature T [kPa]
    - T: air temperature [°C]

    Parameters
    ----------
    t_air : Air temperature at 2 m height [°C]

    Returns
    -------
    - Saturation vapour pressure at the given air temperature [kPa]
    """

    return 0.6108 * math.exp((17.27 * t_air) / (t_air + 237.3))


def saturation_vapour_pressure(t_max: float, t_min: float) -> float:
    """
    Equation (12) in Allen et al. (1998).

    e_s = (e^o(T_max) + e^o(T_min)) / 2

    - e_s: Mean saturation vapour pressure [kPa]
    - e^o(T_max): saturation vapour pressure at the maximum air temperature
      [kPa]
    - e^o(T_min): saturation vapour pressure at the minimum air temperature
      [kPa]

    Parameters
    ----------
    t_min : Minimum air temperature at 2 m height [°C]
    t_max : Maximum air temperature at 2 m height [°C]

    Returns
    -------
    - Mean saturation vapour pressure [kPa]
    """

    return (
        saturation_vapour_pressure_temp(t_air=t_max) +
        saturation_vapour_pressure_temp(t_air=t_min)
    ) / 2


def slope_vapour_pressure_curve(t_mean: float, e_s: float) -> float:
    """
    The slope of the relationship between saturation vapour pressure and
    temperature.

    Equation (13) in Allen et al. (1998).

    Δ = (4098(0.6108exp(17.27T / (T + 237.3)))) / (T + 237.3)²
      = (4098e_s) / (T + 237.3)²

    - Δ: slope vapour pressure curve [kPa °C⁻¹]
    - T: mean air temperature at 2 m height [°C]
    - e_s: mean saturation vapour pressure [kPa]

    For standardisation, the mean temperature for 24-hour periods is defined
    as the mean of the daily maximum and minimum temperatures rather than as
    the average of hourly temperature measurements.

    Parameters
    ----------
    mean_t : Mean air temperature at 2 m height [°C]
    e_s : Mean saturation vapour pressure [kPa]

    Returns
    -------
    - Slope vapour pressure curve [kPa °C⁻¹]
    """

    return (4098 * e_s) / math.pow((t_mean + 273.3), 2)


def net_radiation(r_ns: float, r_nl: float) -> float:
    """
    The net radiation is the difference between the incoming net shortwave
    radiation and the outgoing net longwave radiation.

    Equation (40) in Allen et al. (1998)

    R_n = R_ns - R_nl

    - R_n: net radiation at the crop surface [MJ m⁻² day⁻¹]
    - R_ns: incoming net shortwave radiation [MJ m⁻² day⁻¹]
    - R_nl: outgoing net longwave radiation [MJ m⁻² day⁻¹]

    Parameters
    ----------
    r_ns : Incoming net shortwave radiation [MJ m⁻² day⁻¹]
    r_nl : Outgoing net longwave radiation [MJ m⁻² day⁻¹]

    Returns
    -------
    - Net radiation at the crop surface [MJ m⁻² day⁻¹]
    """

    return r_ns - r_nl


# def soil_heat_flux_density(shf: float) -> float:
#     """
#     The soil heat flux is small compared to net radiation, particularly when
#     the surface is covered by vegetation. As the magnitude of the day heat
#     flux beneath the grass reference surface is relatively small, it may be
#     ignored. This is shown in Equation (42) in Allen et al. (1998).

#     Based on Nolan and Flanagan (2020)'s derivation of evapotranspiration
#     using the FAO Penman-Monteith equation, the soil heat flux density is
#     equivalent to the mean surface sensible heat flux.

#     G = 0.0864SH

#     - G: soil heat flux density [MJ m⁻² day⁻¹]
#     - SH: mean surface sensible heat flux [W m⁻²]

#     Parameters
#     ----------
#     shf : Mean surface sensible heat flux [J m⁻² day⁻¹]

#     Returns
#     -------
#     - Soil heat flux density [MJ m⁻² day⁻¹]
#     """

#     return shf / 1e6


# def atmospheric_pressure(elevation: float) -> float:
#     """
#     The atmospheric pressure based on Equation (7) in Allen et al. (1998).

#     P = 101.3((293 - 0.0065z) / 293)^5.26

#     - P: atmospheric pressure [kPa]
#     - z: elevation above sea level [m]

#     Parameters
#     ----------
#     elevation : Elevation above sea level [m]

#     Returns
#     -------
#     - Atmospheric pressure [kPa]
#     """

#     return 101.3 * math.pow(((293 - 0.0065 * elevation) / 293), 5.26)


def psychrometric_constant(p_atm: float) -> float:
    """
    The psychrometric constant is the ratio of specific heat of moist air at
    constant pressure to latent heat of vaporisation of water. It can be
    estimated from atmospheric pressure.

    This method assumes that the air is saturated with water vapour at the
    minimum daily temperature.

    Equation (8) in Allen et al (1998).

    γ = 0.000665P

    - γ: psychrometric constant [kPa °C⁻¹]
    - P: atmospheric pressure [kPa]

    Parameters
    ----------
    p_atm : Atmospheric pressure [kPa]

    Returns
    -------
    - Psychrometric constant [kPa °C⁻¹]
    """

    return 0.665 / 1e3 * p_atm


def wind_speed(u_z: float, h_z: float) -> float:
    """
    If not measured at 2 m height, convert wind speed measured at different
    heights above the soil surface to wind speed at 2 m above the surface,
    assuming a short grass surface.

    Equation (47) in Allen et al (1998).

    u₂ = u_z * (4.87 / ln(67.8z - 5.42))

    - u₂: wind speed at 2 m height [m s⁻¹]
    - u_z: measured wind speed at z m above ground surface [m s⁻¹]
    - z: height of measurement above ground surface [m]

    Parameters
    ----------
    u_z : Measured wind speed at z m above ground surface [m s⁻¹]
    h_z : Height of measurement above ground surface [m]

    Returns
    -------
    - Wind speed at 2 m height [m s⁻¹]
    """

    return u_z * (4.87 / math.log((67.8 * h_z) - 5.42))


def actual_vapour_pressure(rh_mean: float, e_s: float) -> float:
    """
    In the absence of dewpoint temperature, or psychrometric data (i.e. wet
    and dry bulb temperatures), or maximum (and minimum) relative humidity,
    the mean relative humidity is used to calculate the actual vapour
    pressure.

    Equation (19) in Allen et al. (1998).

    e_a = (RH_mean / 100)((e^o(T_max) + e^o(T_min)) / 2)
        = (RH_mean / 100)e_s

    - e_a: actual vapour pressure [kPa]
    - RH_mean: mean relative humidity
    - e_s: saturation vapour pressure [kPa]

    Parameters
    ----------
    rh_mean : Relative humidity at 2 m height
    e_s : Saturation vapour pressure [kPa]

    Returns
    -------
    - Actual vapour pressure [kPa]
    """

    return rh_mean / 100 * e_s


def saturation_vapour_pressure_deficit(e_s: float, e_a: float) -> float:
    """
    Equivalent to e_s - e_a

    - e_s: saturation vapour pressure [kPa]
    - e_a: actual vapour pressure [kPa]

    Parameters
    ----------
    e_s : Saturation vapour pressure [kPa]
    e_a : Actual vapour pressure [kPa]

    Returns
    -------
    - Saturation vapour pressure deficit [kPa]
    """

    return e_s - e_a


def fao_penman_monteith(
    r_n: float, t_mean: float, u_2: float, e_deficit: float,
    delta: float, gamma: float
) -> float:
    """
    The FAO Penman-Monteith equation.

    Equation (6) in Allen et al. (1998)

    ETo = (
        (0.408Δ(R_n - G) + γ(900 / (T + 273))u₂(e_s - e_a)) /
        (Δ + γ(1 + 0.34u₂))
    )

    ETo: reference evapotranspiration [mm day⁻¹]
    Δ: slope vapour pressure curve [kPa °C⁻¹]
    R_n: net radiation at the crop surface [MJ m⁻² day⁻¹]
    G: soil heat flux density [MJ m⁻² day⁻¹]
    γ: psychrometric constant [kPa °C⁻¹]
    T: mean daily air temperature at 2 m height [°C]
    u₂: wind speed at 2 m height [m s⁻¹]
    e_s: saturation vapour pressure [kPa]
    e_a: actual vapour pressure [kPa]

    The soil heat flux is small compared to net radiation, particularly when
    the surface is covered by vegetation. As the magnitude of the day heat
    flux beneath the grass reference surface is relatively small, it may be
    ignored. This is shown in Equation (42) in Allen et al. (1998).

    Parameters
    ----------
    gamma : Slope vapour pressure curve [kPa °C⁻¹]
    r_n : Net radiation at the crop surface [MJ m⁻² day⁻¹]
    gamma : Psychrometric constant [kPa °C⁻¹]
    t_mean : Mean daily air temperature at 2 m height [°C]
    u_2 : Wind speed at 2 m height [m s⁻¹]
    e_deficit : Saturation vapour pressure deficit [kPa]

    Returns
    -------
    - Reference evapotranspiration [mm day⁻¹]
    """

    return (
        (
            (0.408 * delta * r_n) +
            (gamma * (900 / (t_mean + 273)) * u_2 * e_deficit)
        ) /
        (delta + (gamma * (1 + 0.34 * u_2)))
    )


def fao_penman_monteith_mera(data):
    """
    Apply the FAO Penman-Monteith equation on the MÉRA dataset.

    Parameters
    ----------
    data : Xarray dataset
    """

    data = data.assign(
        t_mean=mean_air_temperature(
            t_min=data["t_min"] - 273.15, t_max=data["t_max"] - 273.15
        )
    )

    data = data.assign(
        e_s=saturation_vapour_pressure(
            t_max=data["t_max"] - 273.15, t_min=data["t_min"] - 273.15
        )
    )

    data = data.assign(
        delta=slope_vapour_pressure_curve(
            t_mean=data["t_mean"], e_s=data["e_s"]
        )
    )

    data = data.assign(
        r_n=net_radiation(r_ns=data["r_ns"] / 1e6, r_nl=data["r_nl"] / 1e6)
    )

    data = data.assign(
        gamma=psychrometric_constant(p_atm=data["p_atm"] / 1e3)
    )

    data = data.assign(u_2=wind_speed(u_z=data["u_10"], h_z=10))

    data = data.assign(
        e_a=actual_vapour_pressure(rh_mean=data["rh_mean"], e_s=data["e_s"])
    )

    data = data.assign(
        e_deficit=saturation_vapour_pressure_deficit(
            e_s=data["e_s"], e_a=data["e_a"]
        )
    )

    data = data.assign(
        PET=fao_penman_monteith(
            delta=data["delta"], gamma=data["gamma"], r_n=data["r_n"],
            t_mean=data["t_mean"], u_2=data["u_2"],
            e_deficit=data["e_deficit"]
        )
    )

    data = data.drop_vars(
        [
            "t_min", "t_max", "t_mean", "e_s", "delta", "r_ns", "r_nl", "r_n",
            "gamma", "p_atm", "u_10", "u_2", "e_a", "rh_mean", "e_deficit"
        ]
    )

    return data
