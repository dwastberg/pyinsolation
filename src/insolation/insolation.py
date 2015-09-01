import numpy as np
from math import radians, cos, sin, exp


def insolation(zenith, jd, height, visibility, RH, tempK, O3, alphag):
    """Computes direct and diffuse solar irradiance perpendicular to the beam,
    for a given zenith angle, Julian Day, altitude and atmospheric conditions.

    :param zenith:  Zenith angle in degrees.
    :param jd: Julian Day.
    :param height: Altitude above sea level.
    :param visibility:  Visibility [km].
    :param RH: Relative humidity [%].
    :param tempK: Air temperature [K].
    :param O3: Ozone thickness [m].
    :param alphag: Albedo of the surrounding terrain [0 to 1].

    >>> print map(lambda f: round(f,4),insolation(30,0,3200,28,60,278.15,0.02,0.2))
    [952.9052, 134.0539]
    >>> print map(lambda f: round(f,4),insolation(0,0,0,10,20,278.15,0.1,0.5))
    [771.9925, 333.1886]

    """
    if tempK <= 153:
        raise ValueError("temperature should be in Kelvin")

    Isc = 1361.0
    theta = radians(zenith)
    ssctalb = 0.9  # single scattering albedo (aerosols)(Iqbal, 1983)
    Fc = 0.84  # ratio of forward to total energy scattered (Iqbal, 1983)
    Pz = z2p(height)
    Mr = 1.0 / (cos(theta) + 0.15 * ((93.885 - zenith) ** (-1.253)))
    Ma = Mr * Pz / 1013.25
    wvap_s = wvapsat(tempK)
    Wprec = 46.5 * (RH / 100.0) * wvap_s / tempK  # Prata 1996
    rho2 = sunr(jd)
    TauR = exp((-.09030 * (Ma ** 0.84)) * (1.0 + Ma - (Ma ** 1.01)))
    TauO = 1.0 - ((0.1611 * (O3 * Mr) * (1.0 + 139.48 * (O3 * Mr)) ** (-0.3035)) - 0.002715 * (O3 * Mr) * (
    1.0 + 0.044 * (O3 * Mr) + 0.0003 * (O3 * Mr) ** 2) ** (-1))
    TauG = exp(-0.0127 * (Ma ** 0.26))
    TauW = 1.0 - 2.4959 * (Wprec * Mr) * ((1.0 + 79.034 * (Wprec * Mr)) ** 0.6828 + 6.385 * (Wprec * Mr)) ** (-1)
    TauA = (0.97 - 1.265 * (visibility ** (-0.66))) ** (Ma ** 0.9)  # Machler, 1983
    TauTotal = TauR * TauO * TauG * TauW * TauA
    In = 0.9751 * (1 / rho2) * Isc * TauTotal
    tauaa = 1.0 - (1.0 - ssctalb) * (1.0 - Ma + Ma ** 1.06) * (1.0 - TauA)
    Idr = 0.79 * (1 / rho2) * Isc * cos(theta) * TauO * TauG * TauW * tauaa * 0.5 * (1.0 - TauR) / (
    1.0 - Ma + Ma ** (1.02))
    tauas = (TauA) / tauaa
    Ida = 0.79 * (1 / rho2) * Isc * cos(theta) * TauO * TauG * TauW * tauaa * Fc * (1.0 - tauas) / (
    1.0 - Ma + Ma ** 1.02)
    alpha_atmos = 0.0685 + (1.0 - Fc) * (1.0 - tauas)
    Idm = (In * cos(theta) + Idr + Ida) * alphag * alpha_atmos / (1.0 - alphag * alpha_atmos)
    Id = Idr + Ida + Idm
    return [In, Id]


def z2p(z, P0=101325, T0=288.15):
    """
    Computes air pressure for a given altitude according to the standard atmosphere.

    :param z: altitude above sea level in metres [0:10000]
    :param P0: Pressure at sea level.
    :param T0: Temperature at sea level
    :return: Pressure in hPa

    >>> z2p(0)
    1013.25
    >>> print round(z2p(1000),4)
    898.7592
    >>> print round(z2p(20,100000,250),4)
    997.2699

    """

    earth_G = 9.80665  # Acceleration due to gravity (m s-2)
    earth_R = 6.3756766E6  # Average earths radius (m)
    Md = 28.966  # Molecular weight of dry air
    r_star = 8.3145  # Universal gas constant J/molK
    stlapse = -0.0065  # standard lapse rate K/m
    H1 = (earth_R * z) / (earth_R + z)
    HB = 0.0
    zp = P0 * (T0 / (T0 + stlapse * (H1 - HB))) ** ((earth_G * Md) / (r_star * stlapse * 1000))
    zp = zp / 100.0
    return zp


def wvapsat(tempk, ice=False):
    """
    Computes the saturation pressure of water vapour in air over water or ice

    :param tempk: Air temperature [K]
    :param ice: Over ice or water [True,False]
    :return: Partial pressure of water vapour [hPa].

    >>> print round(wvapsat(273),4)
    6.0371
    >>> print round(wvapsat(300),4)
    35.3136
    """
    if ice:
        tempcl = tempk - 273.15
        a0 = 6.109177956
        a1 = 5.03469897e-1
        a2 = 1.886013408e-2
        a3 = 4.176223716e-4
        a4 = 5.824720280e-6
        a5 = 4.838803174e-8
        a6 = 1.838826904e-10
    else:
        tempcl = tempk
        a0 = 6984.505294
        a1 = -188.9039310
        a2 = 2.133357675
        a3 = -1.288580973e-2
        a4 = 4.393587233e-5
        a5 = -8.023923082e-8
        a6 = 6.136820929e-11
    return a0 + tempcl * (a1 + tempcl * (a2 + tempcl * (a3 + tempcl * (a4 + tempcl * (a5 + tempcl * a6)))))


def sunr(jd):
    """
    Calculates the Earth radius vector

    :param jd: Julian Day
    :return: Earth Radius Vector in Astronomical Units (AU)

    >>> print round(sunr(0),4)
    0.9949
    >>> print round(sunr(10000),4)
    1.0166
    """
    jdc = (jd - 2451545.0) / 36525.0
    ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)
    gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
    gmasr = radians(gmas)
    seqc = sin(gmasr) * (1.914602 - jdc * (0.004817 + 0.000014 * jdc)) + sin(2 * gmas) * (
    0.019993 - 0.000101 * jdc) + sin(3 * gmasr) * 0.000289
    sta = gmas + seqc
    sunrv = (1.000001018 * (1 - ecc ** 2)) / (1 + ecc * cos(radians(sta)))
    return sunrv


if __name__ == "__main__":
    import doctest

    doctest.testmod()
