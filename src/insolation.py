from __future__ import division

from math import radians, degrees, cos, sin, tan, asin, acos, atan2, exp, floor, pi


def jday_from_datetime(year_or_datetime, month=None, day=None, hour=12, minute=0, sec=0):
    """computes the Julian Day from from a date and time.  Valid for dates between 1901 and 2099

    :param year_or_datetime: Either the year as an int or a python datetime object. If first argument is a python
    datetime object, the rest of the arguments are ignored otherwise.
    :param month: month of year as int in the range 1-12
    :param day: day of the month
    :param hour: hour of day
    :param minute: minute
    :param second: second
    :return:  Julian day of the supplied date

    >>> print round(jday_from_datetime(2000,1,1,12,0,0),4)
    2451545.0
    >>> print round(jday_from_datetime(2016,3,9,16,10,30),4)
    2457457.1944
    >>> print round(jday_from_datetime(datetime.datetime(2016,3,9,16,10,30)),4)
    2457457.1944
    """
    try:
        year = year_or_datetime.year
        month = year_or_datetime.month
        day = year_or_datetime.day
        hour = year_or_datetime.hour
        minute = year_or_datetime.minute
        sec = year_or_datetime.second
    except AttributeError:
        year = year_or_datetime
        if month is None or day is None:
            raise ValueError("month and day must be set")

    hour = hour + minute / 60 + sec / 60
    jd = 367 * year - (7 * (year + (month + 9) // 12)) // 4 + (275 * month) // 9 + day + 1721013.5 + hour / 24
    return jd


def sunposition(jd, lat, lon, timezone=0):
    """calculate the suns position in the sky at a given time and place. No refraction, center of disc

    :param jd: time in Julian days
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :param timezone: Time zone, west is negative.
    :return: azimuth and zenith angles of the sun at the given time and place

    >>> print map(lambda x:round(x,4),sunposition(2457457.16667,54,12,1))
    [237.7328, 73.5839]
    >>> print map(lambda x:round(x,4),sunposition(2452460,0,0,0))
    [2.5879, 22.8911]
    """
    sunv = sunvector(jd, lat, lon, timezone)
    azimuth = degrees(pi - atan2(sunv[0], sunv[1]))
    zenith = degrees(acos(sunv[2]))

    return (azimuth, zenith)


def sunvector(jd, lat, lon, timezone=0):
    """Calculates a unit vector in the direction of the sun from the observer position.
    :param jd: time in Julian days
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :param timezone: Time zone, west is negative.
    :return: unit vector in the direction of the sun from the observer position

    >>> print map(lambda x:round(x,4),sunvector(2457457.16667,54,12,1))
    [-0.8111, 0.5121, 0.2826]
    >>> print map(lambda x:round(x,4),sunvector(2457000,12,54,-8))
    [-0.0639, 0.1867, -0.9803]
    """
    omegar = _hourangle(jd, lon, timezone)
    deltar = radians(declination(jd))
    lambdar = radians(lat)
    svx = -sin(omegar) * cos(deltar)
    svy = sin(lambdar) * cos(omegar) * cos(deltar) - cos(lambdar) * sin(deltar)
    svz = cos(lambdar) * cos(omegar) * cos(deltar) + sin(lambdar) * sin(deltar)
    return (svx, svy, svz)


def _hourangle(jd, longitude, timezone):
    """
    internal function for solar position

    >>> print round(_hourangle(2457457.16667,54,12),4)
    -1.1971
    >>> print round(_hourangle(2457457.16667,12,54),4)
    -12.9257

    """
    hour = ((jd - floor(jd)) * 24 + 12) % 24
    eqt = eqtime(jd)
    stndmeridian = timezone * 15
    deltalontime = longitude - stndmeridian
    deltalontime = deltalontime * 24.0 / 360.0
    omegar = pi * (((hour + deltalontime + eqt / 60) / 12.0) - 1.0)
    return omegar


def eqtime(jd):
    """Computes the equation of time for a given Julian Day

    >>> print round(eqtime(2457561.06944),4)
    -1.9076
    >>> print round(eqtime(2457457.16667),4)
    -10.3537
    """
    jdc = (jd - 2451545.0) / 36525.0
    sec = 21.448 - jdc * (46.8150 + jdc * (0.00059 - jdc * (0.001813)))
    e0 = 23.0 + (26.0 + (sec / 60.0)) / 60.0
    ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)
    oblcorr = e0 + 0.00256 * cos(radians(125.04 - 1934.136 * jdc))
    y = (tan(radians(oblcorr) / 2)) ** 2
    l0 = 280.46646 + jdc * (36000.76983 + jdc * (0.0003032))
    l0 = (l0 - 360 * (l0 // 360)) % 360
    rl0 = radians(l0)
    gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
    gmas = radians(gmas)
    EqTime = y * sin(2 * rl0) - 2.0 * ecc * sin(gmas) + 4.0 * ecc * y * sin(gmas) * cos(2 * rl0) - \
             0.5 * y ** 2 * sin(4 * rl0) - 1.25 * ecc ** 2 * sin(2 * gmas)
    return (degrees(EqTime) * 4)


def declination(jd):
    """Computes the declination of the Sun for a given Julian Day

    >>> print round(declination(2457457.16667),4)
    -4.1502
    >>> print round(declination(2457561.06944),4)
    23.4333
    >>>
    """
    jdc = (jd - 2451545.0) / 36525.0
    sec = 21.448 - jdc * (46.8150 + jdc * (0.00059 - jdc * (0.001813)))
    e0 = 23.0 + (26.0 + (sec / 60.0)) / 60.0
    oblcorr = e0 + 0.00256 * cos(radians(125.04 - 1934.136 * jdc))
    l0 = 280.46646 + jdc * (36000.76983 + jdc * (0.0003032))
    l0 = (l0 - 360 * (l0 // 360)) % 360
    gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
    gmas = radians(gmas)
    seqcent = sin(gmas) * (1.914602 - jdc * (0.004817 + 0.000014 * jdc)) + \
              sin(2 * gmas) * (0.019993 - 0.000101 * jdc) + sin(3 * gmas) * 0.000289
    suntl = l0 + seqcent
    sal = suntl - 0.00569 - 0.00478 * sin(radians(125.04 - 1934.136 * jdc))
    delta = asin(sin(radians(oblcorr)) * sin(radians(sal)))
    return (degrees(delta))


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
    :return: [Direct irradiance, Diffuse irradiance]

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
    import datetime

    doctest.testmod()
