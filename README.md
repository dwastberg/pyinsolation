Python functions for computing the direct and diffuse solar irradition on a surface using the Bird and Hulstrom
insolation algorithm.  Includes helper functions for calculating julian days and sun positions.

Direct port of insolation and related function from the R insol package.  See
http://www.meteoexploration.com/R/insol/index.html for references and more details

## Usage

Calculate Julian day of 15:30 on March 27th, 2016

    julianday=jday_from_datetime(2016,3,27,15,30)
or

    dt=datetime.datetime(2016,3,27,15,30)
    julianday=jday_from_datetime(dt)

Find the sun's position at a given time and position on earth

    #latitude and longitude in decimal degrees
    azimuth,zenith=sunposition(julianday,latitude,longitude,timezone=0)

Calculate direct and diffuse irradiance

    height = 50 #meters above sea level
    visibility = 30 #kilometers
    RH= 40  #relative humidity (%)
    tempK = 278 #Air temperature (K)
    O3 = 0.02 #Ozone thickness (m)
    albedo=0.2 #Albedo of the surrounding terrain [0 to 1]

    direct,diffuse=insolation(zenith, julianday, height, visibility, RH, tempK, O3, albedo):




