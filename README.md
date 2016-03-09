Python functions for computing the direct and diffuse solar irradition on a surface using the Bird and Hulstrom
insolation algorithm.  Includes helper functions for calculating julian days and sun positions.

## Usage

Calculate Julian day of 15:30 on March 27th, 2016

    jday_from_datetime(2016,3,27,15,30)
or

    dt=datetime.datetime(2016,3,27,15,30)
    jday_from_datetime(dt)

Find the sun's position at a given time and position on earth

    azimuth,zenith=sunposition(julianday,latitude,longitude,timezone=0)



Direct port of insolation and related function from the R insol package
(http://www.meteoexploration.com/R/insol/index.html)
