module constants
    !! Predefined constants used in DGSWEM.
    !!
    !! Introduced August 27, 2026 by Ashton Cole.

    use sizes, only : sz
    
    implicit none

    real(sz), parameter :: prbckgrnd = 1013.0
    !! Background atmospheric pressure, in mb. Used in wind calculations.
    real(sz), parameter :: rhoair = 1.15d0
    !! Density of air, in kg/m^3. Used in wind and vortex calculations.
    real(sz), parameter :: windreduction = 0.9d0
    !! Scaling factor for wind velocity. Used in wind and vortex calculations.
    real(sz), parameter :: one2ten = 0.8928d0
    !! Conversion factor from one-minute averaged wind velocity to ten-minute
    !! averaged wind velocity. This accounts for short-duration gusts. Used in
    !! vortex calculations.
    real(sz), parameter :: ten2one = 1.0/0.8928d0
    !! Conversion factor from ten-minute averaged wind velocity to one-minute
    !! averaged wind velocity. This accounts for short-duration gusts. Used in
    !! vortex calculations.
    real(sz), parameter :: wavewindmultiplier = 1.0
    !! UNUSED?
    !! Scaling factor on winds sent to SWAN.
    real(sz), parameter ::  pi = 3.141592653589793d0
    !! Mathematical constant pi.
    real(sz), parameter ::  deg2rad = pi / 180.0
    !! Conversion factor from degrees to radians.
    real(sz), parameter ::  rad2deg = 180.0 / pi
    !! Conversion factor from radians to degrees.
    real(sz), parameter ::  e = 2.718281828459045d0
    !! Euler's constant.
    real(sz), parameter :: day2sec = 24.0*3600.0
    !! UNUSED?
    !! Conversion factor from days to seconds.
    real(sz), parameter :: sec2day = 1.0 / day2sec
    !! UNUSED?
    !! Conversion factor from seconds to days.
    real(sz), parameter :: rearth = 6378206.4d0
    !! Radius of Earth, in m.
    real(sz), parameter :: nm2m = 1852.0
    !! Conversion factor from nautical miles to meters.
    real(sz), parameter :: m2nm = 1.0 / nm2m
    !! Conversion factor from meters to nautical miles.
    real(sz), parameter :: kt2ms = nm2m / 3600.0
    !! Conversion factor from knots to m/s.
    real(sz), parameter :: ms2kt = 1.0 / kt2ms
    !! Conversion factor from m/s to knots.
    real(sz), parameter :: omega = 2.0 * pi / 86164.2d0
    !! Earth's angular velocity, in radians/s. Used to compute Coriolis force.
    real(sz), parameter :: mb2pa = 100.0
    !! Conversion factor from milibars to Pascals.
    real(sz), parameter :: barmin = 0.01d0
    !! Minimum height threshold or tolerance for barrier calculations. This
    !! filters out negligible height differences.
    real(sz), parameter :: depavg = 0.0d0
    !! UNUSED?
    real(sz), parameter :: depmax = 0.0d0
    !! UNUSED?
    real(sz), parameter :: depmin = 0.0d0
    !! UNUSED?
    
end module constants