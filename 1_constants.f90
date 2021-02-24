! CONSTANTS -----------------------------------------------------------

module constants
    integer, parameter :: dp = kind(0.0d0)
    real(dp), parameter :: pi = 4.0*atan(1.0)
    real(dp), parameter :: ne = exp(1.0)
    complex(dp), parameter :: Im = (0.0,1.0)

end module constants

