module debye
contains
double precision function cheval ( n, a, t )

!*****************************************************************************80
!
!! CHEVAL evaluates a Chebyshev series.
!
!  Discussion:
!
!    This function evaluates a Chebyshev series, using the
!    Clenshaw method with Reinsch modification, as analysed
!    in the paper by Oliver.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon 
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    J Oliver,
!    An error analysis of the modified Clenshaw method for
!    evaluating Chebyshev and Fourier series,
!    Journal of the IMA, 
!    Volume 20, 1977, pages 379-391.
!
!  Parameters:
!
!    Input, integer N, the number of terms in the sequence.
!
!    Input, real ( kind = 8 ) A(0:N), the coefficients of the Chebyshev series.
!
!    Input, real ( kind = 8 ) T, the value at which the series is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CHEVAL, the value of the Chebyshev series at T.
!
  implicit none

  integer n

  real ( kind = 8 ) :: a(0:n)
  real ( kind = 8 ) :: d1
  real ( kind = 8 ) :: d2
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  real ( kind = 8 ) :: t
  real ( kind = 8 ), parameter :: test = 0.6D+00
  real ( kind = 8 ) :: tt
  real ( kind = 8 ), parameter :: two = 2.0D+00
  real ( kind = 8 ) :: u0
  real ( kind = 8 ) :: u1
  real ( kind = 8 ) :: u2
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  u1 = zero
!
!  T <= -0.6, Reinsch modification.
!
  if ( t <= -test ) then

    d1 = zero
    tt = ( t + half ) + half
    tt = tt + tt

    do i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) - d2
      u1 = d1 - u2
    end do

    cheval = ( d1 - d2 ) / two
!
!  -0.6 < T < 0.6, Standard Clenshaw method.
!
  else if ( t < test ) then

    u0 = zero
    tt = t + t

    do i = n, 0, -1
      u2 = u1
      u1 = u0
      u0 = tt * u1 + a(i) - u2
    end do

    cheval = ( u0 - u2 ) / two
!
!  0.6 <= T, Reinsch modification.
!
  else

    d1 = zero
    tt = ( t - half ) - half
    tt = tt + tt

    do i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) + d2
      u1 = d1 + u2
    end do

    cheval = ( d1 + d2 ) / two

  end if

  return
end function cheval

double precision function debye3 ( xvalue )

!*****************************************************************************80
!
!! DEBYE3 calculates the Debye function of order 3.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE3(x) = 3 / x^3 * Integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon 
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE3, the value of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: eight = 8.0D+00
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer i
  integer nexp
  integer, parameter :: nterms = 16
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xvalue
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  real ( kind = 8 ) adeb3(0:18),debinf,expmx, &
       pt375,rk,sevp5,sum1,t,twenty, &
       xk,xki,xlim1,xlim2,xlow,xupper
  data pt375/0.375d0/
  data sevp5,twenty/7.5d0 , 20.0d0/
  data debinf/0.51329911273421675946d-1/
  data adeb3/2.70773706832744094526d0, &
             0.34006813521109175100d0, &
            -0.1294515018444086863d-1, &
             0.79637553801738164d-3, &
            -0.5463600095908238d-4, &
             0.392430195988049d-5, &
            -0.28940328235386d-6, &
             0.2173176139625d-7, &
            -0.165420999498d-8, &
             0.12727961892d-9, &
            -0.987963459d-11, &
             0.77250740d-12, &
            -0.6077972d-13, &
             0.480759d-14, &
            -0.38204d-15, &
             0.3048d-16, &
            -0.244d-17, &
             0.20d-18, &
            -0.2d-19/
!
!   Machine-dependent constants
!
  data xlow,xupper/0.298023d-7,35.35051d0/
  data xlim1,xlim2/708.39642d0,0.9487163d103/
!
  x = xvalue

  if ( x < zero ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DEBYE3 - Fatal error!'
    write ( *, '(a)' ) '  Argument X < 0.'
    debye3 = zero
    return
  end if

  if ( x < xlow ) then
    debye3 = ( ( x - sevp5 ) * x + twenty ) / twenty
  else if ( x <= 4 ) then
    t = ( ( x * x / eight ) - half ) - half
    debye3 = cheval ( nterms, adeb3, t ) - pt375 * x
  else
!
!   Code for x > 4.0
!
     if ( xlim2 < x ) then
        debye3 = zero
     else
        debye3 = one / ( debinf * x * x * x )
        if ( x < xlim1 ) then
           expmx = exp ( -x )
           if ( xupper < x ) then
              sum1 = ((( x + three ) * x + six ) * x + six ) / ( x * x * x )
           else
              sum1 = zero
              rk = aint ( xlim1 / x )
              nexp = int ( rk )
              xk = rk * x
              do i = nexp, 1, -1
                 xki = one / xk
                 t =  ((( six * xki + six ) * xki + three ) * xki + one ) / rk
                 sum1 = sum1 * expmx + t
                 rk = rk - one
                 xk = xk - x
              end do
           end if
           debye3 = debye3 - three * sum1 * expmx
        end if
     end if
  end if

  return
end function debye3
end module debye
