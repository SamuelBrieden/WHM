module HyperGeo
implicit none
public cchg
contains

subroutine gamma ( x, ga )

!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program that the copyright
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = rk ) GA, the value of the Gamma function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), dimension ( 26 ) :: g = (/ &
    1.0D+00, &
    0.5772156649015329D+00, &
   -0.6558780715202538D+00, &
   -0.420026350340952D-01, &
    0.1665386113822915D+00, &
   -0.421977345555443D-01, &
   -0.96219715278770D-02, &
    0.72189432466630D-02, &
   -0.11651675918591D-02, &
   -0.2152416741149D-03, &
    0.1280502823882D-03, &
   -0.201348547807D-04, &
   -0.12504934821D-05, &
    0.11330272320D-05, &
   -0.2056338417D-06, &
    0.61160950D-08, &
    0.50020075D-08, &
   -0.11812746D-08, &
    0.1043427D-09, &
    0.77823D-11, &
   -0.36968D-11, &
    0.51D-12, &
   -0.206D-13, &
   -0.54D-14, &
    0.14D-14, &
    0.1D-15 /)
  real ( kind = rk ) ga
  real ( kind = rk ) gr
  integer k
  integer m
  integer m1
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r
  real ( kind = rk ) x
  real ( kind = rk ) z

  if ( x == aint ( x ) ) then

    if ( 0.0D+00 < x ) then
      ga = 1.0D+00
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0D+00 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0D+00
      do k = 1, m
        r = r * ( z - real ( k, kind = rk ) )
      end do
      z = z - real ( m, kind = rk )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0D+00 / ( gr * z )

    if ( 1.0D+00 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0D+00 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end

subroutine cchg ( a, b, z, chg )

!*****************************************************************************80
!
!! cchg() computes the confluent hypergeometric function.
!
!  Discussion:
!
!    This function computes the confluent hypergeometric function
!    M(a,b,z) with real parameters a, b and complex argument z.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Input:
!
!    real ( kind = rk ) A, B, parameter values.
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) CHG, the value of M(a,b,z).
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a0
  real ( kind = rk ) a1
  real ( kind = rk ) b
  real ( kind = rk ) ba
  complex ( kind = ck ) cfac
  complex ( kind = ck ) chg
  complex ( kind = ck ) chg1
  complex ( kind = ck ) chg2
  complex ( kind = ck ) chw
  complex ( kind = ck ) ci
  complex ( kind = ck ) cr
  complex ( kind = ck ) cr1
  complex ( kind = ck ) cr2
  complex ( kind = ck ) crg
  complex ( kind = ck ) cs1
  complex ( kind = ck ) cs2
  complex ( kind = ck ) cy0
  complex ( kind = ck ) cy1
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) g3
  integer i
  integer j
  integer k
  integer la
  integer m
  integer n
  integer nl
  integer ns
  real ( kind = rk ) phi
  real ( kind = rk ) pi
  real ( kind = rk ) x
  real ( kind = rk ) x0
  real ( kind = rk ) y
  complex ( kind = ck ) z
  complex ( kind = ck ) z0

  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = ck )
  a0 = a
  a1 = a
  z0 = z

  if ( b == 0.0D+00 .or. b == - int ( abs ( b ) ) ) then
    chg = cmplx ( 1.0D+30, 0.0D+00, kind = ck )
  else if ( a == 0.0D+00 .or. z == 0.0D+00 ) then
    chg = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
  else if ( a == -1.0D+00 ) then
    chg = 1.0D+00 - z / b
  else if ( a == b ) then
    chg = exp ( z )
  else if ( a - b == 1.0D+00 ) then
    chg = ( 1.0D+00 + z / b ) * exp ( z )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
    chg = ( exp ( z ) - 1.0D+00 ) / z
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
    m = int ( - a )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
    chg = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
    do k = 1, m
      cr = cr * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * z
      chg = chg + cr
    end do
  else

    x0 = real ( z, kind = rk )
    if ( x0 < 0.0D+00 ) then
      a = b - a
      a0 = a
      z = - z
    end if

    if ( a < 2.0D+00 ) then
      nl = 0
    else
      nl = 1
      la = int ( a )
      a = a - la - 1.0D+00
    end if

    do n = 0, nl

      if ( 2.0D+00 <= a0 ) then
        a = a + 1.0D+00
      end if

      if ( cdabs ( z ) < 20.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

        chg = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
        crg = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
        do j = 1, 500
          crg = crg * ( a + j - 1.0D+00 ) / ( j * ( b + j - 1.0D+00 ) ) * z
          chg = chg + crg
          if ( abs ( ( chg - chw ) / chg ) < 1.0D-15 ) then
            exit
          end if
          chw = chg
        end do

      else

        call gamma ( a, g1 )
        call gamma ( b, g2 )
        ba = b - a
        call gamma ( ba, g3 )
        cs1 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
        cs2 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
        cr1 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
        cr2 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )

        do i = 1, 8
          cr1 = - cr1 * (     a + i - 1.0D+00 ) * ( a - b + i ) / ( z * i )
          cr2 =   cr2 * ( b - a + i - 1.0D+00 ) * ( i - a ) / ( z * i )
          cs1 = cs1 + cr1
          cs2 = cs2 + cr2
        end do

        x = real ( z, kind = rk )
        y = imag ( z )

        if ( x == 0.0D+00 .and. 0.0D+00 <= y ) then
          phi = 0.5D+00 * pi
        else if ( x == 0.0D+00 .and. y <= 0.0D+00 ) then
          phi = -0.5D+00 * pi
        else
          phi = atan ( y / x )
        end if

        if ( -1.5D+00 * pi < phi .and. phi <= -0.5 * pi ) then
          ns = -1
        else if ( -0.5D+00 * pi < phi .and. phi < 1.5D+00 * pi ) then
          ns = 1
        end if

        if ( y == 0.0D+00 ) then
          cfac = cos ( pi * a )
        else
          cfac = exp ( ns * ci * pi * a )
        end if

        chg1 = g2 / g3 * z ** ( - a ) * cfac * cs1
        chg2 = g2 / g1 * exp ( z ) * z ** ( a - b ) * cs2
        chg = chg1 + chg2

      end if

      if ( n == 0 ) then
        cy0 = chg
      else if ( n == 1 ) then
        cy1 = chg
      end if

    end do

    if ( 2.0D+00 <= a0 ) then
      do i = 1, la - 1
        chg = ( ( 2.0D+00 * a - b + z ) * cy1 + ( b - a ) * cy0 ) / a
        cy0 = cy1
        cy1 = chg
        a = a + 1.0D+00
      end do
    end if

    if ( x0 < 0.0D+00 ) then
      chg = chg * exp ( - z )
    end if

  end if

  !a = a1
  !z = z0

  return
end
end module HyperGeo