module mgh

  use set_precision, only : rk
  
  private

  ! REAL CONSTANTS
  real(kind=rk), parameter, private :: PI = 4.0_rk*atan( 1.0_rk )
  ! real(kind=rk), parameter, private :: macheps   = 1.0e-16_rk
  ! real(kind=rk), parameter, private :: macheps12 = sqrt( macheps )
  ! real(kind=rk), parameter, private :: macheps13 = macheps ** ( 1.0_rk / 3.0_rk )

  ! GLOBAL SCALARS
  integer, private :: global_m    ! Number of equations
  integer, private :: global_n    ! Number of variables
  integer, private :: problem = 0 ! Number of the problem (1 to 35)

  ! PROBLEM SPECIFIC PARAMETER DATA FOR COMPUTING OBJECTIVE FUNCTIONS
  real(kind=rk), parameter, dimension(3), private :: y5 = [ &
       1.5_rk,  &
       2.25_rk, &
       2.625_rk ]

  real(kind=rk), parameter, dimension(15), private :: y8 = [ &
       0.14_rk, &
       0.18_rk, &
       0.22_rk, &
       0.25_rk, &
       0.29_rk, &
       0.32_rk, &
       0.35_rk, &
       0.39_rk, &
       0.37_rk, &
       0.58_rk, &
       0.73_rk, &
       0.96_rk, &
       1.34_rk, &
       2.10_rk, &
       4.39_rk ]

  real(kind=rk), parameter, dimension(15), private :: y9 = [ &
       0.0009_rk, &
       0.0044_rk, &
       0.0175_rk, &
       0.0540_rk, &
       0.1295_rk, &
       0.2420_rk, &
       0.3521_rk, &
       0.3989_rk, &
       0.3521_rk, &
       0.2420_rk, &
       0.1295_rk, &
       0.0540_rk, &
       0.0175_rk, &
       0.0044_rk, &
       0.0009_rk ]

  real(kind=rk), parameter, dimension(16), private :: y10 = [ &
       34780.0_rk, &
       28610.0_rk, &
       23650.0_rk, &
       19630.0_rk, &
       16370.0_rk, &
       13720.0_rk, &
       11540.0_rk, &
       9744.0_rk,  &
       8261.0_rk,  &
       7030.0_rk,  &
       6005.0_rk,  &
       5147.0_rk,  &
       4427.0_rk,  &
       3820.0_rk,  &
       3307.0_rk,  &
       2872.0_rk ]

  real(kind=rk), parameter, dimension(11), private :: y15 = [ &
       0.1957_rk, &
       0.1947_rk, &
       0.1735_rk, &
       0.1600_rk, &
       0.0844_rk, &
       0.0627_rk, &
       0.0456_rk, &
       0.0342_rk, &
       0.0323_rk, &
       0.0235_rk, &
       0.0246_rk ]

  real(kind=rk), parameter, dimension(11), private :: u15 = [ &
       4.0000_rk, &
       2.0000_rk, &
       1.0000_rk, &
       0.5000_rk, &
       0.2500_rk, &
       0.1670_rk, &
       0.1250_rk, &
       0.1000_rk, &
       0.0833_rk, &
       0.0714_rk, &
       0.0625_rk ]

  real(kind=rk), parameter, dimension(33), private :: y17 = [ &
       8.44E-1_rk, &
       9.08E-1_rk, &
       9.32E-1_rk, &
       9.36E-1_rk, &
       9.25E-1_rk, &
       9.08E-1_rk, &
       8.81E-1_rk, &
       8.5E-1_rk,  &
       8.18E-1_rk, &
       7.84E-1_rk, &
       7.51E-1_rk, &
       7.18E-1_rk, &
       6.85E-1_rk, &
       6.58E-1_rk, &
       6.28E-1_rk, &
       6.03E-1_rk, &
       5.8E-1_rk,  &
       5.58E-1_rk, &
       5.38E-1_rk, &
       5.22E-1_rk, &
       5.06E-1_rk, &
       4.9E-1_rk,  &
       4.78E-1_rk, &
       4.67E-1_rk, &
       4.57E-1_rk, &
       4.48E-1_rk, &
       4.38E-1_rk, &
       4.31E-1_rk, &
       4.24E-1_rk, &
       4.2E-1_rk,  &
       4.14E-1_rk, &
       4.11E-1_rk, &
       4.06E-1_rk ]

  real(kind=rk), parameter, dimension(65), private :: y19 = [      &
       1.366e0_rk, 1.191e0_rk, 1.112e0_rk, 1.013e0_rk, 9.91e-1_rk, &
       8.85e-1_rk, 8.31e-1_rk, 8.47e-1_rk, 7.86e-1_rk, 7.25e-1_rk, &
       7.46e-1_rk, 6.79e-1_rk, 6.08e-1_rk, 6.55e-1_rk, 6.16e-1_rk, &
       6.06e-1_rk, 6.02e-1_rk, 6.26e-1_rk, 6.51e-1_rk, 7.24e-1_rk, &
       6.49e-1_rk, 6.49e-1_rk, 6.94e-1_rk, 6.44e-1_rk, 6.24e-1_rk, &
       6.61e-1_rk, 6.12e-1_rk, 5.58e-1_rk, 5.33e-1_rk, 4.95e-1_rk, &
       5.00e-1_rk, 4.23e-1_rk, 3.95e-1_rk, 3.75e-1_rk, 3.72e-1_rk, &
       3.91e-1_rk, 3.96e-1_rk, 4.05e-1_rk, 4.28e-1_rk, 4.29e-1_rk, &
       5.23e-1_rk, 5.62e-1_rk, 6.07e-1_rk, 6.53e-1_rk, 6.72e-1_rk, &
       7.08e-1_rk, 6.33e-1_rk, 6.68e-1_rk, 6.45e-1_rk, 6.32e-1_rk, &
       5.91e-1_rk, 5.59e-1_rk, 5.97e-1_rk, 6.25e-1_rk, 7.39e-1_rk, &
       7.10e-1_rk, 7.29e-1_rk, 7.20e-1_rk, 6.36e-1_rk, 5.81e-1_rk, &
       4.28e-1_rk, 2.92e-1_rk, 1.62e-1_rk, 9.8e-2_rk, 5.4e-2_rk ]

  public :: mgh_evalf, mgh_evalg, mgh_evalh, mgh_evalt, mgh_get_dims, &
            mgh_get_name, mgh_get_x0, mgh_set_dims, mgh_set_problem, rk
  
contains

  ! ------------------------------------------------------------------

  subroutine mgh_get_name( name )
    ! SCALAR ARGUMENT
    character(len=60), intent(out) :: name

    ! This subroutine:
    ! Returns a string with the problem name

    select case ( problem )

    case ( 1 )
       name = "Rosenbrock"

    case ( 2 )
       name = "Freudenstein and Roth"

    case ( 3 )
       name = "Powell badly scaled"

    case ( 4 )
       name = "Brown badly scaled"

    case ( 5 )
       name = "Beale"

    case ( 6 )
       name = "Jennrich and Sampson"

    case ( 7 )
       name = "Helical valley"

    case ( 8 )
       name = "Bard"

    case ( 9 )
       name = "Gaussian"

    case ( 10 )
       name = "Meyer"

    case ( 11 )
       name = "Gulf research and development"

    case ( 12 )
       name = "Box three-dimensional"

    case ( 13 )
       name = "Powell singular"

    case ( 14 )
       name = "Wood"

    case ( 15 )
       name = "Kowalik and Osborne"

    case ( 16 )
       name = "Brown and Dennis"

    case ( 17 )
       name = "Osborne 1"

    case ( 18 )
       name = "Biggs EXP6"

    case ( 19 )
       name = "Osborne 2"

    case ( 20 )
       name = "Watson"

    case ( 21 )
       name = "Extended Rosenbrock"

    case ( 22 )
       name = "Extended Powell singular"

    case ( 23 )
       name = "Penalty I"

    case ( 24 )
       name = "Penalty II"

    case ( 25 )
       name = "Variably dimensioned"

    case ( 26 )
       name = "Trigonometric"

    case ( 27 )
       name = "Brown almost-linear"

    case ( 28 )
       name = "Discrete boundary value"

    case ( 29 )
       name = "Discrete integral equation"

    case ( 30 )
       name = "Broyden tridiagonal"

    case ( 31 )
       name = "Broyden banded"

    case ( 32 )
       name = "Linear - full rank"

    case ( 33 )
       name = "Linear - rank 1"

    case ( 34 )
       name = "Linear - rank 1 with zero columns and rows"

    case ( 35 )
       name = "Chebyquad"

    end select
    
  end subroutine mgh_get_name
  
  ! ------------------------------------------------------------------

  subroutine mgh_set_problem( user_problem, flag )
    ! SCALAR ARGUMENT
    integer, intent(in)  :: user_problem
    integer, intent(out) :: flag

    ! Set problem to compute objective function and derivatives.
    ! Set the global variable **problem**.
    ! If problem is not between 1 and 35, returns flag = - 1.

    ! LOCAL SCALAR
    integer :: n, m

    if ( user_problem .ge. 1 .and. user_problem .le. 35 ) then
       problem = user_problem
       call mgh_get_default_dims( n, m )
       flag = 0
    else
       flag = - 1
    end if

  end subroutine mgh_set_problem

  ! ------------------------------------------------------------------

  subroutine mgh_get_default_dims( n, m )
    ! SCALAR ARGUMENTS
    integer, intent(out), optional :: m, n

    ! This subroutine set and return the default number of variable
    ! and/or equations for each problem

    ! LOCAL SCALARS
    integer :: local_m, local_n

    select case ( problem )

    case ( 1 )
       local_n = 2
       local_m = 2

    case ( 2 )
       local_n = 2
       local_m = 2

    case ( 3 )
       local_n = 2
       local_m = 2

    case ( 4 )
       local_n = 2
       local_m = 3

    case ( 5 )
       local_n = 2
       local_m = 3

    case ( 6 )
       local_n = 2
       local_m = 10

    case ( 7 )
       local_n = 3
       local_m = 3

    case ( 8 )
       local_n = 3
       local_m = 15

    case ( 9 )
       local_n = 3
       local_m = 15

    case ( 10 )
       local_n = 3
       local_m = 16

    case ( 11 )
       local_n = 3
       local_m = 99

    case ( 12 )
       local_n = 3
       local_m = 10

    case ( 13 )
       local_n = 4
       local_m = 4

    case ( 14 )
       local_n = 4
       local_m = 6

    case ( 15 )
       local_n = 4
       local_m = 11

    case ( 16 )
       local_n = 4
       local_m = 20

    case ( 17 )
       local_n = 5
       local_m = 33

    case ( 18 )
       local_n = 6
       local_m = 13

    case ( 19 )
       local_n = 11
       local_m = 65

    case ( 20 )
       local_n = 6
       local_m = 31

    case ( 21 )
       local_n = 10
       local_m = local_n

    case ( 22 )
       local_n = 12
       local_m = local_n

    case ( 23 )
       local_n = 4
       local_m = local_n + 1

    case ( 24 )
       local_n = 4
       local_m = 2*local_n

    case ( 25 )
       local_n = 10
       local_m = local_n + 2

    case ( 26 )
       local_n = 10
       local_m = local_n

    case ( 27 )
       local_n = 40
       local_m = local_n

    case ( 28 )
       local_n = 10
       local_m = local_n

    case ( 29 )
       local_n = 10
       local_m = local_n

    case ( 30 )
       local_n = 10
       local_m = local_n

    case ( 31 )
       local_n = 10
       local_m = local_n

    case ( 32 )
       local_n = 10
       local_m = 10

    case ( 33 )
       local_n = 10
       local_m = 10

    case ( 34 )
       local_n = 10
       local_m = 10

    case ( 35 )
       local_n = 8
       local_m = local_n

    case default
       local_n = 0
       local_m = 0

    end select

    if ( present( m ) ) then
       m        = local_m
       global_m = local_m
    end if
    
    if ( present( n ) ) then
       n        = local_n
       global_n = local_n
    end if

  end subroutine mgh_get_default_dims

  ! ------------------------------------------------------------------

  subroutine mgh_get_dims( n, m )
    ! SCALAR ARGUMENTS
    integer, intent(out), optional :: m, n

    if ( present( m ) ) m = global_m
    if ( present( n ) ) n = global_n
  end subroutine mgh_get_dims

  ! ------------------------------------------------------------------

  subroutine mgh_set_dims( n, m, flag )
    ! SCALAR ARGUMENTS
    integer, intent(in),  optional :: m, n
    integer, intent(out), optional :: flag

    ! This subroutine:
    !
    ! (1) Set global variables ** global_n ** and ** global_m ** from
    !     user input values
    !
    ! (2) Validate user input values.
    !     If    n is  not valid, flag = -1
    !     If    m is  not valid, flag = -2
    !     If both are not valid, flag = -3

    select case ( problem )

    case ( 6 )
       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 11 )
       if ( present ( m ) ) then
          if ( m .ge. global_n .and. m .le. 100 ) then
             global_m = m
          end if
       end if

    case ( 12 )
       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 16 )
       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 18 )
       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 20 )
       if ( present( n ) ) then
          if ( n .ge. 2 .and. n .le. 31 ) then
             global_n = n
          end if
       end if

    case ( 21 )
       if ( present( n ) ) then
          if ( n .gt. 0 .and. modulo( n, 2 ) .eq. 0 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 22 )
       if ( present( n ) ) then
          if ( n .gt. 0 .and. modulo( n, 4 ) .eq. 0 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 23 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n + 1
          end if
       end if

    case ( 24 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = 2*global_n
          end if
       end if

    case ( 25 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n + 2
          end if
       end if

    case ( 26 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 27 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 28 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 29 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 30 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 31 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n
             global_m = global_n
          end if
       end if

    case ( 32 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n

             if ( .not. present( m ) ) global_m = n
          end if
       end if

       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 33 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n

             if ( .not. present( m ) ) global_m = n
          end if
       end if

       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 34 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n

             if ( .not. present( m ) ) global_m = n
          end if
       end if

       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    case ( 35 )
       if ( present( n ) ) then
          if ( n .ge. 1 ) then
             global_n = n

             if ( .not. present( m ) ) global_m = n
          end if
       end if

       if ( present( m ) ) then
          if ( m .ge. global_n ) then
             global_m = m
          end if
       end if

    end select

    ! Validate input n and m
    flag = 0

    if ( present( n ) ) then
       if ( n .ne. global_n ) flag = - 1
    end if
    
    if ( present( m ) ) then
       if ( m .ne. global_m ) flag = flag - 2
    end if

  end subroutine mgh_set_dims
  
  ! ------------------------------------------------------------------

  subroutine mgh_get_x0( x0, factor )
    ! SCALAR ARGUMENT
    real(kind=rk), intent(in), optional :: factor
    
    ! ARRAY ARGUMENT
    real(kind=rk), dimension(global_n), intent(out) :: x0

    ! This subroutine:
    ! Returns the initial point for each problem, scaled by constant
    ! ** factor **

    ! LOCAL SCALAR
    integer       :: j
    real(kind=rk) :: h

    select case ( problem )
    case ( 1 )
       x0 = [ - 1.2_rk, 1.0_rk ]

    case ( 2 )
       x0 = [ 0.5_rk, - 2.0_rk ]

    case ( 3 )
       x0 = [ 0.0_rk, 1.0_rk ]

    case ( 4 )
       x0 = [ 1.0_rk, 1.0_rk ]

    case ( 5 )
       x0 = [ 1.0_rk, 1.0_rk ]

    case ( 6 )
       x0 = [ 0.3_rk, 0.4_rk ]

    case ( 7 )
       x0 = [ - 1.0_rk, 0.0_rk, 0.0_rk ]

    case ( 8 )
       x0 = [ 1.0_rk, 1.0_rk, 1.0_rk ]

    case ( 9 )
       x0 = [ 0.4_rk, 1.0_rk, 0.0_rk ]

    case ( 10 )
       x0 = [ 2.0e-02_rk, 4.0e+03_rk, 2.5e+02_rk ]

    case ( 11 )
       x0 = [ 5.0_rk, 2.5_rk, 1.5e-01_rk ]

    case ( 12 )
       x0 = [ 0.0_rk, 10.0_rk, 20.0_rk ]

    case ( 13 )
       x0 = [ 3.0_rk, - 1.0_rk, 0.0_rk, 1.0_rk ]

    case ( 14 )
       x0 = [ - 3.0_rk, - 1.0_rk, - 3.0_rk, - 1.0_rk ]

    case ( 15 )
       x0 = [ 0.25_rk, 0.39_rk, 0.415_rk, 0.39_rk ]

    case ( 16 )
       x0 = [ 25.0_rk, 5.0_rk, - 5.0_rk, - 1.0_rk ]

    case ( 17 )
       x0 = [ 0.5_rk, 1.5_rk, - 1.0_rk, 0.01_rk, 0.02_rk ]

    case ( 18 )
       x0 = [ 1.0_rk, 2.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk ]

    case ( 19 )
       x0 = [ 1.3_rk, 0.65_rk, 0.65_rk, 0.7_rk, 0.6_rk, 3.0_rk, &
            5.0_rk, 7.0_rk, 2.0_rk, 4.5_rk, 5.5_rk ]

    case ( 20 )
       x0 = 0.0_rk

    case ( 21 )
       do j = 1, global_n, 2
          x0(j  ) = - 1.2_rk
          x0(j+1) =   1.0_rk
       end do

    case ( 22 )
       do j = 1, global_n, 4
          x0(j  ) =   3.0_rk
          x0(j+1) = - 1.0_rk
          x0(j+2) =   0.0_rk
          x0(j+3) =   1.0_rk
       end do

    case ( 23 )
       x0 = [ ( j, j = 1, global_n ) ]

    case ( 24 )
       x0 = 0.5_rk

    case ( 25 )
       x0 = [ ( 1.0_rk - ( real( j, rk ) / real( global_n, rk ) ), j = 1, global_n ) ]

    case ( 26 )
       x0 = 1.0_rk / real( global_n, rk )

    case ( 27 )
       x0 = 0.5_rk

    case ( 28 )
       h = 1.0_rk / ( real( global_n, rk ) + 1.0_rk )
       x0 = [ ( j*h*( j*h - 1.0_rk ), j = 1, global_n ) ]

    case ( 29 )
       h = 1.0_rk / ( real( global_n, rk ) + 1.0_rk )
       x0 = [ ( j*h*( j*h - 1.0_rk ), j = 1, global_n ) ]

    case ( 30 )
       x0 = - 1.0_rk

    case ( 31 )
       x0 = - 1.0_rk

    case ( 32 )
       x0 = 1.0_rk

    case ( 33 )
       x0 = 1.0_rk

    case ( 34 )
       x0 = 1.0_rk

    case ( 35 )
       x0 = [ ( real( j, rk ) / ( real( global_n, rk ) + 1.0_rk ), j = 1, global_n ) ]

    end select

    if ( present( factor ) .and. factor .ne. 1.0_rk ) then
       if ( problem .ne. 20 ) then
          x0(1:global_n) = factor*x0(1:global_n)
       else
          x0(1:global_n) = factor
       end if
    end if

  end subroutine mgh_get_x0
  
  ! ------------------------------------------------------------------
  
  subroutine mgh_evalf( x, f, flag )
    ! SCALAR ARGUMENTS
    integer,       intent(out) :: flag
    real(kind=rk), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=rk), dimension(:), intent(in) :: x

    ! This subroutine computes the objective function value for
    !                  \sum_{i=1}^n f_i^2(x)
    !
    ! Returning flag
    !  0    successful
    ! -1    problem is not between 1 and 35
    ! -3    division by zero
    
    ! LOCAL SCALARS
    integer       :: i, j
    real(kind=rk) :: arg, d1, d2, r, s1, s2, s3, t, t1, t2, t3, t4, &
                     th, tpi

    ! LOCAL ARRAYS
    real(kind=rk), dimension(global_n)   :: fvec, w1
    real(kind=rk), dimension(global_n+1) :: w2
    
    flag = 0

    select case ( problem )

    case ( 1 ) ! Rosenbrock
       f = ( 10.0_rk*( x(2) - x(1)**2 ) )**2 + ( 1.0_rk - x(1) )**2

    case ( 2 )  ! Freudenstein and Roth
       t1 = - 13.0_rk + x(1) + ( ( 5.0_rk - x(2) )*x(2) - 2.0_rk )*x(2)
       t2 = - 29.0_rk + x(1) + ( ( x(2) + 1.0_rk )*x(2) - 14.0_rk )*x(2)
       f = t1**2 + t2**2

    case ( 3 ) ! Powell badly scaled
       t1 = 1e+4_rk*x(1)*x(2) - 1.0_rk
       s1 = exp( -x(1) )
       s2 = exp( -x(2) )
       t2 = s1 + s2 - 1.0001_rk
       f = t1**2 + t2**2

    case ( 4 ) ! Brown badly-scaled
       t1 = x(1) - 1.0e+6_rk
       t2 = x(2) - 2.0e-6_rk
       t3 = x(1)*x(2) - 2.0_rk
       f = t1**2 + t2**2 + t3**2

    case ( 5 ) ! Beale
       f = 0.0_rk
       do i = 1, 3
          s1 = 1.0_rk - x(2)**i
          t1 = y5(i) - x(1)*s1
          f = f + t1**2
       end do

    case ( 6 ) ! Jennrich and Sampson
       t1 = exp( x(1) )
       t2 = exp( x(2) )
       d1 = t1
       d2 = t2
       f = 0.0_rk
       do i = 1, global_m
          f = f + ( 2.0_rk + 2.0_rk*i - ( d1 + d2 ) )**2
          d1 = d1*t1
          d2 = d2*t2
       end do

    case ( 7 ) ! Helical valley
       tpi = 8.0_rk*atan( 1.0_rk )
       th = sign( 0.25_rk, x(2) )

       if ( x(1) .gt. 0.0_rk ) th = atan( x(2)/x(1) )/tpi
       if ( x(1) .lt. 0.0_rk ) th = atan( x(2)/x(1) )/tpi + 0.5_rk

       r = sqrt( x(1)**2 + x(2)**2 )
       t = x(3) - 10.0_rk*th
       f = 100.0_rk*( t**2 + ( r - 1.0_rk )**2 ) + x(3)**2

    case ( 8 ) ! Bard
       f = 0.0_rk
       do i = 1, global_m
          t1 = real( i, rk )
          t2 = 16.0_rk - real( i, rk )
          t3 = min( t1, t2 )

          d1 = t2*x(2) + t3*x(3)
          if ( d1 .ne. 0.0_rk ) then 
             f = f + ( y8(i) - ( x(1) + ( t1 / d1 ) ) )**2
          else
             f = huge( rk )
             flag = - 3
             return
          end if
       end do
       
    case ( 9 ) ! Gaussian
       f = 0.0_rk
       do i = 1, global_m
          d1 = 0.5_rk*real( i-1, rk )
          d2 = 3.5_rk - d1 - x(3)
          r = exp( -0.5_rk*x(2)*d2**2 )
          t = x(1)*r - y9(i)
          f = f + t**2
       end do

    case ( 10 ) ! Meyer
       f = 0.0_rk
       do i = 1, global_m
          t1 = 4.5e+01_rk + 5.0_rk*i

          t2 = t1 + x(3)

          if ( t2 .ne. 0.0_rk ) then
             f = f + ( x(1)*exp( x(2) / t2 ) - y10(i) )**2
          else
             f = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 11 ) ! Gulf research and development
       f = 0.0_rk
       d1 = 2.0_rk/3.0_rk
       do i = 1, global_m
          arg = real( i, rk )/1.0e+2_rk
          r = abs( ( -50.0_rk*log( arg ) )**d1 + 25.0_rk - x(2) )
          t1 = r**x(3)/x(1)
          t2 = exp( -t1 )
          t = t2 - arg
          f = f + t**2
       end do

    case ( 12 ) ! Box three-dimensional
       f = 0.0_rk
       do i = 1, global_m
          d1 = real( i, rk )
          d2 = d1/10.0_rk
          s1 = exp( -d2*x(1) )
          s2 = exp( -d2*x(2) )
          s3 = exp( -d2 ) - exp( -d1 )
          t = s1 - s2 - s3*x(3)
          f = f + t**2
       end do

    case ( 13 ) ! Powell singular
       t1 = x(1) + 10.0_rk*x(2)
       t2 = sqrt( 5.0_rk )*( x(3) - x(4) )
       t3 = ( x(2) - 2.0_rk*x(3) ) ** 2
       t4 = sqrt( 10.0_rk )*( x(1) - x(4) ) ** 2
       f = t1**2 + t2**2 + t3**2 + t4**2

    case ( 14 ) ! Wood
       s1 = x(2) - x(1)**2
       s2 = 1.0_rk - x(1)
       s3 = x(2) - 1.0_rk
       t1 = x(4) - x(3)**2
       t2 = 1.0_rk - x(3)
       t3 = x(4) - 1.0_rk
       f = 100.0_rk*s1**2 + s2**2 + 90.0_rk*t1**2 + t2**2 &
            + 10.0_rk*( s3 + t3 )**2 + ( s3 - t3 )**2/10.0_rk

    case ( 15 ) ! Kowalik and Osborne
       f = 0.0_rk
       do i = 1, global_m
          t1 = u15(i) ** 2 + u15(i)*x(3) + x(4)

          if ( t1 .ne. 0.0_rk ) then
             f = f + ( y15(i) - ( x(1)*u15(i)*( u15(i) + x(2) ) ) / t1 )**2
          else
             f = huge( rk )
             flag = - 3
             return
          end if
       end do
       
    case ( 16 ) ! Brown and Dennis
       f = 0.0_rk
       do i = 1, global_m
          d1 = real( i, rk )/5.0_rk
          d2 = sin( d1 )
          t1 = x(1) + d1*x(2) - exp( d1 )
          t2 = x(3) + d2*x(4) - cos( d1 )
          t = t1**2 + t2**2
          f = f + t**2
       end do

    case ( 17 ) ! Osborne 1
       f = 0.0_rk
       do i = 1, global_m
          t1 = 10.0_rk*real( i - 1, rk )
          f = f + ( y17(i) - ( x(1) + x(2)*exp( - t1*x(4) ) + x(3)*exp( - t1*x(5) ) ) )**2
       end do

    case ( 18 ) ! Biggs EXP6
       f = 0.0_rk
       do i = 1, global_m
          d1 = real( i, rk )/10.0_rk
          d2 = exp( -d1 ) - 5.0_rk*exp( -10.0_rk*D1 ) + 3.0_rk*exp( -4.0_rk*d1 )
          s1 = exp( -d1*x(1) )
          s2 = exp( -d1*x(2) )
          s3 = exp( -d1*x(5) )
          t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
          f = f + t**2
       end do

    case ( 19 ) ! Osborne 2
       f = 0.0_rk
       do i = 1, global_m
          t1 = real( i - 1, rk ) / 10.0_rk
          f = f + ( y19(i) - ( x(1)*exp( - t1*x(5) ) + x(2)*exp( - ( t1 - x(9) )**2*x(6) ) + &
               x(3)*exp( - ( t1 - x(10) )**2*x(7) ) + x(4)*exp( - ( t1 - x(11) )**2*x(8) ) ) )**2
       end do

    case ( 20 ) ! Watson
       f = 0.0_rk
       do i = 1, 29
          d1 = real( i, rk )/29.0_rk
          s1 = 0.0_rk
          d2 = 1.0_rk
          do j = 2, global_n
             s1 = s1 + real( j-1, rk )*d2*x(j)
             d2 = d1*d2
          end do
          s2 = 0.0_rk
          d2 = 1.0_rk
          do j = 1, global_n
             s2 = s2 + d2*x(j)
             d2 = d1*d2
          end do
          t = s1 - s2**2 - 1.0_rk
          f = f + t**2
       end do
       t1 = x(2) - x(1)**2 - 1.0_rk
       f = f + x(1)**2 + t1**2

    case ( 21 ) ! Extended Rosenbrock
       f = 0.0_rk
       do j = 1, global_n, 2
          t1 = 1.0_rk - x(j)
          t2 = 10.0_rk*( x(j+1) - x(j)**2 )
          f = f + t1**2 + t2**2
       end do

    case ( 22 ) ! Extended Powell singular
       f = 0.0_rk
       do j = 1, global_n, 4
          t = x(j) + 10.0_rk*x(j+1)
          t1 = x(j+2) - x(j+3)
          s1 = 5.0_rk*t1
          t2 = x(j+1) - 2.0_rk*x(j+2)
          s2 = t2**3
          t3 = x(j) - x(j+3)
          s3 = 10.0_rk*t3**3
          f = f + t**2 + s1*t1 + s2*t2 + s3*t3
       end do

    case ( 23 ) ! Penalty I
       t1 = -0.25_rk
       t2 = 0.0_rk
       do j = 1, global_n
          t1 = t1 + x(j)**2
          t2 = t2 + ( x(j) - 1.0_rk )**2
       end do
       f = 1.0e-05_rk*t2 + t1**2

    case ( 24 ) ! Penalty II
       t1 = -1.0_rk
       t2 = 0.0_rk
       t3 = 0.0_rk
       d1 = exp( 0.1_rk )
       d2 = 1.0_rk
       s2 = 0.0_rk
       do j = 1, global_n
          t1 = t1 + real( global_n-j+1, rk )*x(j)**2
          s1 = exp( x(j)/10.0_rk )
          if ( j .ne. 1 ) then
             s3 = s1 + s2 - d2*( d1 + 1.0_rk )
             t2 = t2 + s3**2
             t3 = t3 + ( s1 - 1.0_rk/d1 )**2
          end if
          s2 = s1
          d2 = d1*d2
       end do
       f = 1.0e-05_rk*( t2 + t3 ) + t1**2 + ( x(1) - 0.2_rk )**2

    case ( 25 ) ! Variably dimensioned
       t1 = 0.0_rk
       t2 = 0.0_rk
       do j = 1, global_n
          t1 = t1 + real( j, rk )*( x(j) - 1.0_rk )
          t2 = t2 + ( x(j) - 1.0_rk )**2
       end do
       f = t2 + t1**2*( 1.0_rk + t1**2 )

    case ( 26 ) ! Trigonometric
       s1 = 0.0_rk
       do j = 1, global_n
          s1 = s1 + cos( x(j) )
       end do
       f = 0.0_rk
       do j = 1, global_n
          t = real( global_n+j, rk ) - sin( x(j) ) - s1 - real( j, rk )*cos( x(j) )
          f = f + t**2
       end do

    case ( 27 ) ! Brown almost linear
       s1 = sum( x(1:global_n) )
       f = 0.0_rk
       do i = 1, global_n-1
          t = x(i) + s1 - ( global_n + 1.0_rk )
          f = f + t**2
       end do
       f = f + ( product( x(1:global_n) ) - 1.0_rk )**2

    case ( 28 ) ! Discrete boundary
       d1 = 1.0_rk/( real( global_n, rk ) + 1.0_rk )
       
       f = ( 2.0_rk*x(1) - x(2) + 0.5_rk*d1**2*( x(1) + d1 + 1.0_rk )**3 )**2
       do i = 2, global_n-1
          d2 = real( i, rk )*d1
          f = f + &
               ( 2.0_rk*x(i) - x(i-1) - x(i+1) + 0.5_rk*d1**2*( x(i) + d2 + 1.0_rk )**3 )**2
       end do
       d2 = real( global_n, rk )*d1
       f = f + &
            ( 2.0_rk*x(global_n) - x(global_n-1) + 0.5_rk*d1**2*( x(global_n) + d2 + 1.0_rk )**3 )**2

    case ( 29 ) ! Discrete integral equation
       d1 = 1.0_rk/( real( global_n, rk ) + 1.0_rk )
       f = 0.0_rk

       w1(1)          = d1*( x(1) + d1 + 1.0_rk )**3
       w2(global_n)   = ( 1.0_rk - global_n*d1 )*( x(global_n) + global_n*d1 + 1.0_rk )**3
       w2(global_n+1) = 0.0_rk
       do i = 2, global_n
          t1 = real( i, rk )*d1
          t2 = real( global_n-i+1, rk )*d1
          w1(i)            = w1(i-1)          +              t1*( x(i)            + t1 + 1.0_rk )**3
          w2(global_n-i+1) = w2(global_n-i+2) + ( 1.0_rk - t2 )*( x(global_n-i+1) + t2 + 1.0_rk )**3
       end do
       
       do i = 1, global_n
          t1 = real( i, rk )*d1
          f = f + ( x(i) + 0.5_rk*d1*( ( 1.0_rk - t1 )*w1(i) + t1*w2(i+1) ) )**2
       end do

    case ( 30 ) ! Broyden tridiagonal
       f = ( ( 3.0_rk - 2.0_rk*x(1) )*x(1) - 2.0_rk*x(2) + 1.0_rk )**2
       do i = 2, global_m-1
          f = f + ( ( 3.0_rk - 2.0_rk*x(i) )*x(i) - x(i-1) - 2.0_rk*x(i+1) + 1.0_rk )**2
       end do
       f = f + ( ( 3.0_rk - 2.0_rk*x(global_n) )*x(global_n) - x(global_n-1) + 1.0_rk )**2

    case ( 31 ) ! Broyden banded
       f = 0.0_rk
       do i = 1, global_m
          s1 = x(i)*( 2.0_rk + 5.0_rk*x(i) ** 2 ) + 1.0_rk
          do j = max( 1, i - 5 ), min( global_n, i + 1 )
             if ( j .ne. i ) then
                s1 = s1 - x(j)*( 1.0_rk + x(j) )
             end if
          end do
          f = f + s1**2
       end do

    case ( 32 ) ! Linear Function - full rank
       s1 = sum( x(1:global_n) )
       d1 = 2.0_rk / real( global_m, rk )
       f = 0.0_rk
       do i = 1, global_n
          t = x(i) - d1*s1 - 1.0_rk
          f = f + t**2
       end do
       do i = global_n+1, global_m
          t = - d1*s1 - 1.0_rk
          f = f + t**2
       end do

    case ( 33 ) ! Linear function - rank 1
       s1 = 0.0_rk
       do j = 1, global_n
          s1 = s1 + j*x(j)
       end do
       f = 0.0_rk
       do i = 1, global_m
          t = i*s1 - 1.0_rk
          f = f + t**2
       end do

    case ( 34 ) ! Linear function - rank 1 with zero columns and rows
       s1 = 0.0_rk
       do j = 2, global_n-1
          s1 = s1 + j*x(j)
       end do
       f = 2.0_rk
       do i = 2, global_m-1
          t = ( i - 1 )*s1 - 1.0_rk
          f = f + t**2
       end do

    case ( 35 ) ! Chebyquad
       fvec(1:global_m) = 0.0_rk
       do j = 1, global_n
          t1 = 1.0_rk
          t2 = 2.0_rk*x(j) - 1.0_rk
          t = 2.0_rk*t2
          do i = 1, global_n
             fvec(i) = fvec(i) + t2
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       f = 0.0_rk
       d1 = 1.0_rk/real( global_n, rk )
       do i = 1, global_n
          t = d1*fvec(i)
          if ( modulo( i, 2 ) .eq. 0 ) t = t + 1.0_rk/( real( i, rk )**2 - 1.0_rk )
          f = f + t**2
       end do

    case default
       flag = - 1

    end select

  end subroutine mgh_evalf

  ! ------------------------------------------------------------------

  subroutine mgh_evalg( x, g, flag )
    ! SCALAR ARGUMENTS
    integer, intent(out) :: flag
    
    ! ARRAY ARGUMENTS
    real(kind=rk), dimension(:), intent(in)  :: x
    real(kind=rk), dimension(:), intent(out) :: g

    ! This subroutine computes the first derivative vector for the function
    !                      \sum_{i=1}^n f_i^2(x)
    !
    ! Returning flag
    !  0    successful
    ! -1    problem is not between 1 and 35
    ! -2    i is not between 1 and m
    ! -3    division by zero

    ! LOCAL SCALARS
    integer       :: i, j
    real(kind=rk) :: arg, d1, d2, d3, r, s1, s2, s3, t, t1, t2, t3, &
                     t4, th, tpi

    ! LOCAL WORKING ARRAYS
    real(kind=rk), dimension(global_n)   :: w1
    real(kind=rk), dimension(global_n+1) :: w2
    real(kind=rk), dimension(global_m)   :: w3

    flag = 0

    select case ( problem )

    case ( 1 ) ! Rosenbrock
       t1 = 10.0_rk*( x(2) - x(1)**2 )
       t2 = 1.0_rk - x(1)

       g(1) = 2.0_rk*( - 20.0_rk*x(1)*t1 - t2 )
       g(2) = 20.0_rk*t1

    case ( 2 ) ! Freudenstein and Roth
       t1 = - 13.0_rk + x(1) + ( ( 5.0_rk - x(2) )*x(2) - 2.0_rk )*x(2)
       t2 = - 29.0_rk + x(1) + ( ( x(2) + 1.0_rk )*x(2) - 14.0_rk )*x(2)

       g(1) = 2.0_rk*( t1 + t2 )
       g(2) = 2.0_rk*( ( ( 10.0_rk - 3.0_rk*x(2) )*x(2) - 2.0_rk )*t1 + &
                       ( ( 3.0_rk*x(2) + 2.0_rk )*x(2) - 14.0_rk )*t2 )
       
    case ( 3 ) ! Powell badly scaled
       t1 = 1.0e+4*x(1)*x(2) - 1.0_rk
       s1 = exp( -x(1) )
       s2 = exp( -x(2) )
       t2 = s1 + s2 - 1.0001_rk
       g(1) = 2.0_rk*( 1.0e+4_rk*x(2)*t1 - s1*t2 )
       g(2) = 2.0_rk*( 1.0e+4_rk*x(1)*t1 - s2*t2 )

    case ( 4 ) ! Brown badly-scaled
       t1 = x(1) - 1.0e+6_rk
       t2 = x(2) - 2.0e-6_rk
       t3 = x(1)*x(2) - 2.0_rk
       g(1) = 2.0_rk*( t1 + x(2)*t3 )
       g(2) = 2.0_rk*( t2 + x(1)*t3 )
       
    case ( 5 ) ! Beale
       g(1:2) = 0.0_rk
       do i = 1, 3
          t1 = y5(i) - x(1)*( 1.0_rk - x(2)**i )

          g(1) = g(1) + 2.0_rk*t1*( x(2)**i - 1.0_rk )
          g(2) = g(2) + 2.0_rk*t1*( i*x(1)*x(2)**(i-1) )
       end do

    case ( 6 ) ! Jennrich and Sampson
       g(1:global_n) = 0.0_rk
       t1 = exp( x(1) )
       t2 = exp( x(2) )
       d1 = t1
       d2 = t2
       do i = 1, global_m
          s1 = 2.0_rk + 2.0_rk*i - ( d1 + d2 )
          g(1) = g(1) - 2.0_rk*i*s1*d1
          g(2) = g(2) - 2.0_rk*i*s1*d2
          d1 = d1*t1
          d2 = d2*t2
       end do

    case ( 7 ) ! Helical valley
       tpi = 8.0_rk*atan( 1.0_rk )
       th = sign( 0.25_rk, x(2) )
       if ( x(1) .gt. 0.0_rk ) th = atan( x(2)/x(1) )/tpi
       if ( x(1) .lt. 0.0_rk ) th = atan( x(2)/x(1) )/tpi + 0.5_rk
       arg = x(1)**2 + x(2)**2
       r = sqrt( arg )
       t = x(3) - 10.0_rk*th
       s1 = 10.0_rk*t/(tpi*arg)
       g(1) = 200.0_rk*( x(1) - x(1)/r + x(2)*s1 )
       g(2) = 200.0_rk*( x(2) - x(2)/r - x(1)*s1 )
       g(3) = 2.0_rk*( 100.0_rk*t + x(3) )

    case ( 8 ) ! Bard
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          t1 = real( i, rk )
          t2 = 16.0_rk - real( i, rk )
          t3 = min( t1, t2 )

          d1 = t2*x(2) + t3*x(3)
          d2 = d1**2
          if ( d1 .ne. 0.0_rk ) then 
             s1 = y8(i) - ( x(1) + ( t1 / d1 ) )
             s2 = t1 / d2
          else
             g(1:3) = huge( rk )
             flag = - 3
             return
          end if

          g(1) = g(1) - 2.0_rk*s1
          g(2) = g(2) + 2.0_rk*s1*s2*t2
          g(3) = g(3) + 2.0_rk*s1*s2*t3
       end do

    case ( 9 ) ! Gaussian
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          d1 = 0.5_rk*real( i-1, rk )
          d2 = 3.5_rk - d1 - x(3)
          arg = -0.5_rk*x(2)*d2**2
          r = exp( arg )
          t = x(1)*r - y9(i)
          s1 = r*t
          s2 = d2*s1
          g(1) = g(1) + s1
          g(2) = g(2) - d2*s2
          g(3) = g(3) + s2
       end do
       g(1) = 2.0_rk*g(1)
       g(2) = x(1)*g(2)
       g(3) = 2.0_rk*x(1)*x(2)*g(3)

    case ( 10 ) ! Meyer
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          t1 = 4.5e+01_rk + 5.0_rk*i
          t2 = t1 + x(3)
          if ( t2 .ne. 0.0_rk ) then
             t3 = exp( x(2) / t2 )
             s1 = x(1)*t3 - y10(i)
             g(1) = g(1) + 2.0_rk*s1*t3
             g(2) = g(2) + 2.0_rk*s1*( t3*x(1)/t2 )
             g(3) = g(3) - 2.0_rk*s1*( t3*x(1)*x(2)/( t2**2 ) )
          else
             g(1:3) = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 11 ) ! Gulf research and development
       g(1:global_n) = 0.0_rk
       d1 = 2.0_rk/3.0_rk
       do i = 1, global_m
          arg = real( i, rk )/1.0e+2_rk
          r = abs( ( -50.0_rk*log( arg ) )**d1 + 25.0_rk - x(2) )
          t1 = r**x(3)/x(1)
          t2 = exp( -t1 )
          t = t2 - arg
          s1 = t1*t2*t
          g(1) = g(1) + s1
          g(2) = g(2) + s1/r
          g(3) = g(3) - s1*log( r )
       end do
       g(1) = 2.0_rk*g(1)/x(1)
       g(2) = 2.0_rk*x(3)*g(2)
       g(3) = 2.0_rk*g(3)

    case ( 12 ) ! Box three-dimensional
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          d1 = real( i, rk )
          d2 = d1/10.0_rk
          s1 = exp( -d2*x(1) )
          s2 = exp( -d2*x(2) )
          s3 = exp( -d2 ) - exp( -d1 )
          t = s1 - s2 - s3*x(3)
          th = d2*t
          g(1) = g(1) - s1*th
          g(2) = g(2) + s2*th
          g(3) = g(3) - s3*t
       end do
       g(1:global_n) = 2.0_rk*g(1:global_n)

    case ( 13 ) ! Powell singular
       t1 = x(1) + 10.0_rk*x(2)
       t2 = sqrt( 5.0_rk )*( x(3) - x(4) )
       t3 = ( x(2) - 2.0_rk*x(3) ) ** 2
       t4 = sqrt( 10.0_rk )*( x(1) - x(4) ) ** 2

       d1 = sqrt( 5.0_rk )
       d2 = 2.0_rk*( x(2) - 2.0_rk*x(3) )
       d3 = 2.0_rk*sqrt( 10.0_rk )*( x(1) - x(4) )

       g(1) = 2.0_rk*( t1 + t4*d3 )
       g(2) = 2.0_rk*( 10.0_rk*t1 + t3*d2 )
       g(3) = 2.0_rk*( t2*d1 - 2.0_rk*t3*d2 )
       g(4) = 2.0_rk*( - t2*d1 - t4*d3 )

    case ( 14 ) ! Wood
       s1 = x(2) - x(1)**2
       s2 = 1.0_rk - x(1)
       s3 = x(2) - 1.0_rk
       t1 = x(4) - x(3)**2
       t2 = 1.0_rk - x(3)
       t3 = x(4) - 1.0_rk
       g(1) = -2.0_rk*( 2.0e+2_rk*x(1)*s1 + s2 )
       g(2) = 2.0e+2_rk*s1 + 20.2_rk*s3 + 19.8_rk*t3
       g(3) = -2.0_rk*( 1.8e+2_rk*x(3)*t1 + t2 )
       g(4) = 1.8e+2_rk*t1 + 20.2_rk*t3 + 19.8_rk*s3

    case ( 15 ) ! Kowalik and Osborne
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          t1 = u15(i)**2 + u15(i)*x(3) + x(4)

          if ( t1 .ne. 0.0_rk ) then
             d1 = y15(i) - ( x(1)*u15(i)*( u15(i) + x(2) ) ) / t1
             g(1) = g(1) - 2.0_rk*d1*u15(i)*( u15(i) + x(2) )/t1
             g(2) = g(2) - 2.0_rk*d1*x(1)*u15(i)/t1
             g(3) = g(3) + 2.0_rk*d1*x(1)*u15(i)**2*( u15(i) + x(2) )/t1**2
             g(4) = g(4) + 2.0_rk*d1*x(1)*u15(i)*( u15(i) + x(2) )/t1**2
          else
             g(1:4) = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 16 ) ! Brown and Dennis
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          d1 = real( i, rk )/5.0_rk
          s1 = exp( d1 )
          s2 = sin( d1 )
          s3 = cos( d1 )
          t = ( x(1) + d1*x(2) - s1 )**2 + ( x(3) + x(4)*s2 - s3 )**2
          g(1) = g(1) + 4.0_rk*t*( x(1) + d1*x(2) - s1 )
          g(2) = g(2) + 4.0_rk*t*( d1*( x(1) + d1*x(2) - s1 ) )
          g(3) = g(3) + 4.0_rk*t*( x(3) + s2*x(4) - s3 )
          g(4) = g(4) + 4.0_rk*t*( s2*( x(3) + s2*x(4) - s3 ) )
       end do

    case ( 17 ) ! Osborne 1
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          t1 = 10.0_rk*real( i - 1, rk )
          d1 = exp( - t1*x(4) )
          d2 = exp( - t1*x(5) )
          s1 = y17(i) - ( x(1) + x(2)*d1 + x(3)*d2 )
          g(1) = g(1) - 2.0_rk*s1
          g(2) = g(2) - 2.0_rk*s1*d1
          g(3) = g(3) - 2.0_rk*s1*d2
          g(4) = g(4) + 2.0_rk*s1*t1*x(2)*d1
          g(5) = g(5) + 2.0_rk*s1*t1*x(3)*d2
       end do

    case ( 18 ) ! Biggs EXP6
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          d1 = real( i, rk )/10.0_rk
          d2 = exp( -d1 ) - 5.0_rk*exp( -10.0_rk*d1 ) + 3.0_rk*exp( -4.0_rk*d1 )
          s1 = exp( -d1*x(1) )
          s2 = exp( -d1*x(2) )
          s3 = exp( -d1*x(5) )
          t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
          th = d1*t
          g(1) = g(1) - s1*th
          g(2) = g(2) + s2*th
          g(3) = g(3) + s1*t
          g(4) = g(4) - s2*t
          g(5) = g(5) - s3*th
          g(6) = g(6) + s3*t
       end do
       g(1) = 2.0_rk*x(3)*g(1)
       g(2) = 2.0_rk*x(4)*g(2)
       g(3) = 2.0_rk*g(3)
       g(4) = 2.0_rk*g(4)
       g(5) = 2.0_rk*x(6)*g(5)
       g(6) = 2.0_rk*g(6)

    case ( 19 ) ! Osborne 2
       g(1:global_n) = 0.0_rk
       do i = 1, global_m
          t1 = real( i - 1, rk ) / 10.0_rk
          s1 = y19(i) - ( x(1)*exp( - t1*x(5) ) + x(2)*exp( - ( t1 - x(9) )**2*x(6) ) + &
               x(3)*exp( - ( t1 - x(10) )**2*x(7) ) + x(4)*exp( - ( t1 - x(11) )**2*x(8) ) )
          g( 1) = g( 1) - 2.0_rk*s1*exp( - t1*x(5) )
          g( 2) = g( 2) - 2.0_rk*s1*exp( - ( t1 - x( 9) )**2*x(6) )
          g( 3) = g( 3) - 2.0_rk*s1*exp( - ( t1 - x(10) )**2*x(7) )
          g( 4) = g( 4) - 2.0_rk*s1*exp( - ( t1 - x(11) )**2*x(8) )
          g( 5) = g( 5) + 2.0_rk*s1*t1*x(1)*exp( - t1*x(5) )
          g( 6) = g( 6) + 2.0_rk*s1*( t1 - x( 9) )**2*x(2)*exp( - ( t1 - x( 9) )**2*x(6) )
          g( 7) = g( 7) + 2.0_rk*s1*( t1 - x(10) )**2*x(3)*exp( - ( t1 - x(10) )**2*x(7) )
          g( 8) = g( 8) + 2.0_rk*s1*( t1 - x(11) )**2*x(4)*exp( - ( t1 - x(11) )**2*x(8) )
          g( 9) = g( 9) - 4.0_rk*s1*( t1 - x( 9) )*x(6)*x(2)*exp( - ( t1 - x( 9) )**2*x(6) )
          g(10) = g(10) - 4.0_rk*s1*( t1 - x(10) )*x(7)*x(3)*exp( - ( t1 - x(10) )**2*x(7) )
          g(11) = g(11) - 4.0_rk*s1*( t1 - x(11) )*x(8)*x(4)*exp( - ( t1 - x(11) )**2*x(8) )
       end do

    case ( 20 ) ! Watson
       g(1:global_n) = 0.0_rk
       do i = 1, 29
          d1 = real( i, rk )/29.0_rk
          s1 = 0.0_rk
          d2 = 1.0_rk
          do j = 2, global_n
             s1 = s1 + real( j-1, rk )*d2*x(j)
             d2 = d1*d2
          end do
          s2 = 0.0_rk
          d2 = 1.0_rk
          do j = 1, global_n
             s2 = s2 + d2*x(j)
             d2 = d1*d2
          end do
          t = s1 - s2**2 - 1.0_rk
          s3 = 2.0_rk*d1*s2
          d2 = 2.0_rk/d1
          do j = 1, global_n
             g(j) = g(j) + d2*( real( j-1, rk ) - s3 )*t
             d2 = d1*d2
          end do
       end do
       t1 = x(2) - x(1)**2 - 1.0_rk
       g(1) = g(1) + x(1)*( 2.0_rk - 4.0_rk*t1 )
       g(2) = g(2) + 2.0_rk*t1

    case ( 21 ) ! Extended Rosenbrock
       do j = 1, global_n, 2
          t1 = 1.0_rk - x(j)
          g(j+1) = 2.0e+2_rk*( x(j+1) - x(j)**2 )
          g(j)   = -2.0_rk*( x(j)*g(j+1) + t1 )
       end do

    case ( 22 ) ! Extended Powell singular
       do j = 1, global_n, 4
          t = x(j) + 10.0_rk*x(j+1)
          t1 = x(j+2) - x(j+3)
          s1 = 5.0_rk*t1
          t2 = x(j+1) - 2.0_rk*x(j+2)
          s2 = 4.0_rk*t2**3
          t3 = x(j) - x(j+3)
          s3 = 20.0_rk*t3**3
          g(j) = 2.0_rk*( t + s3 )
          g(j+1) = 20.0_rk*t + s2
          g(j+2) = 2.0_rk*( s1 - s2 )
          g(j+3) = -2.0_rk*( s1 + s3 )
       end do

    case ( 23 ) ! Penalty I
       t1 = sum( x(1:global_n)**2 ) - 0.25_rk
       d1 = 2.0e-5_rk
       th = 4.0_rk*t1
       g(1:global_n) = d1*( x(1:global_n) - 1.0_rk ) + x(1:global_n)*th

    case ( 24 ) ! Penalty II
       t1 = - 1.0_rk
       do j = 1, global_n
          t1 = t1 + real( global_n-j+1, rk )*x(j)**2
       end do
       d1 = exp( 0.1_rk )
       d2 = 1.0_rk
       s2 = 0.0_rk
       th = 4.0_rk*t1
       do j = 1, global_n
          g(j) = real( global_n-j+1, rk )*x(j)*th
          s1 = exp(x(j)/10.0_rk)
          if ( j .ne. 1 ) then
             s3 = s1 + s2 - d2*( d1 + 1.0_rk )
             g(j) = g(j) + 1.0e-5_rk*s1*( s3 + s1 - 1.0_rk/d1 )/5.0_rk
             g(j-1) = g(j-1) + 1.0e-5_rk*s2*s3/5.0_rk
          end if
          s2 = s1
          d2 = d1*d2
       end do
       g(1) = g(1) + 2.0_rk*( x(1) - 0.2_rk )

    case ( 25 ) ! Variably dimensioned
       t1 = 0.0_rk
       do j = 1, global_n
          t1 = t1 + real( j, rk )*( x(j) - 1.0_rk )
       end do
       t = t1*( 1.0_rk + 2.0_rk*t1**2 )
       do j = 1, global_n
          g(j) = 2.0_rk*( x(j) - 1.0_rk + real( j, rk )*t )
       end do

    case ( 26 ) ! Trigonometric
       s1 = 0.0_rk
       do j = 1, global_n
          g(j) = cos( x(j) )
          s1 = s1 + g(j)
       end do
       s2 = 0.0_rk
       do j = 1, global_n
          th = sin( x(j) )
          t = real( global_n+j, rk ) - th - s1 - real( j, rk )*g(j)
          s2 = s2 + t
          g(j) = ( real( j, rk )*th - g(j) )*t
       end do
       g(1:global_n) = 2.0_rk*( g(1:global_n) + sin( x(1:global_n) )*s2 )

    case ( 27 ) ! Brown almost linear
       g(1:global_n) = 0.0_rk

       s1 = sum( x(1:global_n) )
       do i = 1, global_n-1
          t = x(i) + s1 - ( global_n + 1 )

          g(1:global_n) = g(1:global_n) + 2.0_rk*t
          g(i) = g(i) + 2.0_rk*t
       end do

       w1(1) = 1.0_rk
       w2(global_n) = 1.0_rk
       do i = 2, global_n
          w1(i) = w1(i-1)*x(i-1)
          w2(global_n-i+1) = w2(global_n-i+2)*x(global_n-i+2)
       end do

       g(1:global_n) = g(1:global_n) + 2.0_rk        &
            *( w1(global_n)*x(global_n) - 1.0_rk ) &
            *w1(1:global_n)*w2(1:global_n)

    case ( 28 ) ! Discrete boundary value
       d1 = 1.0_rk / ( real( global_n, rk ) + 1.0_rk )

       ! For m = 1
       d2 = d1
       t  = 2.0_rk*x(1) - x(2) + d1**2*( x(1) + d2 + 1.0_rk )**3/2.0_rk
       g(1) =   2.0_rk*t*( 2.0_rk + 1.5_rk*d1**2*( x(1) + d2 + 1.0_rk )**2 )
       g(2) = - 2.0_rk*t
       
       do i = 2, global_n-1
          d2 = real( i, rk )*d1
          t  = 2.0_rk*x(i) - x(i-1) - x(i+1) + d1**2*( x(i) + d2 + 1.0_rk )**3/2.0_rk
          g(i-1) = g(i-1) - 2.0_rk*t
          g(i)   = g(i)   + 2.0_rk*t*( 2.0_rk + 1.5_rk*d1**2*( x(i) + d2 + 1.0_rk )**2 )
          g(i+1) =        - 2.0_rk*t
       end do

       ! For m = n
       d2 = real( global_n, rk )*d1
       t  = 2.0_rk*x(global_n) - x(global_n-1) + d1**2*( x(global_n) + d2 + 1.0_rk )**3/2.0_rk
       g(global_n-1) = g(global_n-1) - 2.0_rk*t
       g(global_n)   = g(global_n)   + 2.0_rk*t*( 2.0_rk + 1.5_rk*d1**2*( x(global_n) + d2 + 1.0_rk )**2 )

    case ( 29 ) ! Discrete integral equation
       g(1:global_n) = 0.0_rk
       
       d1 = 1.0_rk/( real( global_n, rk ) + 1.0_rk )
       w1(1)          = d1*( x(1) + d1 + 1.0_rk )**3
       w2(global_n)   = ( 1.0_rk - global_n*d1 )*( x(global_n) + global_n*d1 + 1.0_rk )**3
       w2(global_n+1) = 0.0_rk
       do i = 2, global_n
          t1 = real( i, rk )*d1
          t2 = real( global_n-i+1, rk )*d1
          w1(i)            = w1(i-1)          +              t1*( x(i)            + t1 + 1.0_rk )**3
          w2(global_n-i+1) = w2(global_n-i+2) + ( 1.0_rk - t2 )*( x(global_n-i+1) + t2 + 1.0_rk )**3
       end do
       
       do i = 1, global_n
          t1 = real( i, rk )*d1
          t  = x(i) + 0.5_rk*d1*( ( 1.0_rk - t1 )*w1(i) + t1*w2(i+1) )

          do j = 1, i
             t2 = real( j, rk )*d1
             g(j) = g(j) + 2.0_rk*t*( 1.5_rk*d1*( 1.0_rk - t1 )*t2*( x(j) + t2 + 1.0_rk ) ** 2 )
          end do

          g(i) = g(i) + 2.0_rk*t

          do j = i+1, global_n
             t2 = real( j, rk )*d1
             g(j) = g(j) + 2.0_rk*t*( 1.5_rk*d1*( 1.0_rk - t2 )*t1*( x(j) + t2 + 1.0_rk ) ** 2 )
          end do
       end do

    case ( 30 ) ! Broyden tridiagonal
       ! For i = 1
       t = ( 3.0_rk - 2.0_rk*x(1) )*x(1) - 2.0_rk*x(2) + 1.0_rk
       g(1) =   2.0_rk*t*( 3.0_rk - 4.0_rk*x(1) )
       g(2) = - 4.0_rk*t

       do i = 2, global_n-1
          t = ( 3.0_rk - 2.0_rk*x(i) )*x(i) - x(i-1) - 2.0_rk*x(i+1) + 1.0_rk
          g(i-1) = g(i-1) - 2.0_rk*t
          g(i)   = g(i)   + 2.0_rk*t*( 3.0_rk - 4.0_rk*x(i) )
          g(i+1) = - 4.0_rk*t
       end do

       ! For i = n
       t = ( 3.0_rk - 2.0_rk*x(global_n) )*x(global_n) - x(global_n-1) + 1.0_rk
       g(global_n-1) = g(global_n-1) - 2.0_rk*t
       g(global_n)   = g(global_n)   + 2.0_rk*t*( 3.0_rk - 4.0_rk*x(global_n) )

    case ( 31 ) ! Broyden banded
       g(1:global_n) = 0.0_rk
       do i = 1, global_n
          s1 = 0.0_rk
          do j = max( 1, i-5 ), i-1
             s1 = s1 + x(j)*( 1.0_rk + x(j) )
          end do
          if ( i .ne. global_n ) s1 = s1 + x(i+1)*( 1.0_rk + x(i+1) )
          
          t = x(i)*( 2.0_rk + 5.0_rk*x(i)**2 ) + 1.0_rk - s1
          do j = max( 1, i-5 ), i-1
             g(j) = g(j) - 2.0_rk*t*( 1.0_rk + 2.0_rk*x(j) )
          end do
          g(i) = g(i) + 2.0_rk*t*( 2.0_rk + 15.0_rk*x(i)**2 )
          if ( i .ne. global_n ) g(i+1) = g(i+1) - 2.0_rk*t*( 1.0_rk + 2.0_rk*x(i+1) )
       end do

    case ( 32 ) ! Linear function - full rank
       arg = 2.0_rk/real( global_m, rk )
       s1  = sum( x(1:global_n) )

       w3(         1:global_n) = x(1:global_n) - arg*s1 - 1.0_rk
       w3(global_n+1:global_m) =               - arg*s1 - 1.0_rk

       s2 = sum( w3(1:global_m) )
       g(1:global_n) = 2.0_rk*( w3(1:global_n) - arg*s2 )
       
    case ( 33 ) ! Linear function - rank 1
       s1 = 0.0_rk
       do i = 1, global_n
          s1 = s1 + i*x(i)
       end do

       w3(1:global_m) = [ ( i*( i*s1 - 1.0_rk ), i=1, global_m ) ]
       
       d1 = 2.0_rk*sum( w3(1:global_m) )
       g(1:global_n) = [ ( d1*i, i=1, global_n ) ]
  
    case ( 34 ) ! Linear function - rank 1 with zero columns and rows
       s1 = 0.0_rk
       do i = 2, global_n-1
          s1 = s1 + i*x(i)
       end do

       w3(2:global_n-1) = [ ( ( i-1.0_rk )*( ( i-1.0_rk )*s1 - 1.0_rk ), i=2, global_n-1 ) ]

       d1 = 2.0_rk*sum( w3(2:global_n-1) )
  
       g(1) = 0.0_rk
       g(global_n) = 0.0_rk
       g(2:global_n-1) = [ ( d1*i, i=2, global_n-1 ) ]

    case ( 35 ) ! Chebyquad
       w1(1:global_n) = 0.0_rk
       do j = 1, global_n
          t1 = 1.0_rk
          t2 = 2.0_rk*x(j) - 1.0_rk
          t = 2.0_rk*t2
          do i = 1, global_n
             w1(i) = w1(i) + t2
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       d1 = 1.0_rk/real( global_n, rk )
       do i = 1, global_n
          w1(i) = d1*w1(i)
          if ( modulo( i, 2 ) .eq. 0 ) then
             w1(i) = w1(i) + 1.0_rk/( real( i, rk )**2 - 1.0_rk )
          end if
       end do

       do j = 1, global_n
          g(j) = 0.0_rk
          t1 = 1.0_rk
          t2 = 2.0_rk*x(j) - 1.0_rk
          t = 2.0_rk*t2
          s1 = 0.0_rk
          s2 = 2.0_rk
          do i = 1, global_n
             g(j) = g(j) + w1(i)*s2
             th = 4.0_rk*t2 + t*s2 - s1
             s1 = s2
             s2 = th
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       d2 = 2.0_rk*d1
       g(1:global_n) = d2*g(1:global_n)

    case default
       flag = - 1

    end select
    
  end subroutine mgh_evalg
  
  ! ------------------------------------------------------------------

  subroutine mgh_evalh( x, h, flag )
    ! SCALAR ARGUMENT
    integer, intent(out) :: flag
    
    ! ARRAY ARGUMENTS
    real(kind=rk), dimension(:),   intent(in)  :: x
    real(kind=rk), dimension(:,:), intent(out) :: h

    ! This subroutine computes the second derivative matrix for the function
    !                      \sum_{i=1}^n f_i^2(x)
    ! Returning flag
    !  0    successful
    ! -1    problem is not between 1 and 35
    ! -3    division by zero

    ! LOCAL SCALARS
    integer       :: i, j, k, l
    real(kind=rk) :: arg, d1, d2, d3, d4, logr, p1, p2, piarg,       &
                     piarg2, r, r3inv, s1, s2, s3, s1s2, s1s3, s2s3, &
                     ss1, ss2, t, t1, t2, t3, t4, th, tt, tt1, tt2,  &
                     tth

    ! LOCAL WORKING ARRAYS
    real(kind=rk), dimension(global_n)                :: fvec, gvec, w1
    real(kind=rk), dimension(global_n+1)              :: w2
    real(kind=rk), dimension(0:global_n,1:global_n+1) :: prod

    flag = 0

    select case ( problem )

    case ( 1 ) ! Rosenbrock
       t1 = 10.0_rk*( x(2) - x(1)**2 )
       h(1,1) = 2.0_rk*( 400.0_rk*x(1)**2 - 20.0_rk*t1 + 1.0_rk )
       h(1,2) = - 400.0_rk*x(1)
       h(2,2) = 200.0_rk

    case ( 2 ) ! Freudenstein and Roth
       t1 = - 13.0_rk + x(1) + ( ( 5.0_rk - x(2) )*x(2) - 2.0_rk )*x(2)
       t2 = - 29.0_rk + x(1) + ( ( x(2) + 1.0_rk )*x(2) - 14.0_rk )*x(2)
       h(1,1) = 4.0_rk
       h(1,2) = 8.0_rk*( 3.0_rk*x(2) - 4.0_rk )
       h(2,2) = 2.0_rk*(                                                                      &
            t1*( - 6.0_rk*x(2) + 10.0_rk ) + ( ( 10.0_rk - 3.0_rk*x(2) )*x(2) - 2.0_rk )**2 + &
            t2*(   6.0_rk*x(2) +  2.0_rk ) + ( ( 3.0_rk*x(2) + 2.0_rk )*x(2) - 14.0_rk )**2 )

    case ( 3 ) ! Powell badly scaled
       t1 = 1.0e+04_rk*x(1)*x(2) - 1
       t2 = exp( - x(1) ) + exp( - x(2) ) - 1.0001_rk
       h(1,1) = 2.0_rk*( 1.0e+8_rk*x(2)**2 + t2*exp( - x(1) ) + exp( - x(1) )**2 )
       h(1,2) = 2.0_rk*( t1*1.0e+04_rk + 1.0e+8_rk*x(1)*x(2) + exp( - x(1) - x(2) ) )
       h(2,2) = 2.0_rk*( 1.0e+8_rk*x(1)**2 + t2*exp( - x(2) ) + exp( - x(2) )**2 )

    case ( 4 ) ! Brown badly scaled
       h(1,1) = 2.0_rk*( 1.0_rk + x(2)**2 )
       h(1,2) = 4.0_rk*( x(1)*x(2) - 1.0_rk )
       h(2,2) = 2.0_rk*( 1.0_rk + x(1)**2 )

    case ( 5 ) ! Beale
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       
       do i = 1, global_m
          t1 = y5(i) - x(1)*( 1.0_rk - x(2)**i )
          h(1,1) = h(1,1) + 2.0_rk*( x(2)**i - 1.0_rk )**2
          h(1,2) = h(1,2) + 2.0_rk*( t1*i*x(2)**(i-1) + i*x(1)*x(2)**(i-1)*( x(2)**i - 1.0_rk ) )
          h(2,2) = h(2,2) + 2.0_rk*( t1*( i - 1 )*i*x(1)*x(2)**( i - 2 ) + ( i*x(1)*x(2)**( i - 1 ) )**2 )
       end do

    case ( 6 ) ! Jennrich and Sampson
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          d1 = exp( i*x(1) )
          d2 = exp( i*x(2) )
          t1 = 2.0_rk + 2.0_rk*i - ( d1 + d2 )
          h(1,1) = h(1,1) + 2.0_rk*( (i*d1)**2 - t1*i**2*d1 )
          h(1,2) = h(1,2) + 2.0_rk*i**2*d1*d2
          h(2,2) = h(2,2) + 2.0_rk*( (i*d2)**2 - t1*i**2*d2 )
       end do

    case ( 7 ) ! Helical Valley
       if (x(1) .eq. 0.0_rk) then
          th = sign( 0.25_rk, x(2) )
       else
          th = atan( x(2)/x(1) ) / ( 2.0_rk*pi )
          if (x(1) .lt. 0.0_rk) th = th + 0.5_rk
       end if
       arg = x(1)**2 + x(2)**2
       piarg = pi*arg
       piarg2 = piarg*arg
       r3inv = 1.0_rk / sqrt( arg )**3
       t = x(3) - 10.0_rk*th
       s1 = 5.0_rk*t / piarg
       p1 = 2.0e+3_rk*x(1)*x(2)*t / piarg2
       p2 = ( 5.0_rk/piarg )**2

       h(1,1) = 2.0e+2_rk - 2.0e+2_rk*(r3inv-p2)*x(2)**2 - p1
       h(1,2) = 2.0e+2_rk*x(1)*x(2)*r3inv + &
            1.0e+3_rk/piarg2*( t*(x(1)**2-x(2)**2) - 5.0_rk*x(1)*x(2)/pi )
       h(2,2) = 2.0e+2_rk - 2.0e+2_rk*(r3inv-p2)*x(1)**2 + p1
       h(1,3) =  1.0e+3_rk*x(2) / piarg
       h(2,3) = -1.0e+3_rk*x(1) / piarg
       h(3,3) = 202.0_rk

    case ( 8 ) ! Bard
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       
       do i = 1, global_m
          d1 = real( i, rk )
          d2 = 16.0_rk - real( i, rk )
          d3 = min( d1, d2 )

          s1 = d2*x(2) + d3*x(3)
          if ( s1 .ne. 0.0_rk ) then
             t1 = y8(i) - ( x(1) + ( d1 / s1 ) )
             t2 = d1 / s1**2
             t3 = - 2.0_rk*d1 / s1**3
             s2 = t1*t3 + t2**2
             h(1,1) = h(1,1) + 2.0_rk
             h(1,2) = h(1,2) - 2.0*d2*t2
             h(2,2) = h(2,2) + 2.0_rk*s2*d2**2
             h(1,3) = h(1,3) - 2.0_rk*t2*d3
             h(2,3) = h(2,3) + 2.0_rk*s2*d2*d3
             h(3,3) = h(3,3) + 2.0_rk*s2*d3**2
          else
             h(1:global_n,1:global_n) = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 9 ) ! Gaussian
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       do i = 1, global_m
          d1 = 0.5_rk*real( i-1, rk )
          d2 = 3.5_rk - d1 - x(3)
          arg = 0.5_rk*x(2)*d2**2
          r = exp( - arg )
          t = x(1)*r - y9(i)
          t1 = 2.0_rk*x(1)*r - y9(i)
          h(1,1) = h(1,1) + r**2
          h(1,2) = h(1,2) - r*t1*d2**2
          h(2,2) = h(2,2) + r*t1*d2**4
          h(1,3) = h(1,3) + d2*r*t1
          h(2,3) = h(2,3) + d2*r*( t - arg*t1 )
          h(3,3) = h(3,3) + r*( x(2)*t1*d2**2 - t )
       end do
       h(1,1) = 2.0_rk*h(1,1)
       h(2,2) = 0.5_rk*x(1)*h(2,2)
       h(1,3) = 2.0_rk*x(2)*h(1,3)
       h(2,3) = 2.0_rk*x(1)*h(2,3)
       h(3,3) = 2.0_rk*x(1)*x(2)*h(3,3)

    case ( 10 ) ! Meyer
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          d1 = 4.5e+01_rk + 5.0_rk*real( i, rk )
          d2 = d1 + x(3)
          if ( d2 .ne. 0.0_rk ) then
             t1 = exp( x(2) / d2 )
             t2 = x(1)*t1 - y10(i)
             s1 = t2 + t1*x(1)
             h(1,1) = h(1,1) + 2.0_rk*t1**2
             h(1,2) = h(1,2) + 2.0_rk*s1*t1/d2
             h(2,2) = h(2,2) + 2.0_rk*t1*s1*x(1)/d2**2
             h(1,3) = h(1,3) - 2.0_rk*t1*s1*x(2)/d2**2
             h(2,3) = h(2,3) - 2.0_rk*t1*x(1)/d2**2*( t2 + s1*x(2)/d2 )
             h(3,3) = h(3,3) + 2.0_rk*t1*x(1)*x(2)/d2**3*( 2.0_rk*t2 + s1*x(2)/d2 )
          else
             h(1:global_n,1:global_n) = 0.0_rk
             flag = - 3
             return
          end if
       end do

    case ( 11 ) ! Gulf Research and Development
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       d1 = 2.0_rk/3.0_rk
       do i = 1, global_m
          arg = real( i, rk )/1.0e+2_rk
          r = ( -50.0_rk*log( arg ) )**d1 + 25.0_rk - x(2)
          t1 = abs( r )**x(3)/x(1)
          t2 = exp( -t1 )
          t3 = t1*t2*( t1*t2 + ( t1 - 1.0_rk )*( t2 - arg ) )
          t = t1*t2*( t2 - arg )
          logr = log( abs( r ) )
          h(1,1) = h(1,1) + t3 - t
          h(1,2) = h(1,2) + t3/r
          h(2,2) = h(2,2) + ( t + x(3)*t3 )/r**2
          h(1,3) = h(1,3) - t3*logr
          h(2,3) = h(2,3) + (t-x(3)*t3*logr)/r
          h(3,3) = h(3,3) + t3*logr**2
       end do

       h(1,1) = h(1,1) / x(1)**2
       h(1,2) = h(1,2)*x(3)/x(1)
       h(2,2) = h(2,2)*x(3)
       h(1,3) = h(1,3) / x(1)

       do j = 1, global_n
          do i = 1, j
             h(i,j) = 2.0_rk*h(i,j)
          end do
       end do

    case ( 12 ) ! Box-three dimensional
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          d1 = 0.1_rk*real( i, rk )
          t1 = exp( - d1*x(1) )
          t2 = exp( - d1*x(2) )
          t3 = exp( - 10.0_rk*d1 )
          t4 = exp( - d1 )
          s1 = t1 - t2 - x(3)*( t4 - t3 )
          h(1,1) = h(1,1) + 2.0_rk*t1*d1**2*( s1 + t1 )
          h(1,2) = h(1,2) - 2.0_rk*t1*t2*d1**2
          h(2,2) = h(2,2) + 2.0_rk*t2*d1**2*( t2 - s1 )
          h(1,3) = h(1,3) - 2.0_rk*d1*t1*( t3 - t4 )
          h(2,3) = h(2,3) + 2.0_rk*d1*t2*( t3 - t4 )
          h(3,3) = h(3,3) + 2.0_rk*( t3 - t4 )**2
       end do

    case ( 13 ) ! Powell singular
       d3 = ( x(2) - 2.0_rk*x(3) )**2
       d4 = sqrt( 10.0_rk )*( x(1) - x(4) )**2
       t1 = 2.0_rk*( x(2) - 2.0_rk*x(3) )
       t2 = 2.0_rk*sqrt( 10.0_rk )*( x(1) - x(4) )
       t3 = 2.0_rk*sqrt( 10.0_rk )
       h(1,1) = 2.0_rk*( 1.0_rk + d4*t3 + t2**2 )
       h(1,2) = 20.0_rk
       h(2,2) = 2.0_rk*( 100.0_rk + 2.0_rk*d3 + t1**2 )
       h(1,3) = 0.0_rk
       h(2,3) = - 4.0_rk*( 2.0_rk*d3 + t1**2 )
       h(3,3) = 8.0_rk*( 1.25_rk + 2.0_rk*d3 + t1**2 )
       h(1,4) = - 2.0_rk*( d4*t3 + t2**2 )
       h(2,4) = 0.0_rk
       h(3,4) = - 10.0_rk
       h(4,4) = 2.0_rk*( 5.0_rk + d4*t3 + t2**2 )

    case ( 14 ) ! Wood
       h(1,1) = 1.2e+3_rk*x(1)**2 - 4.0e+2_rk*x(2) + 2.0_rk
       h(2,2) = 2.202e+2_rk
       h(3,3) = 1.08e+3_rk*x(3)**2 - 3.6e+2_rk*x(4) + 2.0_rk
       h(4,4) = 2.002e+2_rk
       h(1,2) = -4.0e+2_rk*x(1)
       h(1,3) = 0.0_rk
       h(2,3) = 0.0_rk
       h(1,4) = 0.0_rk
       h(2,4) = 1.98e+1_rk
       h(3,4) = -3.6e+2_rk*x(3)

    case ( 15 ) ! Kowalik and Osborne
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          s1 = u15(i) + x(2)
          t1 = u15(i)**2 + u15(i)*x(3) + x(4)

          if ( t1 .ne. 0.0_rk ) then
             d1 = y15(i) - x(1)*u15(i)*s1/t1

             h(1,1) = h(1,1) + 2.0_rk*( u15(i)*s1/t1 )**2
             h(1,2) = h(1,2) + 2.0_rk*u15(i)/t1*( x(1)*u15(i)*s1/t1 - d1 )
             h(2,2) = h(2,2) + 2.0_rk*( x(1)*u15(i)/t1 )**2
             h(1,3) = h(1,3) + 2.0_rk*s1*u15(i)**2/t1**2*( d1 - u15(i)*s1*x(1)/t1 )
             h(2,3) = h(2,3) + 2.0_rk*x(1)*u15(i)**2/t1**2*( d1 - x(1)*u15(i)*s1/t1 )
             h(3,3) = h(3,3) + 2.0_rk*u15(i)**3*x(1)*s1/t1**3*( x(1)*u15(i)*s1/t1 - 2.0_rk*d1 )
             h(1,4) = h(1,4) + 2.0_rk*u15(i)*s1/t1**2*( d1 - u15(i)*s1*x(1)/t1 )
             h(2,4) = h(2,4) + 2.0_rk*x(1)*u15(i)/t1**2*( d1 - x(1)*u15(i)*s1/t1 )
             h(3,4) = h(3,4) + 2.0_rk*x(1)*u15(i)**2*s1/t1**3*( x(1)*u15(i)*s1/t1 - 2.0_rk*d1 )
             h(4,4) = h(4,4) + 2.0_rk*x(1)*u15(i)*s1/t1**3*( x(1)*u15(i)*s1/t1 - 2.0_rk*d1 )
          else
             h(1:global_n,1:global_n) = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 16 ) ! Brown and Dennis
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          d1 = real( i, rk )/5.0_rk
          d2 = sin( d1 )
          t1 = x(1) + d1*x(2) - exp( d1 )
          t2 = x(3) + d2*x(4) - cos( d1 )
          t = 8.0_rk * t1 * t2
          s1 = 12.0_rk*t1**2 + 4.0_rk*t2**2
          s2 = 12.0_rk*t2**2 + 4.0_rk*t1**2
          h(1,1) = h(1,1) + s1
          h(1,2) = h(1,2) + s1*d1
          h(2,2) = h(2,2) + s1*d1**2
          h(1,3) = h(1,3) + t
          h(2,3) = h(2,3) + t*d1
          h(3,3) = h(3,3) + s2
          h(1,4) = h(1,4) + t*d2
          h(2,4) = h(2,4) + t*d1*d2
          h(3,4) = h(3,4) + s2*d2
          h(4,4) = h(4,4) + s2*d2**2
       end do

    case ( 17 ) ! Osborne 1
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          t1 = 10.0_rk*real( i-1, rk )
          d1 = exp( - t1*x(4) )
          d2 = exp( - t1*x(5) )
          s1 = y17(i) - ( x(1) + x(2)*d1 + x(3)*d2 )
          h(1,1) = h(1,1) + 2.0_rk
          h(1,2) = h(1,2) + 2.0_rk*d1
          h(2,2) = h(2,2) + 2.0_rk*d1**2
          h(1,3) = h(1,3) + 2.0_rk*d2
          h(2,3) = h(2,3) + 2.0_rk*d1*d2
          h(3,3) = h(3,3) + 2.0_rk*d2**2
          h(1,4) = h(1,4) - 2.0_rk*t1*x(2)*d1
          h(2,4) = h(2,4) + 2.0_rk*t1*d1*( s1 - x(2)*d1 )
          h(3,4) = h(3,4) - 2.0_rk*t1*x(2)*d1*d2
          h(4,4) = h(4,4) + 2.0_rk*t1**2*x(2)*d1*( x(2)*d1 - s1 )
          h(1,5) = h(1,5) - 2.0_rk*t1*x(3)*d2
          h(2,5) = h(2,5) - 2.0_rk*t1*x(3)*d1*d2
          h(3,5) = h(3,5) + 2.0_rk*t1*d2*( s1 - x(3)*d2 )
          h(4,5) = h(4,5) + 2.0_rk*t1**2*d1*d2*x(2)*x(3)
          h(5,5) = h(5,5) + 2.0_rk*t1**2*x(3)*d2*( x(3)*d2 - s1 )
       end do

    case ( 18 ) ! Biggs EXP6
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do i = 1, global_m
          d1 = real( i, rk )/10.0_rk
          d2 = exp( - d1 ) - 5.0_rk*exp( - 10.0_rk*d1 ) + 3.0_rk*exp( - 4.0_rk*d1 )
          s1 = exp( - d1*x(1) )
          s2 = exp( - d1*x(2) )
          s3 = exp( - d1*x(5) )
          t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
          d2 = d1**2
          s1s2 = s1 * s2
          s1s3 = s1 * s3
          s2s3 = s2 * s3
          h(1,1) = h(1,1) + d2*s1*( t + x(3)*s1 )
          h(1,2) = h(1,2) - d2*s1s2
          h(2,2) = h(2,2) - d2*s2*( t - x(4)*s2 )
          h(1,3) = h(1,3) - d1*s1*( t + x(3)*s1 )
          h(2,3) = h(2,3) + d1*s1s2
          h(3,3) = h(3,3) + s1**2
          h(1,4) = h(1,4) + d1*s1s2
          h(2,4) = h(2,4) + d1*s2*( t - x(4)*s2 )
          h(3,4) = h(3,4) - s1s2
          h(4,4) = h(4,4) + s2**2
          h(1,5) = h(1,5) + d2*s1s3
          h(2,5) = h(2,5) - d2*s2s3
          h(3,5) = h(3,5) - d1*s1s3
          h(4,5) = h(4,5) + d1*s2s3
          h(5,5) = h(5,5) + d2*s3*( t + x(6)*s3 )
          h(1,6) = h(1,6) - d1*s1s3
          h(2,6) = h(2,6) + d1*s2s3
          h(3,6) = h(3,6) + s1s3
          h(4,6) = h(4,6) - s2s3
          h(5,6) = h(5,6) - d1*s3*( t + x(6)*s3 )
          h(6,6) = h(6,6) + s3**2
       end do
       h(1,1) = x(3)*h(1,1)
       h(2,2) = x(4)*h(2,2)
       h(5,5) = x(6)*h(5,5)
       h(1,2) = x(3)*x(4)*h(1,2)
       h(2,3) = x(4)*h(2,3)
       h(1,4) = x(3)*h(1,4)
       h(1,5) = x(3)*x(6)*h(1,5)
       h(2,5) = x(4)*x(6)*h(2,5)
       h(3,5) = x(6)*h(3,5)
       h(4,5) = x(6)*h(4,5)
       h(1,6) = x(3)*h(1,6)
       h(2,6) = x(4)*h(2,6)

       do j = 1, global_n
          do i = 1, j
             h(i,j) = 2.0_rk*h(i,j)
          end do
       end do

    case ( 19 ) ! Osborne 2
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       
       do i = 1, global_m
          t1 = real( i - 1, rk ) / 10.0_rk
          d1 = exp( - t1*x(5) )
          d2 = exp( - ( t1 - x(9) )**2*x(6) )
          d3 = exp( - ( t1 - x(10) )**2*x(7) )
          d4 = exp( - ( t1 - x(11) )**2*x(8) )
          s1 = y19(i) - ( x(1)*d1 + x(2)*d2 + x(3)*d3 + x(4)*d4 )

          w1(1:global_n) = - [ d1, d2, d3, d4,       & ! 1, 2, 3, 4
                 - t1*x(1)*d1,                       & ! 5
                 - ( t1 - x(9) )**2*x(2)*d2,         & ! 6
                 - ( t1 - x(10) )**2*x(3)*d3,        & ! 7
                 - ( t1 - x(11) )**2*x(4)*d4,        & ! 8
                 2.0_rk*( t1 - x(9) )*x(6)*x(2)*d2,  & ! 9
                 2.0_rk*( t1 - x(10) )*x(7)*x(3)*d3, & ! 10
                 2.0_rk*( t1 - x(11) )*x(8)*x(4)*d4 ]  ! 11

          do k = 1, global_n
             do j = 1, k
                h(j,k) = h(j,k) + 2.0_rk*w1(j)*w1(k)
             end do
          end do

          h( 1, 5) = h( 1, 5) + 2.0_rk*s1*t1*d1
          h( 5, 5) = h( 5, 5) - 2.0_rk*s1*t1**2*x(1)*d1
          h( 2, 6) = h( 2, 6) + 2.0_rk*s1*( t1 - x(9) )**2*d2
          h( 6, 6) = h( 6, 6) - 2.0_rk*s1*( t1 - x(9) )**4*x(2)*d2
          h( 3, 7) = h( 3, 7) + 2.0_rk*s1*( t1 - x(10) )**2*d3
          h( 7, 7) = h( 7, 7) - 2.0_rk*s1*( t1 - x(10) )**4*x(3)*d3
          h( 4, 8) = h( 4, 8) + 2.0_rk*s1*( t1 - x(11) )**2*d4
          h( 8, 8) = h( 8, 8) - 2.0_rk*s1*( t1 - x(11) )**4*x(4)*d4
          h( 2, 9) = h( 2, 9) - 4.0_rk*s1*( t1 - x(9) )*x(6)*d2
          h( 6, 9) = h( 6, 9) + 4.0_rk*s1*( t1 - x(9) )*x(2)*d2*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          h( 9, 9) = h( 9, 9) - 4.0_rk*s1*x(6)*x(2)*d2*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          h( 3,10) = h( 3,10) - 4.0_rk*s1*( t1 - x(10) )*x(7)*d3
          h( 7,10) = h( 7,10) + 4.0_rk*s1*x(3)*( t1 - x(10) )*d3*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          h(10,10) = h(10,10) - 4.0_rk*s1*x(7)*x(3)*d3*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          h( 4,11) = h( 4,11) - 4.0_rk*s1*( t1 - x(11) )*x(8)*d4
          h( 8,11) = h( 8,11) + 4.0_rk*s1*x(4)*( t1 - x(11) )*d4*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          h(11,11) = h(11,11) - 4.0_rk*s1*x(8)*x(4)*d4*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
       end do

    case ( 20 ) ! Watson
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       
       do i = 1, 29
          d1 = real( i, rk )/29.0_rk
          d2 = 1.0_rk
          s1 = 0.0_rk
          s2 = x(1)
          do j = 2, global_n
             s1 = s1 + real( j-1, rk )*d2*x(j)
             d2 = d1*d2
             s2 = s2 + d2*x(j)
          end do
          t = 2.0_rk*( s1 - s2**2 - 1.0_rk )*d1**2
          s3 = 2.0_rk*d1*s2
          d2 = 1.0_rk/d1
          do j = 1, global_n
             t1 = real( j-1, rk ) - s3
             h(j,j) = h(j,j) + ( t1**2 - t )*d2**2
             d3 = 1.0_rk/d1
             do k = 1, j-1
                h(k,j) = h(k,j) + ( t1*( real( k-1, rk ) - s3 ) - t )*d2*d3
                d3 = d1*d3
             end do
             d2 = d1*d2
          end do
       end do
       t3 = x(2) - x(1)**2 - 1.0_rk
       h(1,1) = h(1,1) + 1.0_rk - 2.0_rk*( t3 - 2.0_rk*x(1)**2 )
       h(2,2) = h(2,2) + 1.0_rk
       h(1,2) = h(1,2) - 2.0_rk*x(1)

       do j = 1, global_n
          do i = 1, j
             h(i,j) = 2.0_rk*h(i,j)
          end do
       end do

    case ( 21 ) ! Extended Rosenbrock
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do j = 1, global_n, 2
          h(j  ,j  ) = 1.2e+3_rk*x(j)**2 - 4.0e+2_rk*x(j+1) + 2.0_rk
          h(j  ,j+1) = -4.0e+2_rk*x(j)
          h(j+1,j+1) = 2.0e+2_rk
       end do

    case ( 22 ) ! Extended Powell singular function
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       do j = 1, global_n, 4
          t2 = x(j+1) - 2.0_rk*x(j+2)
          t3 = x(j) - x(j+3)
          s1 = 12.0_rk*t2**2
          s2 = 120.0_rk*t3**2

          h(j  ,j  ) = 2.0_rk + s2
          h(j  ,j+1) = 2.0_rk*10.0_rk
          h(j+1,j+1) = 2.0e+2_rk + s1

          h(j  ,j+2) = 0.0_rk
          h(j+1,j+2) = -2.0_rk*s1
          h(j+2,j+2) = 10.0_rk + 4.0_rk*s1

          h(j  ,j+3) = -s2
          h(j+1,j+3) = 0.0_rk
          h(j+2,j+3) = -10.0_rk
          h(j+3,j+3) = 10.0_rk + s2
       end do

    case ( 23 ) ! Penalty I
       t1 = -0.25_rk
       do j = 1, global_n
          t1 = t1 + x(j)**2
       end do
       d1 = 2.0e-5_rk
       th = 4.0_rk*t1
       do j = 1, global_n
          do k = 1, j-1
             h(k,j) = 8.0_rk*x(j)*x(k)
          end do
          h(j,j) = d1 + th + 8.0_rk*x(j)**2
          ! h(j,j) = th + 8.0_rk*x(j)**2 - 1.0_rk
       end do
            
    case ( 24 ) ! Penalty II
       t1 = -1.0_rk
       do j = 1, global_n
          t1 = t1 + real( global_n-j+1, rk )*x(j)**2
       end do
       d1 = exp( 0.1_rk )
       d2 = 1.0_rk
       th = 4.0_rk*t1
       s2 = 0.0_rk
       do j = 1, global_n
          h(j,j) = 8.0_rk*( real( global_n-j+1, rk )*x(j) )**2 + real( global_n-j+1, rk )*th
          s1 = exp( x(j)/10.0_rk )
          if ( j .gt. 1 ) then
             s3 = s1 + s2 - d2*( d1 + 1.0_rk )
             h(j  ,j  ) = h(j  ,j  ) + 1.0e-5_rk*s1*( s3 + s1 - 1.0_rk/d1 + 2.0_rk*s1 ) / 50.0_rk
             h(j-1,j-1) = h(j-1,j-1) + 1.0e-5_rk*s2*( s2 + s3 ) / 50.0_rk
             do k = 1, j-1
                t1 = exp( real( k, rk )/10.0_rk )
                h(k,j) = 8.0_rk*real( global_n-j+1 )*real( global_n-k+1, rk )*x(j)*x(k)
             end do
             h(j-1,j) = h(j-1,j) + 1.0e-5_rk*s1*s2/50.0_rk
          end if
          s2 = s1
          d2 = d1*d2
       end do
       h(1,1) = h(1,1) + 2.0_rk
         
    case ( 25 ) ! Variably dimensioned
       t1 = 0.0_rk
       do j = 1, global_n
          t1 = t1 + real( j, rk )*( x(j)-1.0_rk )
       end do
       t = 1.0_rk + 6.0_rk*t1**2
       do j = 1, global_n
          h(j,j) = 2.0_rk + 2.0_rk*t*real( j, rk )**2
          do k = 1, j-1
             h(k,j) = 2.0_rk*t*real( j*k, rk )
          end do
       end do

    case ( 26 ) ! Trigonometric function
       s1 = 0.0_rk
       do j = 1, global_n
          h(j,j) = sin( x(j) )
          s1 = s1 + cos( x(j) )
       end do
       s2 = 0.0_rk
       do j = 1, global_n
          th = cos( x(j) )
          t = real( global_n+j, rk ) - h(j,j) - s1 - real( j, rk )*th
          s2 = s2 + t
          do k = 1, j-1
             h(k,j) = 2.0_rk*( sin( x(k) )*( real( global_n+j+k, rk )*h(j,j)-th ) - &
                  h(j,j)*cos( x(k) ) )
          end do
          h(j,j) = real( j*(j+2)+global_n, rk )*h(j,j)**2 + &
               th*( th-real( 2*j+2, rk )*h(j,j) ) + t*( real( j, rk )*th + h(j,j) )
       end do

       do j = 1, global_n
          h(j,j) = 2.0_rk*( h(j,j) + cos( x(j) )*s2 )
       end do

    case ( 27 ) ! Brown almost linear
       do j = 1, global_n
          prod(j-1,j) = 1.0_rk
          do k = 1, j
             h(k,j) = 0.0_rk
          end do

          do k = j, global_n
             prod(k,j) = prod(k-1,j)*x(k)
          end do
          prod(j,global_n+1) = 1.0_rk
       end do
       
       do i = 1, global_m
          if ( i .eq. global_n ) then
             do j = 1, global_n
                t1 = prod(global_n,1)
                t = t1 - 1.0_rk
                h(j,j) = h(j,j) + 2.0_rk*( prod(j-1,1)*prod(global_n,j+1) )**2
                
                do l = 1, j-1
                   t2 = prod(l-1,1)*prod(j-1,l+1)*prod(global_n,j+1)
                   h(l,j) = h(l,j) + 2.0_rk*t2*( 2.0_rk*t1 - 1.0_rk )
                end do
             end do
          else
             do j = 1, global_n
                do k = 1, j
                   if ( j .eq. i .and. k .eq. i ) then
                      h(k,j) = h(k,j) + 8.0_rk
                   elseif ( j .eq. i .or. k .eq. i ) then
                      h(k,j) = h(k,j) + 4.0_rk
                   else
                      h(k,j) = h(k,j) + 2.0_rk
                   end if
                end do
             end do
          end if
       end do

    case ( 28 ) ! Discrete boundary value
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       
       d1 = 1.0_rk / ( real( global_n, rk ) + 1.0_rk )

       ! For i = 1
       d2 = d1
       arg = x(1) + d2 + 1.0_rk
       t  = 2.0_rk*x(1) - x(2) + d1**2*arg**3/2.0_rk
       t1 = 2.0_rk + 1.5_rk*d1**2*arg**2
       h(1,1) =   2.0_rk*( 3.0_rk*d1**2*arg*t + t1**2 )
       h(1,2) = - 2.0_rk*t1
       h(2,2) =   2.0_rk
       
       do i = 2, global_n-1
          d2 = real( i, rk )*d1
          arg = x(i) + d2 + 1.0_rk
          t  = 2.0_rk*x(i) - x(i-1) - x(i+1) + d1**2*arg**3/2.0_rk
          t1 = 2.0_rk + 1.5_rk*d1**2*arg**2
          h(i-1,i-1) = h(i-1,i-1) + 2.0_rk
          h(i-1,i  ) = h(i-1,i  ) - 2.0_rk*t1
          h(i,  i  ) = h(i,  i  ) + 2.0_rk*( 3.0_rk*d1**2*arg*t + t1**2 )
          h(i-1,i+1) = h(i-1,i+1) + 2.0_rk
          h(i  ,i+1) = h(i  ,i+1) - 2.0_rk*t1
          h(i+1,i+1) = h(i+1,i+1) + 2.0_rk
       end do

       ! For i = n
       d2 = real( global_n, rk )*d1
       arg = x(global_n) + d2 + 1.0_rk
       t  = 2.0_rk*x(global_n) - x(global_n-1) + d1**2*arg**3/2.0_rk
       t1 = 2.0_rk + 1.5_rk*d1**2*arg**2
       h(global_n-1,global_n-1) = h(global_n-1,global_n-1) + 2.0_rk
       h(global_n-1,global_n  ) = h(global_n-1,global_n  ) - 2.0_rk*t1
       h(global_n,  global_n  ) = h(global_n,  global_n  ) + 2.0_rk*( 3.0_rk*d1**2*arg*t + t1**2 )

    case ( 29 ) ! Discrete integral equation
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do
       
       d1 = 1.0_rk/( real( global_n, rk ) + 1.0_rk )
       w1(1)          = d1*( x(1) + d1 + 1.0_rk )**3
       w2(global_n)   = ( 1.0_rk - global_n*d1 )*( x(global_n) + global_n*d1 + 1.0_rk )**3
       w2(global_n+1) = 0.0_rk
       do i = 2, global_n
          t1 = real( i, rk )*d1
          t2 = real( global_n-i+1, rk )*d1
          w1(i)            = w1(i-1)          +              t1*( x(i)            + t1 + 1.0_rk )**3
          w2(global_n-i+1) = w2(global_n-i+2) + ( 1.0_rk - t2 )*( x(global_n-i+1) + t2 + 1.0_rk )**3
       end do

       do i = 1, global_n
          t1 = real( i, rk )*d1
          t  = x(i) + 0.5_rk*d1*( ( 1.0_rk - t1 )*w1(i) + t1*w2(i+1) )

          do j = 1, i
             t2 = real( j, rk )*d1
             gvec(j) = 1.5_rk*d1*( 1.0_rk - t1 )*t2*( x(j) + t2 + 1.0_rk )**2

             h(j,j) = h(j,j) + 2.0_rk*t*( 3.0_rk*d1*( 1.0_rk - t1 )*t2*( x(j) + t2 + 1.0_rk ) )
          end do

          gvec(i) = gvec(i) + 1.0_rk

          do j = i+1, global_n
             t2 = real( j, rk )*d1
             gvec(j) = 1.5_rk*d1*t1*( 1.0_rk - t2 )*( x(j) + t2 + 1.0_rk )**2

             h(j,j) = h(j,j) + 2.0_rk*t*( 3.0_rk*d1*t1*( 1.0_rk - t2 )*( x(j) + t2 + 1.0_rk ) )
          end do

          do k = 1, global_n
             do j = 1, k
                h(j,k) = h(j,k) + 2.0_rk*gvec(j)*gvec(k)
             end do
          end do
       end do

    case ( 30 ) ! Broyden tridiagonal
       do j = 1, global_n
          do i = 1, j
             h(i,j) = 0.0_rk
          end do
       end do

       ! For i = 1
       t  = ( 3.0_rk - 2.0_rk*x(1) )*x(1) - 2.0_rk*x(2) + 1.0_rk
       t1 = 3.0_rk - 4.0_rk*x(1)
       h(1,1) =   2.0_rk*( t1**2 - 4.0_rk*t )
       h(1,2) = - 4.0_rk*t1
       h(2,2) =   8.0_rk

       do i = 2, global_n-1
          t  = ( 3.0_rk - 2.0_rk*x(i) )*x(i) - x(i-1) - 2.0_rk*x(i+1) + 1.0_rk
          t1 = 3.0_rk - 4.0_rk*x(i)
          h(i-1,i-1) = h(i-1,i-1) + 2.0_rk
          h(i-1,i  ) = h(i-1,i  ) - 2.0_rk*t1
          h(i,  i  ) = h(i,  i  ) + 2.0_rk*( t1**2 - 4.0_rk*t )
          h(i-1,i+1) = h(i-1,i+1) + 4.0_rk
          h(i  ,i+1) = h(i  ,i+1) - 4.0_rk*t1
          h(i+1,i+1) = h(i+1,i+1) + 8.0_rk
       end do

       ! For i = n
       t  = ( 3.0_rk - 2.0_rk*x(global_n) )*x(global_n) - x(global_n-1) + 1.0_rk
       t1 = 3.0_rk - 4.0_rk*x(global_n)
       h(global_n-1,global_n-1) = h(global_n-1,global_n-1) + 2.0_rk
       h(global_n-1,global_n  ) = h(global_n-1,global_n  ) - 2.0_rk*t1
       h(global_n,  global_n  ) = h(global_n,  global_n  ) + 2.0_rk*( t1**2 - 4.0_rk*t )

    case ( 31 ) ! Broyden banded
       do j = 1, global_n
          do k = 1, j
             h(k,j) = 0.0_rk
          end do
       end do

       do i = 1, global_n
          s1 = 0.0_rk
          do j = max( 1, i-5 ), i-1
             s1 = s1 + x(j)*( 1.0_rk + x(j) )
          end do
          if ( i .ne. global_n ) s1 = s1 + x(i+1)*( 1.0_rk + x(i+1) )
          
          t = x(i)*( 2.0_rk + 5.0_rk*x(i)**2 ) + 1.0_rk - s1
          d1 = 2.0_rk + 15.0_rk*x(i)**2
          do j = max(1,i-5), i-1
             d2 = - 1.0_rk - 2.0_rk*x(j)
             h(j,j) = h(j,j) + 2.0_rk*( d2**2 - 2.0_rk*t )
             do l = j+1, i-1
                h(j,l) = h(j,l) + 2.0_rk*d2*( - 1.0_rk - 2.0_rk*x(l) )
             end do
             h(j,i) = h(j,i) + 2.0_rk*d1*d2
             if ( i .ne. global_n ) then
                h(j,i+1) = h(j,i+1) + 2.0_rk*d2*( - 1.0_rk - 2.0_rk*x(i+1) )
             end if
          end do
          h(i,i) = h(i,i) + 2.0_rk*( 30.0_rk*t*x(i) + d1**2 )
          if ( i .ne. global_n ) then
             d2 = - 1.0_rk - 2.0_rk*x(i+1)
             h(i  ,i+1) = h(i  ,i+1) + 2.0_rk*d1*d2
             h(i+1,i+1) = h(i+1,i+1) + 2.0_rk*( d2**2 - 2.0_rk*t )
          end if
       end do

    case ( 32 ) ! Linear function - full rank
       do j = 1, global_n
          h(j,j) = 2.0_rk
          do k = 1, j-1
             h(k,j) = 0.0_rk
          end do
       end do

    case ( 33 ) ! Linear function - rank 1
       s1 = 0.0_rk
       do i = 1, global_m
          s1 = s1 + i**2
       end do
       s1 = 2.0_rk*s1
       
       do j = 1, global_n
          do i = 1, j
             h(i,j) = i*j*s1
          end do
       end do

    case ( 34 ) ! Linear function - rank 1 with zero columns and rows
       s1 = 0.0_rk
       do i = 2, global_m-1
          s1 = s1 + (i-1)**2
       end do
       s1 = 2.0_rk*s1
       
       do j = 1, global_n
          do i = 1, j
             if ( i .eq. 1 .or. i .eq. global_n .or. &
                  j .eq. 1 .or. j .eq. global_n ) then
                h(i,j) = 0.0_rk
             else
                h(i,j) = i*j*s1
             end if
          end do
       end do

    case ( 35 ) ! Chebyquad
       fvec(1:global_n) = 0.0_rk
       do j = 1, global_n
          t1 = 1.0_rk
          t2 = 2.0_rk*x(j) - 1.0_rk
          t = 2.0_rk*t2
          do i = 1, global_n
             fvec(i) = fvec(i) + t2
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       d1 = 1.0_rk/real( global_n, rk )
       do i = 1, global_n
          fvec(i) = d1*fvec(i)
          if ( modulo( i, 2 ) .eq. 0 ) then
             fvec(i) = fvec(i) + 1.0_rk/( real( i, rk )**2 - 1.0_rk )
          end if
       end do
       d2 = 2.0_rk*d1
       do j = 1, global_n
          h(j,j) = 4.0_rk*d1
          t1 = 1.0_rk
          t2 = 2.0_rk*x(j) - 1.0_rk
          t = 2.0_rk*t2
          s1 = 0.0_rk
          s2 = 2.0_rk
          p1 = 0.0_rk
          p2 = 0.0_rk
          gvec(1) = s2
          do i = 2, global_n
             th = 4.0_rk*t2 + t*s2 - s1
             s1 = s2
             s2 = th
             th = t*t2 - t1
             t1 = t2
             t2 = th
             th = 8.0_rk*s1 + t*p2 - p1
             p1 = p2
             p2 = th
             gvec(i) = s2
             h(j,j) = h(j,j) + fvec(i)*th + d1*s2**2
          end do
          h(j,j) = d2*h(j,j)
          do k = 1, j-1
             h(k,j) = 0.0_rk
             tt1 = 1.0_rk
             tt2 = 2.0_rk*x(k) - 1.0_rk
             tt = 2.0_rk*tt2
             ss1 = 0.0_rk
             ss2 = 2.0_rk
             do i = 1, global_n
                h(k,j) = h(k,j) + ss2*gvec(i)
                tth = 4.0_rk*tt2 + tt*ss2 - ss1
                ss1 = ss2
                ss2 = tth
                tth = tt*tt2 - tt1
                tt1 = tt2
                tt2 = tth
             end do
             h(k,j) = d2*d1*h(k,j)
          end do
       end do

    case default
       flag = - 1

    end select
    
  end subroutine mgh_evalh
  
  ! ------------------------------------------------------------------

  subroutine mgh_evalt( x, t, flag )
    ! SCALAR ARGUMENT
    integer, intent(out) :: flag
    
    ! ARRAY ARGUMENTS
    real(kind=rk), dimension(:),     intent(in)  :: x
    real(kind=rk), dimension(:,:,:), intent(out) :: t

    ! This subroutine computes the third derivative tensor for the
    ! function
    !        \sum_{i=1}^n f_i^2(x)
    ! Returning flag
    !  0    successful
    ! -1    problem is not between 1 and 35
    ! -3    division by zero

    ! LOCAL SCALARS
    integer       :: i, j, k, l
    real(kind=rk) :: arg, d1, d2, d3, d4, d5, d6, h1, h2, h3, r, s1, &
                     s2, t1, t2, t3, t4, td1, td2, td3, theta

    ! LOCAL WORKING ARRAYS
    real(kind=rk), dimension(global_n)                :: w1, w3, w4
    real(kind=rk), dimension(global_n+1)              :: w2
    real(kind=rk), dimension(0:global_n,1:global_n+1) :: prod

    flag = 0

    select case ( problem )

    case ( 1 ) ! Rosenbrock
       t(1,1,1) = 2.4e+3_rk*x(1)
       t(1,1,2) = - 400.0_rk
       t(1,2,2) = 0.0_rk
       t(2,2,2) = 0.0_rk

    case ( 2 ) ! Freudenstein and Roth
       t1 = - 13.0_rk + x(1) + ( ( 5.0_rk - x(2) )*x(2) - 2.0_rk )*x(2)
       t2 = - 29.0_rk + x(1) + ( ( x(2) + 1.0_rk )*x(2) - 14.0_rk )*x(2)
       t(1,1,1) = 0.0_rk
       t(1,1,2) = 0.0_rk
       t(1,2,2) = 24.0_rk
       t(2,2,2) = 6.0_rk*( &
            2.0_rk*t2 + ( ( 3.0_rk*x(2) + 2.0_rk )*x(2) - 14.0_rk )*(   6.0_rk*x(2) + 2.0_rk ) - &
            2.0_rk*t1 + ( ( 10.0_rk - 3.0_rk*x(2) )*x(2) - 2.0_rk )*( - 6.0_rk*x(2) + 10.0_rk ) )

    case ( 3 ) ! Powell badly scaled
       d1 = exp( - x(1) )
       d2 = exp( - x(2) )
       t2 = d1 + d2 - 1.0001_rk
       t(1,1,1) = - 2.0_rk*d1*( t2 + 3.0_rk*d1 )
       t(1,1,2) = 2.0_rk*( 2.0e+8_rk*x(2) - d1*d2 )
       t(1,2,2) = 2.0_rk*( 2.0e+8_rk*x(1) - d1*d2 )
       t(2,2,2) = - 2.0_rk*d2*( t2 + 3.0_rk*d2 )

    case ( 4 ) ! Brown badly scaled
       t(1,1,1) = 0.0_rk
       t(1,1,2) = 4.0_rk*x(2)
       t(1,2,2) = 4.0_rk*x(1)
       t(2,2,2) = 0.0_rk

    case ( 5 ) ! Beale
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do
       
       do i = 1, global_m
          t1 = y5(i) - x(1)*( 1.0_rk - x(2)**i )
          t(1,1,2) = t(1,1,2) + 4.0_rk*i*x(2)**(i-1)*( x(2)**i - 1.0_rk )
          t(1,2,2) = t(1,2,2) + 2.0_rk*(            &
               t1*( i - 1 )*i*x(2)**( i - 2 ) +     &
               2.0_rk*i**2*x(1)*x(2)**( 2*i - 2 ) + &
               ( x(2)**i - 1.0_rk )*( i - 1 )*i*x(1)*x(2)**( i - 2 ) )
          t(2,2,2) = t(2,2,2) + 2.0_rk*(                        &
               t1*( i - 2 )*( i - 1 )*i*x(1)*x(2)**( i - 3 ) +  &
               3.0_rk*( i - 1 )*i**2*x(1)**2*x(2)**( 2*i - 3 ) )
       end do

    case ( 6 ) ! Jennrich and Sampson
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = exp( i*x(1) )
          d2 = exp( i*x(2) )
          t1 = 2.0_rk + 2.0_rk*i - ( d1 + d2 )
          t(1,1,1) = t(1,1,1) + 2.0_rk*i**3*d1*( 3.0_rk*d1 - t1 )
          t(1,1,2) = t(1,1,2) + 2.0_rk*i**3*d1*d2
          t(1,2,2) = t(1,2,2) + 2.0_rk*i**3*d1*d2
          t(2,2,2) = t(2,2,2) + 2.0_rk*i**3*d2*( 3.0_rk*d2 - t1 )
       end do

    case ( 7 ) ! Helical valley
       if ( x(1) .ne. 0.0_rk ) then
          theta = 1.0_rk/( 2.0_rk*PI )*atan( x(2)/x(1) )
          if ( x(1) .lt. 0.0_rk ) theta = theta + 0.5_rk
       else
          theta = huge( rk )
          flag = - 3
       end if
       t1 = 10.0_rk*( x(3) - 10.0_rk*theta )
       t2 = 10.0_rk*( sqrt( x(1)**2 + x(2)**2 ) - 1 )

       s1 = ( x(1)**2 + x(2)**2 )
       if ( s1 .ne. 0.0_rk ) then
          d1 = 50.0_rk/( PI*s1 )
          d2 = 10.0_rk/sqrt( s1 )
          d3 = 1.0_rk/( PI*s1**2 )
          d4 = 10.0_rk/s1**1.5
          d5 = 100.0_rk/( PI*s1**3 )
          d6 = 10.0_rk/s1**2.5
       else
          d1 = huge( rk )
          d2 = huge( rk )
          d3 = huge( rk )
          d4 = huge( rk )
          d5 = huge( rk )
          d6 = huge( rk )
       end if

       t(1,1,1) = 2.0_rk*x(2)*(                                              &
            t1*d5*( 3.0_rk*x(1)**2 - x(2)**2 ) - 3.0e+2_rk*d1*x(1)*x(2)*d3 + &
            t2*( - 3.0_rk*d6*x(1)*x(2) ) + 3.0_rk*d2*x(1)*x(2)*d4 )
       t(1,1,2) = 2.0_rk*(                                                                      &
            t1*d5*( 3.0_rk*x(1)*x(2)**2 - x(1)**3 ) + t2*d6*( 2.0_rk*x(1)**2*x(2) - x(2)**3 ) + &
            x(2)*( 2.0_rk*x(1)**2 - x(2)**2 )*( 1.0e+2_rk*d1*d3 - d2*d4 ) )
       t(1,2,2) = 2.0_rk*(                                                                      &
            t2*d6*( 2.0_rk*x(1)*x(2)**2 - x(1)**3 ) - t1*d5*( 3.0_rk*x(1)**2*x(2) - x(2)**3 ) + &
            x(1)*( 2.0_rk*x(2)**2 - x(1)**2 )*( 1.0e+2_rk*d1*d3 - d2*d4 ) )
       t(2,2,2) = 2.0_rk*(                                                          &
            - t1*d5*( 3.0_rk*x(1)*x(2)**2 - x(1)**3 ) - t2*3.0_rk*d6*x(1)**2*x(2) + &
            3.0_rk*x(1)**2*x(2)*( d2*d4 - 1.0e+2_rk*d1*d3 ) )
       t(1,1,3) = - 2.0e+3_rk*d3*x(1)*x(2)
       t(1,2,3) = 1.0e+3_rk*d3*( x(1)**2 - x(2)**2 )
       t(2,2,3) = 2.0e+3_rk*d3*x(1)*x(2)
       t(1,3,3) = 0.0_rk
       t(2,3,3) = 0.0_rk
       t(3,3,3) = 0.0_rk

    case ( 8 ) ! Bard
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = real( i, rk )
          d2 = 16.0_rk - real( i, rk )
          d3 = min( d1, d2 )

          s1 = d2*x(2) + d3*x(3)
          if ( s1 .ne. 0.0_rk ) then
             t1 = y8(i) - ( x(1) + d1/s1 )
             t2 = d1/s1**2
             t3 = - 2.0_rk*d1/s1**3
             t4 = 6.0_rk*d1/s1**4
          else
             t1 = huge( rk )
             t2 = huge( rk )
             t3 = huge( rk )
             t4 = huge( rk )
             flag = - 3
          end if

          t(1,2,2) = t(1,2,2) - 2.0_rk*t3*d2**2
          t(2,2,2) = t(2,2,2) + 2.0_rk*d2**3*( t1*t4 + 3.0_rk*t2*t3 )
          t(1,2,3) = t(1,2,3) - 2.0_rk*t3*d2*d3
          t(2,2,3) = t(2,2,3) + 2.0_rk*d2**2*d3*( t1*t4 + 3.0_rk*t2*t3 )
          t(1,3,3) = t(1,3,3) - 2.0_rk*t3*d3**2
          t(2,3,3) = t(2,3,3) + 2.0_rk*d2*d3**2*( t1*t4 + 3.0_rk*t2*t3 )
          t(3,3,3) = t(3,3,3) + 2.0_rk*d3**3*( t1*t4 + 3.0_rk*t2*t3 )
       end do

    case ( 9 ) ! Gaussian
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = ( 8.0_rk - real( i, rk ) )/2.0_rk
          d2 = d1 - x(3)
          arg = 0.5_rk*x(2)*d2**2
          r = exp( - arg )
          t1 = x(1)*r - y9(i)
          t(1,1,2) = t(1,1,2) - 4.0_rk*( r*d2 )**2/2.0_rk
          t(1,2,2) = t(1,2,2) + 2.0_rk*     r*d2**4/4.0_rk*( t1 + 3.0_rk*x(1)*r )
          t(2,2,2) = t(2,2,2) - 2.0_rk*x(1)*r*d2**6/8.0_rk*( t1 + 3.0_rk*x(1)*r )
          t(1,1,3) = t(1,1,3) + 4.0_rk*r**2*x(2)*d2
          t(1,2,3) = t(1,2,3) + 2.0_rk*( &
               ( t1*d2*r + r**2*x(1)*d2 )*( 1.0_rk - x(2)*d2**2/2.0_rk ) - x(1)*x(2)*d2**3*r**2 )
          t(2,2,3) = t(2,2,3) + 2.0_rk*(                          &
               t1*x(1)*d2**3*r*( x(2)*d2**2 / 4.0_rk - 1.0_rk ) + &
               r**2*x(1)**2*x(2)*d2**5/4.0_rk -                   &
               r**2*x(1)**2*d2**3*( 1.0_rk - x(2)*d2**2 / 2.0_rk ) )
          t(1,3,3) = t(1,3,3) + 2.0_rk*(                                &
               ( t1*x(2)*r + r**2*x(1)*x(2) )*( x(2)*d2**2 - 1.0_rk ) + &
               2.0_rk*x(1)*( x(2)*d2*r )**2 )
          t(2,3,3) = t(2,3,3) + 2.0_rk*(                                            &
               t1*x(1)*r*( x(2)*d2**2*( 2.5_rk - x(2)*d2**2 / 2.0_rk ) - 1.0_rk ) + &
               x(1)**2*x(2)*d2**2*r**2*( 2.5_rk - 1.5_rk*x(2)*d2**2 ) )
          t(3,3,3) = t(3,3,3) + 2.0_rk*(                      &
               t1*x(1)*x(2)**2*d2*r*( x(2)*d2**2 - 3.0_rk ) + &
               3.0_rk*( x(1)*x(2)*r )**2*d2*( x(2)*d2**2 - 1.0_rk ) )
       end do

    case ( 10 ) ! Meyer
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = 4.5e+01_rk + 5.0_rk*real( i, rk )
          d2 = d1 + x(3)
          if ( d2 .ne. 0.0_rk ) then
             d3 = exp( x(2)/d2 )
             t1 = x(1)*d3 - y10(i)
             t(1,1,2) = t(1,1,2) + 4.0_rk*d3**2/d2
             t(1,2,2) = t(1,2,2) + 2.0_rk*d3/d2**2*( t1 + 3.0_rk*x(1)*d3 )
             t(2,2,2) = t(2,2,2) + 2.0_rk*x(1)*d3/d2**3*( t1 + 3.0_rk*x(1)*d3 )
             t(1,1,3) = t(1,1,3) - 4.0_rk*x(2)*( d3/d2 )**2
             t(1,2,3) = t(1,2,3) - 2.0_rk*(                                  &
                  ( t1*d3/d2**2 + x(1)*( d3/d2 )**2 )*( x(2)/d2 + 1.0_rk ) + &
                  2.0_rk*x(1)*x(2)*d3**2/d2**3 )
             t(2,2,3) = t(2,2,3) - 2.0_rk*(                         &
                  t1*d3*x(1)/d2**3*( x(2)/d2 + 2.0_rk ) +           &
                  2.0_rk*x(1)**2*d3**2/d2**3*( x(2)/d2 + 1.0_rk ) + &
                  x(1)**2*x(2)*d3**2/d2**4 )
             t(1,3,3) = t(1,3,3) + 2.0_rk*(                    &
                  t1*d3*x(2)/d2**3*( x(2)/d2 + 2.0_rk ) +      &
                  x(1)*x(2)*d3**2/d2**3*( x(2)/d2 + 2.0_rk ) + &
                  2.0_rk*x(1)*x(2)**2*d3**2/d2**4 )
             t(2,3,3) = t(2,3,3) + 2.0_rk*(                               &
                  t1*d3*x(1)/d2**3*( ( x(2)/d2 + 2.0_rk )**2 - 2.0_rk ) + &
                  x(1)**2*x(2)*d3**2/d2**4*( x(2)/d2 + 2.0_rk ) +         &
                  2.0_rk*x(1)**2*x(2)*d3**2/d2**4*( x(2)/d2 + 1.0_rk ) )
             t(3,3,3) = t(3,3,3) - 2.0_rk*(                                    &
                  t1*d3*x(1)*x(2)/d2**4*( ( x(2)/d2 + 3.0_rk )**2 - 3.0_rk ) + &
                  3.0_rk*( x(1)*x(2)*d3 )**2/d2**5*( x(2)/d2 + 2.0_rk ) )
          else
             t(1:global_n,1:global_n,1:global_n) = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 11 ) ! Gulf
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = real( i, rk ) / 100.0_rk
          d2 = 25.0_rk + ( - 50.0_rk*log( d1 ) )**( 2.0_rk/3.0_rk )
          d3 = log( abs( d2 - x(2) ) )
          d4 = abs( d2 - x(2) )**x(3)
          if ( x(1) .ne. 0.0_rk .and. d2 .ne. x(2) ) then
             t2 = exp( - d4/x(1) )
             t1 = t2 - d1
             t(1,1,1) = t(1,1,1) + 2.0_rk*(                                           &
                  t1*t2*( d4/( x(1)**4 ) )*( ( ( d4/x(1) - 3.0_rk )**2 - 3.0_rk ) ) + &
                  3.0_rk*t2**2*d4/x(1)**2*( d4/( x(1)**3 ) )*( d4/x(1) - 2.0_rk ) )
             t(1,1,2) = t(1,1,2) + 2.0_rk*(                                                      &
                  t1*t2*x(3)*d4/( x(1)**3*( d2 - x(2) ) )*( ( d4/x(1) - 2.0_rk )**2 - 2.0_rk ) + &
                  t2*x(3)*d4/( x(1)*( d2 - x(2) ) )*t2*( d4/( x(1)**3 ) )*( d4/x(1) - 2.0_rk ) + &
                  2.0_rk*t2*d4/( x(1)**2 )*t2*x(3)*d4/( x(1)**2*( d2 - x(2) ) )*( d4/x(1) - 1.0_rk ) )
             t(1,2,2) = t(1,2,2) + 2.0_rk*(                                                &
                  t1*t2*x(3)*d4/( x(1)**2*( d2 - x(2) )**2 )*                              &
                  ( x(3)*( d4/x(1) - 1.0_rk )**2 + d4/x(1)*( 1.0_rk - x(3) ) - 1.0_rk ) +  &
                  t2*d4/( x(1)**2 )*t2*(x(3)*d4/( x(1)*( d2 - x(2) )**2 )*( x(3)*d4/x(1) - &
                  x(3) + 1.0_rk ) ) + 2.0_rk*t2*x(3)*d4/( x(1)*( d2 - x(2) ) )*t2*x(3)*d4/ &
                  ( x(1)**2*( d2 - x(2) ) )*( d4/x(1) - 1.0_rk ) )
             t(2,2,2) = t(2,2,2) + 2.0_rk*(                                &
                  t1*t2*x(3)*d4/( x(1)*( d2 - x(2) )**3 )*(                &
                  ( x(3)*d4/x(1) - x(3) + 1.0_rk )**2 -                    &
                  x(3)*( d4/x(1)*( x(3) - 1.0_rk ) + 1.0_rk ) + 1.0_rk ) + &
                  3.0_rk*t2*x(3)*d4/( x(1)*( d2 - x(2) ) )*t2*(x(3)*d4/    &
                  ( x(1)*( d2 - x(2) )**2 )*( x(3)*d4/x(1) - x(3) + 1.0_rk ) ) )
             t(1,1,3) = t(1,1,3) + 2.0_rk*( &
                  t1*t2*d4/( x(1)**3 )*d3*( - ( d4/x(1) - 2.0_rk )**2 + 2.0_rk ) +     &
                  2.0_rk*t2*d4/( x(1)**2 )*t2*d3*d4/( x(1)**2 )*( 1.0_rk - d4/x(1) ) - &
                  t2*d3*d4/x(1)*t2*( d4/( x(1)**3 ) )*( d4/x(1) - 2.0_rk ) )
             t(1,2,3) = t(1,2,3) + 2.0_rk*(                                                                    &
                  t1*t2*d4/( x(1)**2*( d2 - x(2) ) )*( - x(3)*d3*(                                             &
                  ( d4/x(1) - 1.0_rk )**2 - d4/x(1) ) + d4/x(1) - 1.0_rk ) +                                   &
                  t2*d4/( x(1)**2 )*t2*d4/( x(1)*( d2 - x(2) ) )*( x(3)*d3*( - d4/x(1) + 1.0_rk ) + 1.0_rk ) + &
                  t2*x(3)*d4/( x(1)*( d2 - x(2) ) )*t2*d3*d4/( x(1)**2 )*( 1.0_rk - d4/x(1) ) -                &
                  t2*d3*d4/x(1)*t2*x(3)*d4/( x(1)**2*( d2 - x(2) ) )*( d4/x(1) - 1.0_rk ) )
             t(2,2,3) = t(2,2,3) + 2.0_rk*(                                                                &
                  t1*t2*d4/( x(1)*( d2 - x(2) )**2 )*( x(3)*d4/x(1)*( x(3)*d3 + 1.0_rk ) -                 &
                  x(3) + ( x(3)*d4/x(1) - x(3) + 1.0_rk )*( - x(3)*d3*d4/x(1) + x(3)*d3 + 1.0_rk ) ) -     &
                  t2*d3*d4/x(1)*t2*(x(3)*d4/( x(1)*( d2 - x(2) )**2 )*( x(3)*d4/x(1) - x(3) + 1.0_rk ) ) + &
                  2.0_rk*t2*x(3)*d4/( x(1)*( d2 - x(2) ) )*t2*d4/( x(1)*( d2 - x(2) ) )*(                  &
                  x(3)*d3*( - d4/x(1) + 1.0_rk ) + 1.0_rk ) )
             t(1,3,3) = t(1,3,3) + 2.0_rk*(                                          &
                  t1*t2*d3**2*d4/( x(1)**2 )*( ( d4/x(1) - 1.0_rk )**2 - d4/x(1) ) + &
                  t2*d4/( x(1)**2 )*t2*d4/x(1)*d3**2*( d4/x(1) - 1.0_rk ) -          &
                  2.0_rk*( t2*d3*d4 )**2/x(1)**3*( 1.0_rk - d4/x(1) ) )
             t(2,3,3) = t(2,3,3) + 2.0_rk*(                                                        &
                  t1*t2*d3*d4/( x(1)*( d2 - x(2) ) )*( ( - d4/x(1) + 1.0_rk )*( -                  &
                  x(3)*d3*d4/x(1) + x(3)*d3 + 1.0_rk ) - d4/x(1)*( x(3)*d3 + 1.0_rk ) + 1.0_rk ) + &
                  t2*x(3)*d4/( x(1)*( d2 - x(2) ) )*t2*d4/x(1)*d3**2*( d4/x(1) - 1.0_rk ) -        &
                  2.0_rk*t2*d3*d4/x(1)*t2*d4/( x(1)*( d2 - x(2) ) )*( x(3)*d3*( - d4/x(1) + 1.0_rk ) + 1.0_rk ) )
             t(3,3,3) = t(3,3,3) + 2.0_rk*(                                   &
                  t1*t2*d3**3*d4/x(1)*( d4/x(1) - ( d4/x(1) - 1.0_rk )**2 ) - &
                  3.0_rk*t2*d3*d4/x(1)*t2*d4/x(1)*d3**2*( d4/x(1) - 1.0_rk ) )
          else
             t(1:global_n,1:global_n,1:global_n) = huge( rk )
             flag = - 3
             return
          end if
       end do

    case ( 12 ) ! Box three-dimensional
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = 0.1_rk*real( i, rk )
          t1 = exp( - d1*x(1) )
          t2 = exp( - d1*x(2) )
          t3 = exp( - 10.0_rk*d1 )
          t4 = exp( - d1 )
          s1 = t1 - t2 - x(3)*( t4 - t3 )
          t(1,1,1) = t(1,1,1) - 2.0_rk*d1**3*t1*( s1 + 3.0_rk*t1 )
          t(1,1,2) = t(1,1,2) + 2.0_rk*d1**3*t1*t2
          t(1,2,2) = t(1,2,2) + 2.0_rk*d1**3*t1*t2
          t(2,2,2) = t(2,2,2) + 2.0_rk*d1**3*t2*( s1 - 3.0_rk*t2 )
          t(1,1,3) = t(1,1,3) + 2.0_rk*( t3 - t4 )*d1**2*t1
          t(2,2,3) = t(2,2,3) - 2.0_rk*( t3 - t4 )*d1**2*t2
       end do

    case ( 13 ) ! Powell singular
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do
       t1 = 2.0_rk*( x(2) - 2.0_rk*x(3) )
       t2 = 2.0_rk*sqrt( 10.0_rk )*( x(1) - x(4) )
       t3 = 2.0_rk*sqrt( 10.0_rk )
       d1 =  6.0_rk*t2*t3
       d2 = 12.0_rk*t1
       t(1,1,1) =   d1
       t(2,2,2) =   d2
       t(2,2,3) = - d2*2.0_rk
       t(2,3,3) =   d2*4.0_rk
       t(3,3,3) = - d2*8.0_rk
       t(1,1,4) = - d1
       t(1,4,4) =   d1
       t(4,4,4) = - d1

    case ( 14 ) ! Wood
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do
       t(1,1,1) =   2.40e+3_rk*x(1)
       t(1,1,2) = - 4.00e+2_rk
       t(3,3,3) =   2.16e+3_rk*x(3)
       t(3,3,4) = - 3.6e+2_rk

    case ( 15 ) ! Kowalik and Osborne
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do
       do i = 1, global_m
          s1 = u15(i) + x(2)
          t1 = u15(i)**2 + u15(i)*x(3) + x(4)
          if ( t1 .ne. 0.0_rk ) then
             d1 = y15(i) - x(1)*u15(i)*s1/t1
             t(1,1,2) = t(1,1,2) + 4.0_rk*s1*( u15(i)/t1 )**2
             t(1,2,2) = t(1,2,2) + 4.0_rk*x(1)*( u15(i)/t1 )**2
             t(1,1,3) = t(1,1,3) - 4.0_rk*s1**2*( u15(i)/t1 )**3
             t(1,2,3) = t(1,2,3) + 2.0_rk*( u15(i)/t1 )**2*( d1 - 3.0_rk*x(1)*s1*u15(i)/t1 )
             t(2,2,3) = t(2,2,3) - 4.0_rk*x(1)**2*( u15(i)/t1 )**3
             t(1,3,3) = t(1,3,3) + 4.0_rk*s1*( u15(i)/t1 )**3*( 2.0_rk*x(1)*s1*u15(i)/t1 - d1 )
             t(2,3,3) = t(2,3,3) + 4.0_rk*x(1)*( u15(i)/t1 )**3*( 2.0_rk*x(1)*s1*u15(i)/t1 - d1 )
             t(3,3,3) = t(3,3,3) + 12.0_rk*x(1)*s1*( u15(i)/t1 )**4*( d1 - x(1)*s1*u15(i)/t1 )
             t(1,1,4) = t(1,1,4) - 4.0_rk*s1**2*u15(i)**2/t1**3
             t(1,2,4) = t(1,2,4) + 2.0_rk*u15(i)/t1**2*( d1 - 3.0_rk*x(1)*s1*u15(i)/t1 )
             t(2,2,4) = t(2,2,4) - 4.0_rk*x(1)**2*u15(i)**2/t1**3
             t(1,3,4) = t(1,3,4) + 4.0_rk*s1*u15(i)**2/t1**3*( 2.0_rk*x(1)*s1*u15(i)/t1 - d1 )
             t(2,3,4) = t(2,3,4) + 4.0_rk*x(1)*u15(i)**2/t1**3*( 2.0_rk*x(1)*s1*u15(i)/t1 - d1 )
             t(3,3,4) = t(3,3,4) + 12.0_rk*x(1)*s1*u15(i)**3/t1**4*( d1 - x(1)*s1*u15(i)/t1 )
             t(1,4,4) = t(1,4,4) + 4.0_rk*s1*u15(i)/t1**3*( 2.0_rk*x(1)*s1*u15(i)/t1 - d1 )
             t(2,4,4) = t(2,4,4) + 4.0_rk*x(1)*u15(i)/t1**3*( - d1 + 2.0_rk*x(1)*s1*u15(i)/t1 )
             t(3,4,4) = t(3,4,4) + 12.0_rk*x(1)*s1*( u15(i)/t1**2 )**2*( d1 - x(1)*s1*u15(i)/t1 )
             t(4,4,4) = t(4,4,4) + 12.0_rk*x(1)*s1*u15(i)/t1**4*( d1 - x(1)*s1*u15(i)/t1 )
          end if
       end do

    case ( 16 ) ! Brown and Dennis
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = real( i, rk )/5.0_rk
          d2 = exp( d1 )
          d3 = sin( d1 )
          d4 = cos( d1 )
          t(1,1,1) = t(1,1,1) + 24.0_rk*( x(1) + d1*x(2) - d2 )
          t(1,1,2) = t(1,1,2) + 24.0_rk*( x(1) + d1*x(2) - d2 )*d1
          t(1,2,2) = t(1,2,2) + 24.0_rk*( x(1) + d1*x(2) - d2 )*d1**2
          t(2,2,2) = t(2,2,2) + 24.0_rk*( x(1) + d1*x(2) - d2 )*d1**3
          t(1,1,3) = t(1,1,3) + 8.0_rk*( x(3) + d3*x(4) - d4 )
          t(1,2,3) = t(1,2,3) + 8.0_rk*( x(3) + d3*x(4) - d4 )*d1
          t(2,2,3) = t(2,2,3) + 8.0_rk*( x(3) + d3*x(4) - d4 )*d1**2
          t(1,3,3) = t(1,3,3) + 8.0_rk*( x(1) + d1*x(2) - d2 )
          t(2,3,3) = t(2,3,3) + 8.0_rk*( x(1) + d1*x(2) - d2 )*d1
          t(3,3,3) = t(3,3,3) + 24.0_rk*( x(3) + d3*x(4) - d4 )
          t(1,1,4) = t(1,1,4) + 8.0_rk*( x(3) + d3*x(4) - d4 )*d3
          t(1,2,4) = t(1,2,4) + 8.0_rk*( x(3) + d3*x(4) - d4 )*d1*d3
          t(2,2,4) = t(2,2,4) + 8.0_rk*( x(3) + d3*x(4) - d4 )*d1**2*d3
          t(1,3,4) = t(1,3,4) + 8.0_rk*( x(1) + d1*x(2) - d2 )*d3
          t(2,3,4) = t(2,3,4) + 8.0_rk*( x(1) + d1*x(2) - d2 )*d1*d3
          t(3,3,4) = t(3,3,4) + 24.0_rk*( x(3) + x(4)*d3 - d4 )*d3
          t(1,4,4) = t(1,4,4) + 8.0_rk*( x(1) + d1*x(2) - d2 )*d3**2
          t(2,4,4) = t(2,4,4) + 8.0_rk*( x(1) + d1*x(2) - d2 )*d1*d3**2
          t(3,4,4) = t(3,4,4) + 24.0_rk*( x(3) + x(4)*d3 - d4 )*d3**2
          t(4,4,4) = t(4,4,4) + 24.0_rk*( x(3) + x(4)*d3 - d4 )*d3**3
       end do

    case ( 17 ) ! Osborne 1
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = 10.0_rk*real( i-1, rk )
          d2 = exp( - d1*x(4) )
          d3 = exp( - d1*x(5) )
          t1 = y17(i) - ( x(1) + x(2)*d2 + x(3)*d3 )
          t(1,2,4) = t(1,2,4) - 2.0_rk*d1*d2
          t(2,2,4) = t(2,2,4) - 4.0_rk*d1*d2**2
          t(2,3,4) = t(2,3,4) - 2.0_rk*d1*d2*d3
          t(1,4,4) = t(1,4,4) + 2.0_rk*d1**2*d2*x(2)
          t(2,4,4) = t(2,4,4) + 2.0_rk*d1**2*d2*( 3.0_rk*d2*x(2) - t1 )
          t(3,4,4) = t(3,4,4) + 2.0_rk*d1**2*d2*d3*x(2)
          t(4,4,4) = t(4,4,4) + 2.0_rk*d1**3*d2*x(2)*( t1 - 3.0_rk*d2*x(2) )
          t(1,3,5) = t(1,3,5) - 2.0_rk*d1*d3
          t(2,3,5) = t(2,3,5) - 2.0_rk*d1*d2*d3 
          t(3,3,5) = t(3,3,5) - 4.0_rk*d1*d3**2
          t(2,4,5) = t(2,4,5) + 2.0_rk*d1**2*d2*d3*x(3)
          t(3,4,5) = t(3,4,5) + 2.0_rk*d1**2*d2*d3*x(2)
          t(4,4,5) = t(4,4,5) - 2.0_rk*d1**3*d2*d3*x(2)*x(3)
          t(1,5,5) = t(1,5,5) + 2.0_rk*d1**2*d3*x(3)
          t(2,5,5) = t(2,5,5) + 2.0_rk*d1**2*d2*d3*x(3)
          t(3,5,5) = t(3,5,5) + 2.0_rk*d1**2*d3*( 3.0_rk*d3*x(3) - t1 )
          t(4,5,5) = t(4,5,5) - 2.0_rk*d1**3*d2*d3*x(2)*x(3)
          t(5,5,5) = t(5,5,5) + 2.0_rk*d1**3*d3*x(3)*( t1 - 3.0_rk*d3*x(3) )
       end do

    case ( 18 ) ! Biggs EXP6
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          d1 = 0.1_rk*real( i, rk )
          d2 = exp( - d1*x(1) )
          d3 = exp( - d1*x(2) )
          d4 = exp( - d1*x(5) )
          t1 = exp( - d1 ) - 5.0_rk*exp( - 10.0_rk*d1 ) + 3.0_rk*exp( - 4.0_rk*d1 )
          t2 = x(3)*d2 - x(4)*d3 + x(6)*d4 - t1
          t(1,1,1) = t(1,1,1) - 2.0_rk*d1**3*d2*x(3)*( t2 + 3.0_rk*d2*x(3) )
          t(1,1,2) = t(1,1,2) + 2.0_rk*d1**3*d2*d3*x(3)*x(4)
          t(1,2,2) = t(1,2,2) + 2.0_rk*d1**3*d2*d3*x(3)*x(4)
          t(2,2,2) = t(2,2,2) + 2.0_rk*d1**3*d3*x(4)*( t2 - 3.0_rk*d3*x(4) )
          t(1,1,3) = t(1,1,3) + 2.0_rk*d1**2*d2*( t2 + 3.0_rk*d2*x(3) )
          t(1,2,3) = t(1,2,3) - 2.0_rk*d1**2*d2*d3*x(4)
          t(2,2,3) = t(2,2,3) - 2.0_rk*d1**2*d2*d3*x(4)
          t(1,3,3) = t(1,3,3) - 4.0_rk*d1*d2**2
          t(1,1,4) = t(1,1,4) - 2.0_rk*d1**2*d2*d3*x(3)
          t(1,2,4) = t(1,2,4) - 2.0_rk*d1**2*d2*d3*x(3)
          t(2,2,4) = t(2,2,4) + 2.0_rk*d1**2*d3*( 3.0_rk*d3*x(4) - t2 )
          t(1,3,4) = t(1,3,4) + 2.0_rk*d1*d2*d3
          t(2,3,4) = t(2,3,4) + 2.0_rk*d1*d2*d3
          t(2,4,4) = t(2,4,4) - 4.0_rk*d1*d3**2
          t(1,1,5) = t(1,1,5) - 2.0_rk*d1**3*d2*d4*x(3)*x(6)
          t(2,2,5) = t(2,2,5) + 2.0_rk*d1**3*d3*d4*x(4)*x(6)
          t(1,3,5) = t(1,3,5) + 2.0_rk*d1**2*d2*d4*x(6)
          t(2,4,5) = t(2,4,5) - 2.0_rk*d1**2*d3*d4*x(6)
          t(1,5,5) = t(1,5,5) - 2.0_rk*d1**3*d2*d4*x(3)*x(6)
          t(2,5,5) = t(2,5,5) + 2.0_rk*d1**3*d3*d4*x(4)*x(6)
          t(3,5,5) = t(3,5,5) + 2.0_rk*d1**2*d2*d4*x(6)
          t(4,5,5) = t(4,5,5) - 2.0_rk*d1**2*d3*d4*x(6)
          t(5,5,5) = t(5,5,5) - 2.0_rk*d1**3*d4*x(6)*( t2 + 3.0_rk*d4*x(6) )
          t(1,1,6) = t(1,1,6) + 2.0_rk*d1**2*d2*d4*x(3)
          t(2,2,6) = t(2,2,6) - 2.0_rk*d1**2*d3*d4*x(4)
          t(1,3,6) = t(1,3,6) - 2.0_rk*d1*d2*d4
          t(2,4,6) = t(2,4,6) + 2.0_rk*d1*d3*d4
          t(1,5,6) = t(1,5,6) + 2.0_rk*d1**2*d2*d4*x(3)
          t(2,5,6) = t(2,5,6) - 2.0_rk*d1**2*d3*d4*x(4)
          t(3,5,6) = t(3,5,6) - 2.0_rk*d1*d2*d4
          t(4,5,6) = t(4,5,6) + 2.0_rk*d1*d3*d4
          t(5,5,6) = t(5,5,6) + 2.0_rk*d1**2*d4*( t2 + 3.0_rk*d4*x(6) )
          t(5,6,6) = t(5,6,6) - 4.0_rk*d1*d4**2
       end do

    case ( 19 ) ! Osborne 2
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          t1 = real( i-1, rk ) / 10.0_rk
          d1 = exp( - t1*x(5) )
          d2 = exp( - ( t1 - x(9) )**2*x(6) )
          d3 = exp( - ( t1 - x(10) )**2*x(7) )
          d4 = exp( - ( t1 - x(11) )**2*x(8) )
          s1 = y19(i) - ( x(1)*d1 + x(2)*d2 + x(3)*d3 + x(4)*d4 )
          t( 1, 1, 5) = t( 1, 1, 5) - 4.0_rk*t1*d1**2
          t( 1, 2, 5) = t( 1, 2, 5) - 2.0_rk*t1*d1*d2
          t( 1, 3, 5) = t( 1, 3, 5) - 2.0_rk*t1*d1*d3
          t( 1, 4, 5) = t( 1, 4, 5) - 2.0_rk*t1*d1*d4
          t( 1, 5, 5) = t( 1, 5, 5) - 2.0_rk*t1**2*d1*( s1 - 3.0_rk*d1*x(1) )
          t( 2, 5, 5) = t( 2, 5, 5) + 2.0_rk*t1**2*d1*d2*x(1)
          t( 3, 5, 5) = t( 3, 5, 5) + 2.0_rk*t1**2*d1*d3*x(1)
          t( 4, 5, 5) = t( 4, 5, 5) + 2.0_rk*t1**2*d1*d4*x(1)
          t( 5, 5, 5) = t( 5, 5, 5) + 2.0_rk*t1**3*d1*x(1)*( s1 - 3.0_rk*d1*x(1) )
          t( 1, 2, 6) = t( 1, 2, 6) - 2.0_rk*d1*d2*( t1 - x(9) )**2
          t( 2, 2, 6) = t( 2, 2, 6) - 4.0_rk*d2**2*( t1 - x(9) )**2
          t( 2, 3, 6) = t( 2, 3, 6) - 2.0_rk*d2*d3*( t1 - x(9) )**2
          t( 2, 4, 6) = t( 2, 4, 6) - 2.0_rk*d2*d4*( t1 - x(9) )**2
          t( 1, 5, 6) = t( 1, 5, 6) + 2.0_rk*t1*d1*d2*x(2)*( t1 - x(9) )**2
          t( 2, 5, 6) = t( 2, 5, 6) + 2.0_rk*t1*d1*d2*x(1)*( t1 - x(9) )**2
          t( 5, 5, 6) = t( 5, 5, 6) - 2.0_rk*t1**2*d1*d2*x(1)*x(2)*( t1 - x(9) )**2
          t( 1, 6, 6) = t( 1, 6, 6) + 2.0_rk*d1*d2*x(2)*( t1 - x(9) )**4
          t( 2, 6, 6) = t( 2, 6, 6) - 2.0_rk*d2*( t1 - x(9) )**4*( s1 - 3.0_rk*d2*x(2) )
          t( 3, 6, 6) = t( 3, 6, 6) + 2.0_rk*d2*d3*x(2)*( t1 - x(9) )**4
          t( 4, 6, 6) = t( 4, 6, 6) + 2.0_rk*d2*d4*x(2)*( t1 - x(9) )**4
          t( 5, 6, 6) = t( 5, 6, 6) - 2.0_rk*t1*d1*d2*x(1)*x(2)*( t1 - x(9) )**4
          t( 6, 6, 6) = t( 6, 6, 6) + 2.0_rk*d2*x(2)*( t1 - x(9) )**6*( s1 - 3.0_rk*d2*x(2) )
          t( 1, 3, 7) = t( 1, 3, 7) - 2.0_rk*d1*d3*( t1 - x(10) )**2
          t( 2, 3, 7) = t( 2, 3, 7) - 2.0_rk*d2*d3*( t1 - x(10) )**2
          t( 3, 3, 7) = t( 3, 3, 7) - 4.0_rk*d3**2*( t1 - x(10) )**2
          t( 3, 4, 7) = t( 3, 4, 7) - 2.0_rk*d3*d4*( t1 - x(10) )**2
          t( 1, 5, 7) = t( 1, 5, 7) + 2.0_rk*t1*d1*d3*x(3)*( t1 - x(10) )**2
          t( 3, 5, 7) = t( 3, 5, 7) + 2.0_rk*t1*d1*d3*x(1)*( t1 - x(10) )**2 
          t( 5, 5, 7) = t( 5, 5, 7) - 2.0_rk*t1**2*d1*d3*x(1)*x(3)*( t1 - x(10) )**2
          t( 2, 6, 7) = t( 2, 6, 7) + 2.0_rk*d2*d3*x(3)*( t1 - x(9) )**2*( t1 - x(10) )**2
          t( 3, 6, 7) = t( 3, 6, 7) + 2.0_rk*d2*d3*x(2)*( t1 - x(9) )**2*( t1 - x(10) )**2
          t( 6, 6, 7) = t( 6, 6, 7) - 2.0_rk*d2*d3*x(2)*x(3)*( t1 - x(9) )**4*( t1 - x(10) )**2
          t( 1, 7, 7) = t( 1, 7, 7) + 2.0_rk*d1*d3*x(3)*( t1 - x(10) )**4
          t( 2, 7, 7) = t( 2, 7, 7) + 2.0_rk*d2*d3*x(3)*( t1 - x(10) )**4
          t( 3, 7, 7) = t( 3, 7, 7) + 2.0_rk*d3*( t1 - x(10) )**4*( 3.0_rk*d3*x(3) - s1 )
          t( 4, 7, 7) = t( 4, 7, 7) + 2.0_rk*d3*d4*x(3)*( t1 - x(10) )**4
          t( 5, 7, 7) = t( 5, 7, 7) - 2.0_rk*t1*d1*d3*x(1)*x(3)*( t1 - x(10) )**4
          t( 6, 7, 7) = t( 6, 7, 7) - 2.0_rk*d2*d3*x(2)*x(3)*( t1 - x(9) )**2*( t1 - x(10) )**4
          t( 7, 7, 7) = t( 7, 7, 7) + 2.0_rk*d3*x(3)*( t1 - x(10) )**6*( s1 - 3.0_rk*d3*x(3) )
          t( 1, 4, 8) = t( 1, 4, 8) - 2.0_rk*d1*d4*( t1 - x(11) )**2
          t( 2, 4, 8) = t( 2, 4, 8) - 2.0_rk*d2*d4*( t1 - x(11) )**2
          t( 3, 4, 8) = t( 3, 4, 8) - 2.0_rk*d3*d4*( t1 - x(11) )**2
          t( 4, 4, 8) = t( 4, 4, 8) - 4.0_rk*d4**2*( t1 - x(11) )**2
          t( 1, 5, 8) = t( 1, 5, 8) + 2.0_rk*t1*d1*d4*x(4)*( t1 - x(11) )**2
          t( 4, 5, 8) = t( 4, 5, 8) + 2.0_rk*t1*d1*d4*x(1)*( t1 - x(11) )**2
          t( 5, 5, 8) = t( 5, 5, 8) - 2.0_rk*t1**2*d1*d4*x(1)*x(4)*( t1 - x(11) )**2
          t( 2, 6, 8) = t( 2, 6, 8) + 2.0_rk*d2*d4*x(4)*( t1 - x(9) )**2*( t1 - x(11) )**2
          t( 4, 6, 8) = t( 4, 6, 8) + 2.0_rk*d2*d4*x(2)*( t1 - x(9) )**2*( t1 - x(11) )**2
          t( 6, 6, 8) = t( 6, 6, 8) - 2.0_rk*d2*d4*x(2)*x(4)*( t1 - x(9) )**4*( t1 - x(11) )**2
          t( 3, 7, 8) = t( 3, 7, 8) + 2.0_rk*d3*d4*x(4)*( t1 - x(10) )**2*( t1 - x(11) )**2
          t( 4, 7, 8) = t( 4, 7, 8) + 2.0_rk*d3*d4*x(3)*( t1 - x(10) )**2*( t1 - x(11) )**2
          t( 7, 7, 8) = t( 7, 7, 8) - 2.0_rk*d3*d4*x(3)*x(4)*( t1 - x(10) )**4*( t1 - x(11) )**2
          t( 1, 8, 8) = t( 1, 8, 8) + 2.0_rk*d1*d4*x(4)*( t1 - x(11) )**4
          t( 2, 8, 8) = t( 2, 8, 8) + 2.0_rk*d2*d4*x(4)*( t1 - x(11) )**4
          t( 3, 8, 8) = t( 3, 8, 8) + 2.0_rk*d3*d4*x(4)*( t1 - x(11) )**4
          t( 4, 8, 8) = t( 4, 8, 8) - 2.0_rk*d4*( t1 - x(11) )**4*( s1 - 3.0_rk*d4*x(4) )
          t( 5, 8, 8) = t( 5, 8, 8) - 2.0_rk*t1*d1*d4*x(1)*x(4)*( t1 - x(11) )**4
          t( 6, 8, 8) = t( 6, 8, 8) - 2.0_rk*d2*d4*x(2)*x(4)*( t1 - x(9) )**2*( t1 - x(11) )**4
          t( 7, 8, 8) = t( 7, 8, 8) - 2.0_rk*d3*d4*x(3)*x(4)*( t1 - x(10) )**2*( t1 - x(11) )**4
          t( 8, 8, 8) = t( 8, 8, 8) + 2.0_rk*d4*x(4)*( t1 - x(11) )**6*( s1 - 3.0_rk*d4*x(4) )
          t( 1, 2, 9) = t( 1, 2, 9) + 4.0_rk*d1*d2*x(6)*( t1 - x(9) )
          t( 2, 2, 9) = t( 2, 2, 9) + 8.0_rk*d2**2*x(6)*( t1 - x(9) )
          t( 2, 3, 9) = t( 2, 3, 9) + 4.0_rk*d2*d3*x(6)*( t1 - x(9) )
          t( 2, 4, 9) = t( 2, 4, 9) + 4.0_rk*d2*d4*x(6)*( t1 - x(9) )
          t( 1, 5, 9) = t( 1, 5, 9) - 4.0_rk*t1*d1*d2*x(2)*x(6)*( t1 - x(9) )
          t( 2, 5, 9) = t( 2, 5, 9) - 4.0_rk*t1*d1*d2*x(1)*x(6)*( t1 - x(9) )
          t( 5, 5, 9) = t( 5, 5, 9) + 4.0_rk*t1**2*d1*d2*x(1)*x(2)*x(6)*( t1 - x(9) )
          t( 1, 6, 9) = t( 1, 6, 9) - 4.0_rk*d1*d2*x(2)*( t1 - x(9) )*( x(6)*( t1 - x(9) )**2 - 1.0_rk )
          t( 2, 6, 9) = t( 2, 6, 9) + 4.0_rk*d2*( t1 - x(9) )*( &
               s1*( ( t1 - x(9) )**2*x(6) - 1.0_rk ) -          &
               2.0_rk*d2*x(2)*x(6)*( t1 - x(9) )**2 )
          t( 3, 6, 9) = t( 3, 6, 9) - 4.0_rk*d2*d3*x(2)*( t1 - x(9) )*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 4, 6, 9) = t( 4, 6, 9) - 4.0_rk*d2*d4*x(2)*( t1 - x(9) )*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 5, 6, 9) = t( 5, 6, 9) + 4.0_rk*t1*d1*d2*x(1)*x(2)*( t1 - x(9) )*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 6, 6, 9) = t( 6, 6, 9) - 4.0_rk*d2*x(2)*( t1 - x(9) )**3*( &
               s1*( ( t1 - x(9) )**2*x(6) - 2.0_rk ) -                  &
               d2*x(2)*x(6)*( t1 - x(9) )**2 + 2.0_rk*d2*x(2)*( ( t1 - x(9) )**2*x(6) - 1.0_rk ) )
          t( 2, 7, 9) = t( 2, 7, 9) - 4.0_rk*d2*d3*x(3)*x(6)*( t1 - x(9) )*( t1 - x(10) )**2
          t( 3, 7, 9) = t( 3, 7, 9) - 4.0_rk*d2*d3*x(2)*x(6)*( t1 - x(9) )*( t1 - x(10) )**2
          t( 6, 7, 9) = t( 6, 7, 9) + 4.0_rk*d2*d3*x(2)*x(3)*( t1 - x(9) )*( t1 - x(10) )**2*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 7, 7, 9) = t( 7, 7, 9) + 4.0_rk*d2*d3*x(2)*x(3)*x(6)*( t1 - x(9) )*( t1 - x(10) )**4
          t( 2, 8, 9) = t( 2, 8, 9) - 4.0_rk*d2*d4*x(4)*x(6)*( t1 - x(9) )*( t1 - x(11) )**2
          t( 4, 8, 9) = t( 4, 8, 9) - 4.0_rk*d2*d4*x(2)*x(6)*( t1 - x(9) )*( t1 - x(11) )**2
          t( 6, 8, 9) = t( 6, 8, 9) + 4.0_rk*d2*d4*x(2)*x(4)*( t1 - x(9) )*( t1 - x(11) )**2*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 8, 8, 9) = t( 8, 8, 9) + 4.0_rk*d2*d4*x(2)*x(4)*x(6)*( t1 - x(9) )*( t1 - x(11) )**4
          t( 1, 9, 9) = t( 1, 9, 9) + 4.0_rk*d1*d2*x(2)*x(6)*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 2, 9, 9) = t( 2, 9, 9) - 4.0_rk*d2*x(6)*(                      &
               ( s1 - d2*x(2) )*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk ) - &
               4.0_rk*d2*x(2)*x(6)*( t1 - x(9) )**2 )
          t( 3, 9, 9) = t( 3, 9, 9) + 4.0_rk*d2*d3*x(2)*x(6)*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 4, 9, 9) = t( 4, 9, 9) + 4.0_rk*d2*d4*x(2)*x(6)*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 5, 9, 9) = t( 5, 9, 9) - 4.0_rk*t1*d1*d2*x(1)*x(2)*x(6)*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 6, 9, 9) = t( 6, 9, 9) + 4.0_rk*(                                                                  &
               s1*d2*x(6)*( x(2)*( t1 - x(9) )**2*( 2.0_rk*( t1 - x(9) )**2*x(6) - 3.0_rk ) - ( t1 - x(9) ) ) - &
               ( d2*x(2) )**2*x(6)*( t1 - x(9) )**2*( 6.0_rk*( t1 - x(9) )**2*x(6) - 5.0_rk ) )
          t( 7, 9, 9) = t( 7, 9, 9) - 4.0_rk*d2*d3*x(2)*x(3)*x(6)*( t1 - x(10) )**2*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 8, 9, 9) = t( 8, 9, 9) - 4.0_rk*d2*d4*x(2)*x(4)*x(6)*( t1 - x(11) )**2*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 9, 9, 9) = t( 9, 9, 9) + 8.0_rk*d2*x(2)*x(6)**2*( t1 - x(9) )*( &
               s1*( 3.0_rk - 2.0_rk*( t1 - x(9) )**2*x(6) ) +                &
               3.0_rk*d2*x(2)*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk ) )
          t( 1, 3,10) = t( 1, 3,10) + 4.0_rk*d1*d3*x(7)*( t1 - x(10) )
          t( 2, 3,10) = t( 2, 3,10) + 4.0_rk*d2*d3*x(7)*( t1 - x(10) )
          t( 3, 3,10) = t( 3, 3,10) + 8.0_rk*d3**2*x(7)*( t1 - x(10) )
          t( 3, 4,10) = t( 3, 4,10) + 4.0_rk*d3*d4*x(7)*( t1 - x(10) )
          t( 1, 5,10) = t( 1, 5,10) - 4.0_rk*t1*d1*d3*x(3)*x(7)*( t1 - x(10) )
          t( 3, 5,10) = t( 3, 5,10) - 4.0_rk*t1*d1*d3*x(1)*x(7)*( t1 - x(10) )
          t( 5, 5,10) = t( 5, 5,10) + 4.0_rk*t1**2*d1*d3*x(1)*x(3)*x(7)*( t1 - x(10) )
          t( 2, 6,10) = t( 2, 6,10) - 4.0_rk*d2*d3*x(3)*x(7)*( t1 - x(9) )**2*( t1 - x(10) )
          t( 3, 6,10) = t( 3, 6,10) - 4.0_rk*d2*d3*x(2)*x(7)*( t1 - x(9) )**2*( t1 - x(10) )
          t( 6, 6,10) = t( 6, 6,10) + 4.0_rk*d2*d3*x(2)*x(3)*x(7)*( t1 - x(9) )**4*( t1 - x(10) )
          t( 1, 7,10) = t( 1, 7,10) - 4.0_rk*d1*d3*x(3)*( t1 - x(10) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 2, 7,10) = t( 2, 7,10) - 4.0_rk*d2*d3*x(3)*( t1 - x(10) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 3, 7,10) = t( 3, 7,10) + 4.0_rk*d3*( t1 - x(10) )*( &
               s1*( ( t1 - x(10) )**2*x(7) - 1.0_rk ) -          &
               3.0_rk*d3*x(3)*x(7)*( t1 - x(10) )**2 +           &
               d3*x(3) )
          t( 4, 7,10) = t( 4, 7,10) - 4.0_rk*d3*d4*x(3)*( t1 - x(10) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 5, 7,10) = t( 5, 7,10) + 4.0_rk*t1*d1*d3*x(1)*x(3)*( t1 - x(10) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 6, 7,10) = t( 6, 7,10) + 4.0_rk*d2*d3*x(2)*x(3)*( t1 - x(9) )**2*( t1 - x(10) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 7, 7,10) = t( 7, 7,10) + 4.0_rk*d3*x(3)*( t1 - x(10) )**3*( &
               - s1*( ( t1 - x(10) )**2*x(7) - 2.0_rk ) +                &
               d3*x(3)*x(7)*( t1 - x(10) )**2 +                           &
               2.0_rk*d3*x(3)*( ( t1 - x(10) )**2*x(7) - 1.0_rk ) )
          t( 3, 8,10) = t( 3, 8,10) - 4.0_rk*d3*d4*x(4)*x(7)*( t1 - x(10) )*( t1 - x(11) )**2
          t( 4, 8,10) = t( 4, 8,10) - 4.0_rk*d3*d4*x(3)*x(7)*( t1 - x(10) )*( t1 - x(11) )**2
          t( 7, 8,10) = t( 7, 8,10) + 4.0_rk*d3*d4*x(3)*x(4)*( t1 - x(10) )*( t1 - x(11) )**2*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 8, 8,10) = t( 8, 8,10) + 4.0_rk*d3*d4*x(3)*x(4)*x(7)*( t1 - x(10) )*( t1 - x(11) )**4
          t( 2, 9,10) = t( 2, 9,10) + 8.0_rk*d2*d3*x(3)*x(6)*x(7)*( t1 - x(9) )*( t1 - x(10) )
          t( 3, 9,10) = t( 3, 9,10) + 8.0_rk*d2*d3*x(2)*x(6)*x(7)*( t1 - x(9) )*( t1 - x(10) )
          t( 6, 9,10) = t( 6, 9,10) - 8.0_rk*d2*d3*x(2)*x(3)*x(7)*( t1 - x(9) )*( t1 - x(10) )*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 7, 9,10) = t( 7, 9,10) - 8.0_rk*d2*d3*x(2)*x(3)*x(6)*( t1 - x(9) )*( t1 - x(10) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 9, 9,10) = t( 9, 9,10) + 8.0_rk*d2*d3*x(2)*x(3)*x(6)*x(7)*( t1 - x(10) )*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 1,10,10) = t( 1,10,10) + 4.0_rk*d1*d3*x(3)*x(7)*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 2,10,10) = t( 2,10,10) + 4.0_rk*d2*d3*x(3)*x(7)*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 3,10,10) = t( 3,10,10) - 4.0_rk*d3*x(7)*(                       &
               ( s1 - d3*x(3) )*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk ) - &
               4.0_rk*d3*x(3)*x(7)*( t1 - x(10) )**2 )
          t( 4,10,10) = t( 4,10,10) + 4.0_rk*d3*d4*x(3)*x(7)*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 5,10,10) = t( 5,10,10) - 4.0_rk*t1*d1*d3*x(1)*x(3)*x(7)*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 6,10,10) = t( 6,10,10) - 4.0_rk*d2*d3*x(2)*x(3)*x(7)*( t1 - x(9) )**2*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 7,10,10) = t( 7,10,10) + 2.0_rk*( &
               s1*( 2.0_rk*d3*x(3)*( x(7)*( t1 - x(10) )**2*( 2.0_rk*x(7)*( t1 - x(10) )**2 - 3.0_rk ) - &
               2.0_rk*x(3)*( 2.0_rk*x(7)*( t1 - x(10) )**2 - 1.0_rk ) ) ) -                              &
               2.0_rk*( d3*x(3) )**2*x(7)*( t1 - x(10) )**2*( 6.0_rk*( t1 - x(10) )**2*x(7) - 5.0_rk ) )
          t( 8,10,10) = t( 8,10,10) - 4.0_rk*d3*d4*x(3)*x(4)*x(7)*( t1 - x(11) )**2*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 9,10,10) = t( 9,10,10) + 8.0_rk*d2*d3*x(2)*x(3)*x(6)*x(7)*( t1 - x(9) )*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t(10,10,10) = t(10,10,10) + 8.0_rk*d3*x(3)*x(7)**2*( t1 - x(10) )*( &
               s1*( 3.0_rk - 2.0_rk*( t1 - x(10) )**2*x(7) ) +                &
               3.0_rk*d3*x(3)*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk ) )
          t( 1, 4,11) = t( 1, 4,11) + 4.0_rk*d1*d4*x(8)*( t1 - x(11) )
          t( 2, 4,11) = t( 2, 4,11) + 4.0_rk*d2*d4*x(8)*( t1 - x(11) )
          t( 3, 4,11) = t( 3, 4,11) + 4.0_rk*d3*d4*x(8)*( t1 - x(11) )
          t( 4, 4,11) = t( 4, 4,11) + 8.0_rk*d4**2*x(8)*( t1 - x(11) )
          t( 1, 5,11) = t( 1, 5,11) - 4.0_rk*t1*d1*d4*x(4)*x(8)*( t1 - x(11) )
          t( 4, 5,11) = t( 4, 5,11) - 4.0_rk*t1*d1*d4*x(1)*x(8)*( t1 - x(11) )
          t( 5, 5,11) = t( 5, 5,11) + 4.0_rk*t1**2*d1*d4*x(1)*x(4)*x(8)*( t1 - x(11) )
          t( 2, 6,11) = t( 2, 6,11) - 4.0_rk*d2*d4*x(4)*x(8)*( t1 - x(9) )**2*( t1 - x(11) )
          t( 4, 6,11) = t( 4, 6,11) - 4.0_rk*d2*d4*x(2)*x(8)*( t1 - x(9) )**2*( t1 - x(11) )
          t( 6, 6,11) = t( 6, 6,11) + 4.0_rk*d2*d4*x(2)*x(4)*x(8)*( t1 - x(9) )**4*( t1 - x(11) )
          t( 3, 7,11) = t( 3, 7,11) - 4.0_rk*d3*d4*x(4)*x(8)*( t1 - x(10) )**2*( t1 - x(11) )
          t( 4, 7,11) = t( 4, 7,11) - 4.0_rk*d3*d4*x(3)*x(8)*( t1 - x(10) )**2*( t1 - x(11) )
          t( 7, 7,11) = t( 7, 7,11) + 4.0_rk*d3*d4*x(3)*x(4)*x(8)*( t1 - x(10) )**4*( t1 - x(11) )
          t( 1, 8,11) = t( 1, 8,11) - 4.0_rk*d1*d4*x(4)*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 2, 8,11) = t( 2, 8,11) - 4.0_rk*d2*d4*x(4)*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 3, 8,11) = t( 3, 8,11) - 4.0_rk*d3*d4*x(4)*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 4, 8,11) = t( 4, 8,11) + 4.0_rk*d4*( t1 - x(11) )*( &
               ( s1 - d4*x(4) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )- &
               2.0_rk*d4*x(4)*x(8)*( t1 - x(11) )**2 )
          t( 5, 8,11) = t( 5, 8,11) + 4.0_rk*t1*d1*d4*x(1)*x(4)*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 6, 8,11) = t( 6, 8,11) + 4.0_rk*d2*d4*x(2)*x(4)*( t1 - x(9) )**2*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 7, 8,11) = t( 7, 8,11) + 4.0_rk*d3*d4*x(3)*x(4)*( t1 - x(10) )**2*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 8, 8,11) = t( 8, 8,11) - 4.0_rk*d4*x(4)*( t1 - x(11) )**3*( &
               s1*( x(8)*( t1 - x(11) )**2 - 2.0_rk ) -                  &
               d4*x(4)*x(8)*( t1 - x(11) )**2 -                          &
               2.0_rk*d4*x(4)*( ( t1 - x(11) )**2*x(8) - 1.0_rk ) )
          t( 2, 9,11) = t( 2, 9,11) + 8.0_rk*d2*d4*x(4)*x(6)*x(8)*( t1 - x(9) )*( t1 - x(11) )
          t( 4, 9,11) = t( 4, 9,11) + 8.0_rk*d2*d4*x(2)*x(6)*x(8)*( t1 - x(9) )*( t1 - x(11) )
          t( 6, 9,11) = t( 6, 9,11) - 8.0_rk*d2*d4*x(2)*x(4)*x(8)*( t1 - x(9) )*( t1 - x(11) )*( ( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 8, 9,11) = t( 8, 9,11) - 8.0_rk*d2*d4*x(2)*x(4)*x(6)*( t1 - x(9) )*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 9, 9,11) = t( 9, 9,11) + 8.0_rk*d2*d4*x(2)*x(4)*x(6)*x(8)*( t1 - x(11) )*( 2.0_rk*( t1 - x(9) )**2*x(6) - 1.0_rk )
          t( 3,10,11) = t( 3,10,11) + 8.0_rk*d3*d4*x(4)*x(7)*x(8)*( t1 - x(10) )*( t1 - x(11) )
          t( 4,10,11) = t( 4,10,11) + 8.0_rk*d3*d4*x(3)*x(7)*x(8)*( t1 - x(10) )*( t1 - x(11) )
          t( 7,10,11) = t( 7,10,11) - 8.0_rk*d3*d4*x(3)*x(4)*x(8)*( t1 - x(10) )*( t1 - x(11) )*( ( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 8,10,11) = t( 8,10,11) - 8.0_rk*d3*d4*x(3)*x(4)*x(7)*( t1 - x(10) )*( t1 - x(11) )*( ( t1 - x(11) )**2*x(8) - 1.0_rk )
          t(10,10,11) = t(10,10,11) + 8.0_rk*d3*d4*x(3)*x(4)*x(7)*x(8)*( t1 - x(11) )*( 2.0_rk*( t1 - x(10) )**2*x(7) - 1.0_rk )
          t( 1,11,11) = t( 1,11,11) + 4.0_rk*d1*d4*x(4)*x(8)*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 2,11,11) = t( 2,11,11) + 4.0_rk*d2*d4*x(4)*x(8)*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 3,11,11) = t( 3,11,11) + 4.0_rk*d3*d4*x(4)*x(8)*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 4,11,11) = t( 4,11,11) + 4.0_rk*d4*x(8)*( &
               ( d4*x(4) - s1 )*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk ) + &
               4.0_rk*d4*x(4)*x(8)*( t1 - x(11) )**2 )
          t( 5,11,11) = t( 5,11,11) - 4.0_rk*t1*d1*d4*x(1)*x(4)*x(8)*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk ) 
          t( 6,11,11) = t( 6,11,11) - 4.0_rk*d2*d4*x(2)*x(4)*x(8)*( t1 - x(9) )**2*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 7,11,11) = t( 7,11,11) - 4.0_rk*d3*d4*x(3)*x(4)*x(8)*( t1 - x(10) )**2*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t( 8,11,11) = t( 8,11,11) + 4.0_rk*d4*x(4)*(                                  &
               s1*( x(8)*( t1 - x(11) )**2*( 2.0_rk*x(8)*( t1 - x(11) )**2 - 3.0_rk ) - &
               ( 2.0_rk*x(8)*( t1 - x(11) )**2 - 1.0_rk ) ) -                           &
               d4*x(4)*x(8)*( t1 - x(11) )**2*( 6.0_rk*( t1 - x(11) )**2*x(8) - 5.0_rk ) )
          t( 9,11,11) = t( 9,11,11) + 8.0_rk*d2*d4*x(2)*x(4)*x(6)*x(8)*( t1 - x(9) )*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t(10,11,11) = t(10,11,11) + 8.0_rk*d3*d4*x(3)*x(4)*x(7)*x(8)*( t1 - x(10) )*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk )
          t(11,11,11) = t(11,11,11) + 8.0_rk*d4*x(4)*x(8)**2*( t1 - x(11) )*( &
               s1*( 3.0_rk - 2.0_rk*( t1 - x(11) )**2*x(8) ) +                &
               3.0_rk*d4*x(4)*( 2.0_rk*( t1 - x(11) )**2*x(8) - 1.0_rk ) )
       end do

    case ( 20 ) ! Watson
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do
       do i = 1, 29
          d1 = real( i, rk )/29.0_rk
          s1 = 0.0_rk
          do j = 1, global_n
             s1 = s1 + x(j)*d1**( j - 1 )
          end do
          do j = 1, global_n
             w1(j) = real( j-1, rk )*d1**( j-2 ) - 2.0_rk*d1**( j-1 )*s1 ! g_i(j)
          end do
          do l = 1, global_n
             do k = 1, l
                t2 = - 2.0_rk*d1**( l+k-2 ) ! h_i(k,l)
                do j = 1, k-1
                   t(j,k,l) = t(j,k,l) + 2.0_rk*w1(j)*t2
                end do
                do j = k+1,l-1
                   t(k,j,l) = t(k,j,l) + 2.0_rk*w1(j)*t2
                end do
                do j = l+1, global_n
                   t(k,l,j) = t(k,l,j) + 2.0_rk*w1(j)*t2
                end do
                if ( k .ne. l ) then
                   t(k,k,l) = t(k,k,l) + 4.0_rk*w1(k)*t2
                   t(k,l,l) = t(k,l,l) + 4.0_rk*w1(l)*t2
                else
                   t(l,l,l) = t(l,l,l) + 6.0_rk*w1(l)*t2
                end if
             end do
          end do
       end do
       t(1,1,1) = t(1,1,1) + 24.0_rk*x(1)
       t(1,1,2) = t(1,1,2) - 4.0_rk

    case ( 21 ) ! Extended Rosenbrock
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do j = 1, global_n, 2
          t(j,j,j  ) =   2.4e+3_rk*x(j)
          t(j,j,j+1) = - 4.0e+2_rk
       end do

    case ( 22 ) ! Extended Powell singular
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_n, 4
          t(i  ,i  ,i  ) =   2.40e+2_rk*( x(i) - x(i+3) )
          t(i+1,i+1,i+1) =   2.40e+1_rk*( x(i+1) - 2.0_rk*x(i+2) )
          t(i+1,i+1,i+2) = - 4.80e+1_rk*( x(i+1) - 2.0_rk*x(i+2) )
          t(i+1,i+2,i+2) =   9.60e+1_rk*( x(i+1) - 2.0_rk*x(i+2) )
          t(i+2,i+2,i+2) = - 1.92e+2_rk*( x(i+1) - 2.0_rk*x(i+2) )
          t(i  ,i  ,i+3) = - 2.40e+2_rk*( x(i) - x(i+3) )
          t(i  ,i+3,i+3) =   2.40e+2_rk*( x(i) - x(i+3) )
          t(i+3,i+3,i+3) = - 2.40e+2_rk*( x(i) - x(i+3) )
       end do

    case ( 23 ) ! Penalty I
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do j = 1, global_n
          do k = 1, j
             t(k,j,j) = 8.0_rk*x(k)
             t(k,k,j) = 8.0_rk*x(j)
          end do
          t(j,j,j) = 24.0_rk*x(j)
       end do

    case ( 24 ) ! Penalty II
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       d1 = exp( 0.1_rk )
       d2 = 1.0_rk
       s2 = exp( x(1)/10.0_rk )

       t1 = exp( x(2)/10.0_rk ) + s2 - d2*( d1 + 1.0_rk )
       t(1,1,1) = 2.0_rk*( 1.0e-8_rk*( t1*s2 + 3.0_rk*s2**2 ) + 12.0_rk*real( global_n, rk )**2*x(1) )

       do j = 2, global_n
          s1 = exp( x(j)/10.0_rk )
          t1 = s1 + s2 - d2*( d1 + 1.0_rk ) ! f_{j}
          t2 = s1 - 1.0_rk/d1               ! f_{j+n-1}

          do k = j-1, 2, -1
             t(j-k,j-k,j) = 8.0_rk*real( global_n-j+k+1, rk )*real( global_n-j+1, rk )*x(j)
             t(j-k,j,  j) = 8.0_rk*real( global_n-j+k+1, rk )*real( global_n-j+1, rk )*x(j-k)
          end do

          t(j-1,j-1,j) = 2.0_rk*( 1.0e-8_rk*s1*s2 + 4.0_rk*real( global_n-j+2, rk )*real( global_n-j+1, rk )*x(j) )
          t(j-1,j  ,j) = 2.0_rk*( 1.0e-8_rk*s1*s2 + 4.0_rk*real( global_n-j+2, rk )*real( global_n-j+1, rk )*x(j-1) )
          t(j,  j  ,j) = 2.0_rk*( 1.0e-8_rk*( t1*s1 + 6.0_rk*s1**2 + t2*s1 ) + 12.0_rk*real( global_n-j+1, rk )**2*x(j) )

          s2 = s1
          d2 = d1*d2
       end do

    case ( 25 ) ! Variably dimensioned
       s1 = 0.0_rk
       do i = 1, global_n
          s1 = s1 + i*( x(i) - 1.0_rk )
       end do

       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 24.0_rk*real( l*k*j, rk )*s1
             end do
          end do
       end do

    case ( 26 ) ! Trigonometric
       s1 = 0.0_rk
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
          s1 = s1 + cos( x(j) )
       end do

       do i = 1, global_m
          t3 = real( global_n+i, rk ) - s1 - real( i, rk )*cos( x(i) ) - sin( x(i) )
          do j = 1, global_n
             d1 = sin( x(j) )
             d2 = cos( x(j) )
             if ( j .eq. i ) then
                t(i,i,i) = t(i,i,i) + 2.0_rk*( t3*( d2 - real( 1+i, rk )*d1 ) + &
                     3.0_rk*( real( 1+i, rk )*d1 - d2 )*( real( 1+i, rk )*d2 + d1 ) )
                do k = j+1, global_n
                   t1 = sin( x(k) )
                   t2 = cos( x(k) )
                   t(i,i,k) = t(i,i,k) + 2.0_rk*t1*( real( 1+i, rk )*d2 + d1 )
                   t(i,k,k) = t(i,k,k) + 2.0_rk*t2*( real( 1+i, rk )*d1 - d2 )
                end do
             else
                t(j,j,j) = t(j,j,j) + 2.0_rk*( 3.0_rk*d1*d2 - t3*d1 )
                do k = j+1, global_n
                   t1 = sin( x(k) )
                   t2 = cos( x(k) )

                   if ( k .eq. i ) then
                      t(j,j,k) = t(j,j,k) + 2.0_rk*( real( 1+i, rk )*t1 - t2 )*d2
                      t(j,k,k) = t(j,k,k) + 2.0_rk*( real( 1+i, rk )*t2 + t1 )*d1
                   else
                      t(j,j,k) = t(j,j,k) + 2.0_rk*t1*d2
                      t(j,k,k) = t(j,k,k) + 2.0_rk*t2*d1
                   end if
                end do
             end if
          end do
       end do

    case ( 27 ) ! Brown almost-linear
       do j = 1, global_n
          t(j,j,j) = 0.0_rk
          prod(j-1,j) = 1.0_rk
          do k = j, global_n
             prod(k,j) = prod(k-1,j)*x(k)
          end do
          prod(j,global_n+1) = 1.0_rk
       end do

       t1 = prod(global_n,1) - 1.0_rk
       do k = 1, global_n
          do j = 1, k-1
             t(j,j,k) = 4.0_rk*prod(j-1,1)*prod(global_n,j+1) * prod(j-1,1)*prod(k-1,j+1)*prod(global_n,k+1)
             t(j,k,k) = 4.0_rk*prod(k-1,1)*prod(global_n,k+1) * prod(j-1,1)*prod(k-1,j+1)*prod(global_n,k+1)
             do i = 1, j-1
                t(i,j,k) = 2.0_rk*( &
                     t1*prod(i-1,1)*prod(j-1,i+1)*prod(k-1,j+1)*prod(global_n,k+1) +                 &
                     prod(i-1,1)*prod(global_n,i+1) * prod(j-1,1)*prod(k-1,j+1)*prod(global_n,k+1) + &
                     prod(j-1,1)*prod(global_n,j+1) * prod(i-1,1)*prod(k-1,i+1)*prod(global_n,k+1) + &
                     prod(k-1,1)*prod(global_n,k+1) * prod(i-1,1)*prod(j-1,i+1)*prod(global_n,j+1) )
             end do
          end do
       end do

    case ( 28 ) ! Discrete boundary value
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       d1 = 1.0_rk / real( global_n+1, rk )

       t1 = 2*x(1) - x(2) + 0.5_rk*d1**2*( x(1) + d1 + 1.0_rk )**3
       t(1,1,1) = 6.0_rk*d1**2*( t1 + 6.0_rk*( x(1) + d1 + 1.0_rk ) + 4.5_rk*d1**2*( x(1) + d1 + 1.0_rk)**3 )
       t(1,1,2) = - 6.0_rk*d1**2*( x(1) + d1 + 1.0_rk )

       t1 = 2*x(global_n) - x(global_n-1) + 0.5_rk*d1**2*( x(global_n) + global_n*d1 + 1.0_rk )**3
       t(global_n-1,global_n,global_n) = - 6.0_rk*d1**2*( x(global_n) + global_n*d1 + 1.0_rk )
       t(global_n,  global_n,global_n) = 6.0_rk*d1**2*( &
            t1 + 6.0_rk*( x(global_n) + global_n*d1 + 1.0_rk ) + 4.5_rk*d1**2*( x(global_n) + global_n*d1 + 1.0_rk)**3 )

       do i = 2, global_n-1
          t1 = 2*x(i) - x(i-1) - x(i+1) + 0.5_rk*d1**2*( x(i) + i*d1 + 1.0_rk )**3

          t(i-1,i,i  ) = - 6.0_rk*d1**2*( x(i) + i*d1 + 1.0_rk )
          t(  i,i,i  ) = 6.0_rk*d1**2*( t1 + 6.0_rk*( x(i) + i*d1 + 1.0_rk ) + 4.5_rk*d1**2*( x(i) + i*d1 + 1.0_rk)**3 )
          t(  i,i,i+1) = - 6.0_rk*d1**2*( x(i) + i*d1 + 1.0_rk )
       end do

    case ( 29 ) ! Discrete integral equation
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       d1 = 1.0_rk/( real( global_n, rk ) + 1.0_rk )
       w1(1)          = d1*( x(1) + d1 + 1.0_rk )**3
       w2(global_n)   = ( 1.0_rk - global_n*d1 )*( x(global_n) + global_n*d1 + 1.0_rk )**3
       w2(global_n+1) = 0.0_rk
       do i = 2, global_n
          t1 = real( i, rk )*d1
          t2 = real( global_n-i+1, rk )*d1
          w1(i)            = w1(i-1)          +              t1*( x(i)            + t1 + 1.0_rk )**3
          w2(global_n-i+1) = w2(global_n-i+2) + ( 1.0_rk - t2 )*( x(global_n-i+1) + t2 + 1.0_rk )**3
       end do

       do i = 1, global_n
          t1 = real( i, rk )*d1
          t3 = x(i) + 0.5_rk*d1*( ( 1.0_rk - t1 )*w1(i) + t1*w2(i+1) )

          do j = 1, i
             t2 = real( j, rk )*d1
             w3(j) = 1.5_rk*d1*( 1.0_rk - t1 )*t2*( x(j) + t2 + 1.0_rk )**2
             w4(j) = 3.0_rk*d1*( 1.0_rk - t1 )*t2*( x(j) + t2 + 1.0_rk )
             t(j,j,j) = t(j,j,j) + 2.0_rk*t3*( 3.0_rk*d1*( 1.0_rk - t1 )*t2 )
          end do

          do j = i+1, global_n
             t2 = real( j, rk )*d1
             w3(j) = 1.5_rk*d1*t1*( 1.0_rk - t2 )*( x(j) + t2 + 1.0_rk )**2
             w4(j) = 3.0_rk*d1*t1*( 1.0_rk - t2 )*( x(j) + t2 + 1.0_rk )
             t(j,j,j) = t(j,j,j) + 2.0_rk*t3*( 3.0_rk*d1*t1*( 1.0_rk - t2 ) )
          end do

          w3(i) = w3(i) + 1.0_rk

          do j = 1, global_n
             t(j,j,j) = t(j,j,j) + 6.0_rk*w4(j)*w3(j)
             do l = 1, j-1
                t(l,j,j) = t(l,j,j) + 2.0_rk*w4(j)*w3(l)
                t(l,l,j) = t(l,l,j) + 2.0_rk*w4(l)*w3(j)
             end do
          end do
       end do

    case ( 30 ) ! Broyden tridiagonal
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       t(1,1,1) = 24.0_rk*( 4.0_rk*x(1) - 3.0_rk )
       t(1,1,2) = 16.0_rk

       t(global_n-1,global_n,global_n) = 8.0_rk
       t(global_n,  global_n,global_n) = 24.0_rk*( 4.0_rk*x(global_n) - 3.0_rk )

       do i = 2, global_n-1
          t(i-1,i,i  ) = t(i-1,i,i  ) +  8.0_rk
          t(i,  i,i  ) = t(i,  i,i  ) + 24.0_rk*( 4.0_rk*x(i) - 3.0_rk )
          t(i,  i,i+1) = t(i,  i,i+1) + 16.0_rk
       end do

    case ( 31 ) ! Broyden banded
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_n
          s1 = 0.0_rk
          do j = max( 1, i-5 ), i-1
             s1 = s1 + x(j)*( 1.0_rk + x(j) )
          end do
          if ( i .ne. global_n ) s1 = s1 + x(i+1)*( 1.0_rk + x(i+1) )
          
          do j = max(1,i-5), i-1
             t(j,j,j) = t(j,j,j) + 12.0_rk*( 1.0_rk + 2.0_rk*x(j) )
             do l = j+1, i-1
                t(j,j,l) = t(j,j,l) + 4.0_rk*( 1.0_rk + 2.0_rk*x(l) )
                t(j,l,l) = t(j,l,l) + 4.0_rk*( 1.0_rk + 2.0_rk*x(j) )
             end do
             t(j,j,i) = t(j,j,i) - 4.0_rk*( 2.0_rk + 15.0_rk*x(i)**2 )
             t(j,i,i) = t(j,i,i) - 60.0_rk*x(i)*( 1.0_rk + 2.0_rk*x(j) )
             if ( i .ne. global_n ) then
                t(j,  j,i+1) = t(j,  j,i+1) + 4.0_rk*( 1.0_rk + 2.0_rk*x(i+1) )
                t(j,i+1,i+1) = t(j,i+1,i+1) + 4.0_rk*( 1.0_rk + 2.0_rk*x(j) )
             end if
          end do

          t1 = x(i)*( 2.0_rk + 5.0_rk*x(i)**2 ) + 1.0_rk - s1
          d1 = 2.0_rk + 15.0_rk*x(i)**2
          d2 = 30.0_rk*x(i)
          t(i,i,i) = t(i,i,i) + 6.0_rk*( 10.0_rk*t1 + d1*d2 )

          if ( i .ne. global_n ) then
             t(i,i,i+1) = t(i,i,i+1) - 60.0_rk*x(i)*( 1.0_rk + 2.0_rk*x(i+1) )
             t(i,i+1,i+1) = t(i,i+1,i+1) - 4.0_rk*( 2.0_rk + 15.0_rk*x(i)**2 )
             t(i+1,i+1,i+1) = t(i+1,i+1,i+1) + 12.0_rk*( 1.0_rk + 2.0_rk*x(i+1) )
          end if
       end do

    case ( 32 ) ! Linear function - full rank
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

    case ( 33 ) ! Linear function - rank 1
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

    case ( 34 ) ! Linear function - rank 1 with zero columns and rows
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

    case ( 35 ) ! Chebyquad
       do j = 1, global_n
          do k = 1, j
             do l = 1, k
                t(l,k,j) = 0.0_rk
             end do
          end do
       end do

       do i = 1, global_m
          s1 = 0
          do j = 1, global_n
             t1 = 1.0_rk
             t2 = 2.0_rk*x(j) - 1.0_rk

             d1 = 0.0_rk
             d2 = 2.0_rk

             h1 = 0.0_rk
             h2 = 0.0_rk

             td1 = 0.0_rk
             td2 = 0.0_rk

             t4 = 2.0_rk*( 2.0_rk*x(j) - 1.0_rk )

             t3  = t2
             d3  = d2
             h3  = h2
             td3 = td2
             do k = 2, i
                td3 = t4*td2 + 12.0_rk*h2 - td1
                td1 = td2
                td2 = td3

                h3 = t4*h2 + 8.0_rk*d2 - h1
                h1 = h2
                h2 = h3

                d3 = t4*d2 + 4.0_rk*t2 - d1
                d1 = d2
                d2 = d3

                t3 = t4*t2 - t1
                t1 = t2
                t2 = t3
             end do
             
             s1    = s1 + t3                           ! f_i
             w1(j) = 1.0_rk / real( global_n, rk )*d3  ! g_i(j)
             w2(j) = 1.0_rk / real( global_n, rk )*h3  ! h_i(j,j)
             w3(j) = 1.0_rk / real( global_n, rk )*td3 ! t_i(j,j,j)
          end do

          s1 = s1 / real( global_n, rk ) ! f_i
          if ( modulo( i, 2 ) .eq. 0 ) s1 = s1 + 1.0_rk / ( real( i, rk ) ** 2 - 1.0_rk )

          do j = 1, global_n
             t(j,j,j) = t(j,j,j) + 2.0_rk*( s1*w3(j) + 3.0_rk*w1(j)*w2(j) )
             do l = 1, j-1
                t(l,j,j) = t(l,j,j) + 2.0_rk*w1(l)*w2(j)
                t(l,l,j) = t(l,l,j) + 2.0_rk*w1(j)*w2(l)
             end do
          end do
       end do
       
    case default
       flag = - 1
       
    end select

  end subroutine mgh_evalt

end module mgh
