module arp

  use arp_model
  use arp_types

  implicit none

  private

  public :: arp_load, arp_output, arp_param, arp_solve, sparse_matrix, &
            sparse_tensor

contains

  ! ==================================================================
  
  subroutine arp_load( param )
    ! STRUCTURE ARGUMENT
    type(arp_param), intent(out) :: param

    param%alpha      = 1.0e-08_dp
    param%bigJ       = 20
    param%dense      = .true.
    param%epsilon    = 1.0e-08_dp
    param%epsilonrel = 1.0e-09_dp
    param%eta1       = 1.0e+03_dp
    param%eta2       = 3.0_dp
    param%gamma1     = 5.0e-01_dp
    param%gamma2     = 1.0e+01_dp
    param%hnnzmax    = 100000
    param%kmax       = 500
    param%p          = 3
    param%printlevel = 1
    param%sigmaini   = .true.
    param%sigmalow   = 1.0e-08_dp
    param%theta      = 1.0e+02_dp
    param%tnnzmax    = 500000
    param%trackunit  = 0

  end subroutine arp_load

  ! ==================================================================

  subroutine arp_solve( x, param, output )

    ! ARRAY ARGUMENT
    real(kind=dp), dimension(:), intent(inout) :: x

    ! STRUCTURE ARGUMENT
    type(arp_param), intent(in) :: param
    
    ! ROUTINE OUTPUT
    type(arp_output), intent(out) :: output

    ! LOCAL SCALARS
    character(len=80)     :: specfnm, outputfnm
    logical               :: checkder, successful
    integer               :: allocstat, cnt, cons, fevalcnt, geninfo,      &
                             gevalcnt, hevalcnt, hnnzmax, flag, inform,    &
                             j, jcnnzmax, k, n, nvparam, nwstep, tevalcnt, &
                             uflag
    real(kind=4)          :: tfinish, tstart
    real(kind=dp)         :: cnorm, efacc, efstain, eoacc, eostain,   &
                             epsfeas, epsopt, ftrial, ginfnorm, mval, &
                             mgnorm2, nlpsupn, sinfnorm, snorm,       &
                             snorm2, sigmaini, xinfnorm, g0infnorm
    real(kind=dp), target :: f, sigma

    ! LOCAL ARRAYS
    character(len=80), dimension(10) :: vparam
    logical,           dimension(0)  :: equatn, linear
    logical,           dimension(11) :: coded
    real(kind=dp),     dimension(0)  :: lambda

    ! LOCAL ALLOCATABLE ARRAYS
    real(kind=dp), allocatable, dimension(:)             :: l, mg, u, s, xtrial
    real(kind=dp), allocatable, dimension(:),     target :: g
    real(kind=dp), allocatable, dimension(:,:),   target :: h
    real(kind=dp), allocatable, dimension(:,:,:), target :: t

    ! LOCAL SPARSE ARRAYS
    type(sparse_matrix), target :: hs
    type(sparse_tensor), target :: ts

    ! EXTERNAL ROUTINES
    external :: algencan

    ! Initializations
    call cpu_time( tstart )
    
    n = size( x )
    
    k = 0

    fevalcnt = 0
    gevalcnt = 0
    hevalcnt = 0
    tevalcnt = 0

    sigmaini = param%sigmalow

    nwstep = 0

    ! Allocate arrays
    allocate( g(n), l(n), mg(n), u(n), s(n), xtrial(n), stat=allocstat )
    if ( allocstat .ne. 0 ) then
       output%info = - 10
       return
    end if
    
    if ( param%dense ) then
       if ( param%p .ge. 2 ) then
          allocate( h(n,n), stat=allocstat )

          if ( param%p .eq. 3 ) then
             allocate( t(n,n,n), stat=allocstat )
          end if
          
          if ( allocstat .ne. 0 ) then
             output%info = - 11
             return
          end if
       end if
    else
       allocate( hs%row(param%hnnzmax), hs%col(param%hnnzmax), &
            hs%val(param%hnnzmax), ts%id1(param%tnnzmax),      &
            ts%id2(param%tnnzmax), ts%id3(param%tnnzmax),      &
            ts%val(param%tnnzmax), stat=allocstat )
    end if
       
    if ( allocstat .ne. 0 ) then
       output%info = - 12
       return
    end if

    ! Set Algencan parameters
    l(1:n) = - huge( 0.0_dp )
    u(1:n) =   huge( 0.0_dp )

    cons = 0

    coded(1:3)  = .true.
    coded(4:11) = .false.

    jcnnzmax = 0
    hnnzmax  = 5 * n * n

    checkder = .false.

    epsfeas   = 1.0d-08
    epsopt    = huge( 0.0d0 )
    efstain   = sqrt( epsfeas )
    eostain   = epsopt ** 1.5d0
    efacc     = sqrt( epsfeas )
    eoacc     = sqrt( epsopt )
    outputfnm = ''
    specfnm   = 'algencan.dat'

    nvparam = 0

    ! Set pointers for model functional evaluation subroutines
    if ( param%dense ) then
       call set_model_pointers( param%p, sigma, param%dense, f, g, h_in=h, t_in=t )
    else
       call set_model_pointers( param%p, sigma, param%dense, f, g, hs_in=hs, ts_in=ts )
    end if
    
    ! Evaluates the function and its derivatives at the initial point
    call param%uevalf( x, f, uflag )
    fevalcnt = fevalcnt + 1
    if ( uflag .ne. 0 ) then
       output%info = - 90
       return
    end if

    call param%uevalg( x, g, uflag )
    gevalcnt = gevalcnt + 1
    if ( uflag .ne. 0 ) then
       output%info = - 91
       return
    end if
    ginfnorm = maxval( abs( g ) )
    g0infnorm = max( 1.0_dp, ginfnorm )

    if ( param%p .ge. 2 ) then
       if ( param%dense ) then
          call param%uevalh ( x, h,  uflag )
       else
          call param%uevalhs( x, hs, uflag )
       end if
       
       hevalcnt = hevalcnt + 1

       if ( uflag .ne. 0 ) then
          output%info = - 92
          return
       end if
    end if

    if ( param%p .ge. 3 ) then
       if ( param%dense ) then
          call param%uevalt ( x, t,  uflag )
       else
          call param%uevalts( x, ts, uflag )
       end if
       
       tevalcnt = tevalcnt + 1

       if ( uflag .ne. 0 ) then
          output%info = - 93
          return
       end if
    end if

    if ( param%printlevel .gt. 0 ) then
       write(*,010) param%p, param%alpha, param%epsilon, param%eta1,  &
            param%eta2, param%gamma1, param%gamma2, param%hnnzmax,    &
            param%bigJ, param%kmax, param%printlevel, param%sigmaini, &
            param%sigmalow, param%theta, param%tnnzmax,               &
            param%trackunit
    end if


    if ( param%printlevel .eq. 1 ) then
       if ( modulo( k, 20 ) .eq. 0 ) then
          write (*,100) 'iter', 'sigma', 'f', '||g||', 'geninfo', '||mg||', &
               '||s||', 'ftrial', 'dec?'
       end if
    end if

    if ( param % printlevel .gt. 1 ) then
       write (*,fmt='(/,1P,A,D24.16,/)') "f0 = ", f
       write (*,fmt='(I0,A,/,2X,1P,D16.8)',advance='no') &
            min( n, 4 ), " coordinate(s) of x0:", x(1)
       do cnt = 2, min( n, 4 )
          write (*,fmt='(1P,2X,D16.8)',advance='no') x(cnt)
       end do
       write (*,*)
    end if

    outer: do
       ! --------------------------------------------------------------
       ! STEP 1: Initialization
       ! --------------------------------------------------------------
       j     = 0
       sigma = 0.0_dp

       if ( param % printlevel .gt. 1 ) then
          write (*,fmt="(/,136('-'),/,57X,A,I0,/)") &
               "Starting Iteration ", k

          write (*,fmt='(/,A)') "Testing stopping criteria"
          write (*,fmt='(1P,2(2X,A,D24.16,/))')          &
               "             ginfnorm = ", ginfnorm,     &
               "              epsilon = ", param%epsilon
       end if
       ! --------------------------------------------------------------

       ! --------------------------------------------------------------
       ! STOPPING CRITERIA
       ! --------------------------------------------------------------
       ! Norm of the gradient
       if ( ginfnorm .le. param%epsilon ) then
          flag = 0
          exit outer
       end if

       ! Norm of the gradient
       ! if ( ginfnorm .le. param%epsilonrel*g0infnorm ) then
       !    flag = 1
       !    exit outer
       ! end if

       ! Maximum number of iterations
       if ( k .ge. param%kmax ) then
          flag = - 2
          exit outer
       end if
       ! --------------------------------------------------------------

       inner: do
          ! --------------------------------------------------------------
          ! STEP 2: Step computation
          ! --------------------------------------------------------------
          s = 0.0_dp
          if ( param%printlevel .gt. 1 ) then
             write (*,fmt='(/,A)') "Calling Algencan with:"
             write (*,fmt='(2X,A,1P,D11.3)') "sigma                = ", sigma
             write (*,fmt='(2X,I0,A,1P,D16.8)',advance='no') &
                  min( n, 4 ), " coordinate(s) of x = ", x(1)
             do cnt = 2, min( n, 4 )
                write (*,fmt='(1P,2X,D16.8)',advance='no') x(cnt)
             end do
             write (*,*)
             write (*,fmt='(2X,I0,A,1P,D16.8)',advance='no') &
                  min( n, 4 ), " coordinate(s) of s = ", s(1)
             do cnt = 2, min( n, 4 )
                write (*,fmt='(1P,2X,D16.8)',advance='no') s(cnt)
             end do
             write (*,*)
          end if

          call algencan( model_evalf, model_evalg, model_evalh,         &
               model_evalc, model_evaljac, model_evalhc, model_evalfc,  &
               model_evalgjac, model_evalgjacp, model_evalhl,           &
               model_evalhlp, jcnnzmax, hnnzmax, epsfeas, epsopt,       &
               efstain, eostain, efacc, eoacc, outputfnm, specfnm,      &
               nvparam, vparam, param%theta, param%p, n, s, l, u, cons, &
               lambda, equatn, linear, coded, checkder, mval, cnorm,    &
               snorm, nlpsupn, geninfo, inform )
       
          if ( inform .ne. 0 ) then
             flag = 10 * inform
             exit outer
          end if

          if ( param % printlevel .gt. 1 ) then
             write (*,fmt='(A)') "Algencan returned with:"
             write (*,fmt='(2X,A,1P,D24.16)') "mval                 = ", mval
             write (*,fmt='(2X,I0,A,1P,D16.8)',advance='no') &
                  min( n, 4 ), " coordinate(s) of s = ", s(1)
             do cnt = 2, min( n, 4 )
                write (*,fmt='(1P,2X,D16.8)',advance='no') s(cnt)
             end do
             write (*,*)
          end if

          call model_evalg( n, s, mg, flag )
          mgnorm2 = norm2( mg )
          snorm2  = norm2( s )
          ! --------------------------------------------------------------

          ! --------------------------------------------------------------
          ! STEP 3: Step control
          ! --------------------------------------------------------------
          sinfnorm = maxval( abs( s ) )
          xinfnorm = maxval( abs( x ) )
          
          if ( j .lt. param%bigJ .and. (                                   &
               ( f - mval ) / max( 1.0_dp, abs( f ) ) .gt. param%eta1 .or. &
                   sinfnorm / max( 1.0_dp, xinfnorm ) .gt. param%eta2 )    &
             ) then

             if ( param % printlevel .gt. 1 ) then
                write (*,fmt='(/,A,A,1P,2(/,2X,A,D24.16,2X,D24.16))')    &
                     "Increasing sigma because the model minimizer ",    &
                     "does not satisfies the step control conditions",   &
                     "         ( f - mval ) / max( 1.0_dp, |f| ) = ",    &
                     ( f - mval ) / max( 1.0_dp, abs( f ) ), param%eta1, &
                     "              ||s|| / max( 1.0_dp, ||x|| ) = ",    &
                     sinfnorm / max( 1.0_dp, xinfnorm ), param%eta2
             end if

             ! Increase regularization parameter
             if ( param%sigmaini ) then
                sigma = max( sigmaini,       param % gamma2 * sigma )
             else
                sigma = max( param%sigmalow, param % gamma2 * sigma )
             end if
             
             if ( param % printlevel .gt. 1 ) then
                write (*,fmt='(/,A,1P,D9.3,/)') "sigma = ", sigma
             end if

             j = j + 1
             cycle inner

          else
             if ( j .ge. param%bigJ ) then
                if ( param%printlevel .gt. 1 ) then
                   write (*,fmt='(/,A,1X,I0,1X,A)')         &
                        "Step control exceeded maximum of", &
                        param%bigJ, "iterations."
                end if
             else
                if ( param % printlevel .gt. 1 ) then
                   write (*,fmt='(/,A,1P,2(/,2X,A,D24.16,2X,D24.16))')      &
                        "Step control conditions were satisfied",           &
                        "         ( f - mval ) / max( 1.0_dp, |f| ) = ",    &
                        ( f - mval ) / max( 1.0_dp, abs( f ) ), param%eta1, &
                        "              ||s|| / max( 1.0_dp, ||x|| ) = ",    &
                        sinfnorm / max( 1.0_dp, xinfnorm ), param%eta2
                end if
             end if
          end if
          ! --------------------------------------------------------------

          ! --------------------------------------------------------------
          ! STEP 4: Step acceptance
          ! --------------------------------------------------------------
          xtrial = x + s

          call param%uevalf( xtrial, ftrial, uflag )
          fevalcnt = fevalcnt + 1
          if ( uflag .ne. 0 ) then
             output%info = - 90
             return
          end if

          ! Print iteration information
          if ( param % printlevel .gt. 0 ) then
             if ( param % printlevel .gt. 1 .or. ( modulo( k, 15 ) .eq. 0 .and. k .ne. 0 ) ) then
                if ( param % printlevel .gt. 1 ) then
                   write (*,fmt='(/,1P,A,D24.16,/,A,I0,/)') "ftrial = ", ftrial, "fevalcnt = ", fevalcnt

                   write (*,fmt='(I0,A,/,2X,1P,D16.8)',advance='no') &
                        min( n, 4 ), " coordinate(s) of xtrial:", xtrial(1)
                   do cnt = 2, min( n, 4 )
                      write (*,fmt='(1P,2X,D16.8)',advance='no') xtrial(cnt)
                   end do
                   write (*,fmt='(/)')
                end if
             
                write (*,100) 'iter', 'sigma', 'f', '||g||', 'geninfo', '||mg||', &
                     '||s||', 'ftrial', 'dec?'
             end if
          
             write (*,200) k, sigma, f, ginfnorm, geninfo, mgnorm2, &
                  snorm2, ftrial, f-ftrial

             if ( param % printlevel .gt. 1 ) then
                write (*,fmt="(2X,136('='))")

                write (*,700) "Checking if iteration is sucessful or not...",                     &
                     "||mg|| = ", mgnorm2,                                                        &
                     "theta * ||s||^", param%p, " = ", param%theta * snorm2**param%p,             &
                     "ftrial = ", ftrial,                                                         &
                     "f - alpha*||s||^", param%p+1, " = ", f - param%alpha * snorm2**(param%p+1), &
                     "f - 1.0e+04 * epsilon^(", param%p+1, "/", param%p, ") = ",                  &
                     f - 1.0e+04_dp * param%epsilon**((param%p+1.0_dp)/param%p),                  &
                     "f - 1.0e-04 * epsilon^(", param%p+1, "/", param%p, ") = ",                  &
                     f - 1.0e-04_dp * param%epsilon**((param%p+1.0_dp)/param%p)
             end if
          end if

          ! CHECKS IF ITERATION WAS SUCESSFUL IN TERMS OF DESCENCE CONDITIONS

          ! CASE 1 : GENCAN found a point with the model gradient norm
          !          sufficiently small
          if ( mgnorm2 .le. param%theta * snorm2**param%p ) then

             if ( param%printlevel .gt. 1 ) then
                write (*,fmt='(/,A)') "Gencan converged to an optimal solution"
             end if
          
             ! Checks if ftrial meets p+1 descense or a fixed descence
             ! of order epsilon^((p+1)/p)
             if ( ftrial .le. f - param%alpha * snorm2**(param%p+1) .or. &
                  ftrial .le. f - 1.0e+04_dp * param%epsilon**((param%p+1.0_dp)/param%p) ) then
                
                successful = .true.

                if ( param%printlevel .gt. 1 ) then
                   if ( ftrial .le. f - param%alpha * snorm2**(param%p+1) ) then
                      if ( param%p .eq. 2 ) then
                         write (*,*) "   ftrial meets cubic descence"
                      else if ( param%p .eq. 3 ) then
                         write (*,*) "   ftrial meets quartic descence"
                      end if
                   end if

                   if ( ftrial .le. f - 1.0e+04_dp * param%epsilon**((param%p+1.0_dp)/param%p) ) then
                      write (*,*) "   ftrial meets fixed descence"
                   end if
                end if
             
             else
                successful = .false.

                if ( param%printlevel .gt. 1 ) then
                   write (*,*) "   ftrial does not meet any descence condition"
                end if
             end if

          ! CASE 2 : GENCAN did not converge to a point with the model
          !          gradient norm sufficiently small
          else

             if ( param%printlevel .gt. 1 ) then
                write (*,fmt='(/,A)') "Gencan did not converge to an optimal solution"
             end if

             ! Checks if ftrial meets a fixed descence of order epsilon^((p+1)/p)
             if ( ( f - ftrial ) .ge. 1.0e-04_dp * param%epsilon**((param%p+1.0_dp)/param%p) ) then
                successful = .true.

                if ( param%printlevel .gt. 1 ) then
                   write (*,*) "   ftrial meets fixed descence"
                end if
             
             else
                successful = .false.

                if ( param%printlevel .gt. 1 ) then
                   write (*,*) "   ftrial does not meet fixed descence"
                end if

                ! If the required epsilon to gencan is too low, checks if ftrial meets a simple
                ! descence condition and estabilishes a maximum number of iterations for it
                if ( param%theta * snorm2**param%p .le. 1.0e-16 .and. ftrial .le. f ) then
                   if ( param%printlevel .gt. 1 ) then
                      write (*,fmt='(4X,A,1P,D16.8)')                                    &
                           "the required optimality tolerance to gencan was so tight: ", &
                           param%theta * snorm2**param%p
                      write (*,*) "   ftrial meets simple descence"
                   end if

                   ! >>>>>>> FOR GENERALIZED MINIMIZATION PROFILE
                   ! Print the objective function value and fevalcnt
                   if ( param%trackunit .gt. 0 ) then
                      write (param%trackunit, fmt='(F40.16, 2X, I8)') ftrial, fevalcnt
                   end if

                   call param%uevalg( x, g, uflag )
                   gevalcnt = gevalcnt + 1
                   if ( uflag .ne. 0 ) then
                      output%info = - 91
                      return
                   end if
                   ginfnorm = maxval( abs( g ) )

                   ! Check stop criteria for this case.
                   if ( ginfnorm .le. param%epsilon ) then
                      if ( param%printlevel .gt. 1 ) then
                         write (*,fmt='(4X,A,1P,D24.16,/)') &
                              "the failure accepted point is an optimal solution ... ||g|| = ", ginfnorm
                      end if

                      flag = 0
                      exit outer
                   end if

                   if ( param%printlevel .gt. 1 ) then
                      write (*,fmt='(4X,A,1P,D24.16,/)') &
                           "the failure accepted point is not an optimal solution ... ||g|| = ", ginfnorm
                   end if
                
                   flag = - 3
                   exit outer
                end if
             end if
          end if

          ! Unsuccessful Iteration
          if ( .not. successful ) then
             if ( param % printlevel .gt. 1 ) then
                write (*,fmt='(/,A,/,4X,1P,A,D24.16)',advance='no') "UNSUCCESSFUL ITERATION", &
                     "Increasing sigma from ", sigma
             end if

             ! Increase regularization parameter
             if ( param%sigmaini ) then
                sigma = max( sigmaini,       param % gamma2 * sigma )
             else
                sigma = max( param%sigmalow, param % gamma2 * sigma )
             end if

             if ( param%printlevel .gt. 1 ) then
                write (*,fmt='(A,D24.16)') " to ", sigma
             end if

             j = j + 1
             cycle inner
             
          ! Successful iteration
          else
             x = xtrial
             f = ftrial
          
             if ( sigma .eq. 0.0_dp ) then
                nwstep = nwstep + 1
             end if

             ! >>>>>>> FOR GENERALIZED MINIMIZATION PROFILE
             ! Print the objective function value and fevalcnt
             if ( param%trackunit .gt. 0 ) then
                write (param%trackunit, fmt='(F40.16, 2X, I8)') f, fevalcnt
             end if

             ! Compute the derivatives of the objective in the new point
             call param%uevalg( x, g, uflag )
             gevalcnt = gevalcnt + 1
             if ( uflag .ne. 0 ) then
                output%info = - 91
                return
             end if
             ginfnorm = maxval( abs( g ) )

             if ( param%p .ge. 2 ) then
                if ( param%dense ) then
                   call param%uevalh ( x, h,  uflag )
                else
                   call param%uevalhs( x, hs, uflag )
                end if

                hevalcnt = hevalcnt + 1

                if ( uflag .ne. 0 ) then
                   output%info = - 92
                   return
                end if
             end if

             if ( param%p .ge. 3 ) then
                if ( param%dense ) then
                   call param%uevalt ( x, t,  uflag )
                else
                   call param%uevalts( x, ts, uflag )
                end if

                tevalcnt = tevalcnt + 1

                if ( uflag .ne. 0 ) then
                   output%info = - 93
                   return
                end if
             end if

             exit inner
          end if
          ! --------------------------------------------------------------

       end do inner

       ! Updates sigmaini, if it is in use
       if ( param%sigmaini ) then
          if ( param%printlevel .gt. 1 ) then
             write (*,fmt='(/,A,D16.8)',advance='no') &
                  "Updating sigmaini from ", sigmaini
          end if

          if ( sigma .ne. 0.0_dp ) then
             sigmaini = param%gamma1 * sigma
          else
             sigmaini = param%gamma1 * sigmaini
          end if

          if ( sigmaini .le. 1.0e-16_dp ) then
             sigmaini = 1.0e-16_dp
          end if

          if ( param%printlevel .gt. 1 ) then
             write (*,fmt='(/,A,D16.8)') " to ", sigmaini
          end if
       end if

       k = k + 1
    end do outer

    call cpu_time( tfinish )

    ! Report results
    if ( param % printlevel .gt. 0 ) then
       select case ( flag )
       case (-3)
          write (*,497)
          
       case (-2)
          write (*,498)

       case (0)
          write (*,500)

       case (1)
          write (*,550)

       case default
          write (*,600) flag
       end select
    end if

    ! print *, count_rate
    if ( param % printlevel .gt. 0 ) then
       write (*,800) f, ginfnorm, k, fevalcnt, gevalcnt, hevalcnt, &
            tevalcnt, ( tfinish-tstart )
    end if

    ! Set output data
    output%fcnt   = fevalcnt
    output%fval   = f
    output%gcnt   = gevalcnt
    output%gnorm  = ginfnorm
    output%hcnt   = hevalcnt
    output%info   = flag
    output%itcnt  = k
    output%nwstep = nwstep
    output%sigma  = sigma
    output%tcnt   = tevalcnt
    output%time   = tfinish - tstart

    ! Deallocate arrays
    deallocate( g, l, mg, u, s, xtrial, stat=allocstat )

    if ( param%dense ) then
       deallocate( h, t, stat=allocstat )
    else
       deallocate( hs%row, hs%col, hs%val, ts%id1, ts%id2, ts%id3, ts%val, &
            stat=allocstat )
    end if

    ! NONEXECUTABLE STATEMENTS
010 format( /, 2X, 1P, 50('*'), /, 2X,  &
         "Starting Adaptive Regularization with p = ", I0, &
         /, 2X, 50('*'), /, 2X,         &
         "alpha       = ", E9.2, /, 2X, &
         "epsilon     = ", E9.2, /, 2X, &
         "eta1        = ", E9.2, /, 2X, &
         "eta2        = ", E9.2, /, 2X, &
         "gamma1      = ", E9.2, /, 2X, &
         "gamma2      = ", E9.2, /, 2X, &
         "hnnzmax     = ", I9,   /, 2X, &
         "J           = ", I9,   /, 2X, &
         "kmax        = ", I9,   /, 2X, &
         "printlevel  = ", I9,   /, 2X, &
         "sigmaini    = ", L1,   /, 2X, &
         "sigmalow    = ", E9.2, /, 2X, &
         "theta       = ", E9.2, /, 2X, &
         "tnnzmax     = ", I9,   /, 2X, &
         "trackunit   = ", I9,          &
         /, 2X, 50('*'), / )

100 format( 2X, 136('='), /, 2X, A7, 2X, A10, 2(2X, A16), 2X, A7, 4(2X, A16), /, 2X, 136('=') )
110 format( 3(2X, A16), 2X, A4 )
200 format( 2X, I7, 1P, 2X, D10.3, 2(2X, D16.8), 2X, I7, 4(2X, D16.8), 1X, I0 )

497 format( /, 2X, 'EXIT: Gencan failed to solve subproblem and no descence holds in the current point.', / )
498 format( /, 2X, 'EXIT: Maximum iterations exceeded.', / )
500 format( /, 2X, 'EXIT: Optimal solution found (absolute norm).', / )
550 format( /, 2X, 'EXIT: Optimal solution found (relative norm).', / )
600 format( /, 2X, 'EXIT: ARP finished with flag ', I0, '.', / )

700 format( /, A, /, 24X, A, 1P, D24.16, /, 15X, A, I0, A, D24.16, /, 24X, A, D24.16, &
            /, 17X, A, I0, A, D24.16, 2(/, 5X, 2(A, I0), A, D24.16))

800 format(/, 1X, "Objective Function ...............: ", 1P, D24.16, &
           /, 1X, "Gradient Norm ....................: ", 1P, D24.16, &
           /, 1X, "Number of Iterations .............: ", I8,         &
           /, 1X, "Number of Function Evaluations ...: ", I8,         &
           /, 1X, "Number of Gradient Evaluations ...: ", I8,         &
           /, 1X, "Number of Hessian Evaluations ....: ", I8,         &
           /, 1X, "Number of Tensor Evaluations .....: ", I8,         &
          //, 1X, "Total CPU Time (in seconds): ", 0PF12.2 )

  end subroutine arp_solve

  ! ==================================================================
  ! ==================================================================

end module arp
