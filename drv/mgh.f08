program mgh

  use arp
  use mgh

  implicit none

  integer, parameter :: dp = kind(0.0d0)

  ! LOCAL SCALARS
  character(len=30) :: arg, filename
  character(len=60) :: probname
  integer           :: allocstat, bigJ, cnt, firstprob, flag, i, iostat, &
                       lastprob, m, n, p, pl, prob
  logical           :: sigmaini, rmp

  ! ARP PARAMETERS
  type(arp_param) :: param

  ! ARP OUTPUT
  type(arp_output) :: output

  ! LOCAL ARRAYS
  real(kind=dp), allocatable, dimension(:) :: x

  bigJ      = 20
  n         = 0
  m         = 0
  p         = 0
  pl        = 1
  firstprob = 0
  lastprob  = 0
  sigmaini  = .true.
  rmp       = .false.

  i = 1
  do while ( i .le. command_argument_count() )
     call get_command_argument( i, arg )

     select case ( trim( arg ) )

     case ( '-p' )
        i = i + 1
        call get_command_argument( i, arg )
        read( arg, *, iostat=iostat ) p
        if ( iostat .ne. 0 ) p = 3

     case ( '-pl' )
        i = i + 1
        call get_command_argument( i, arg )
        read( arg, *, iostat=iostat ) pl
        if ( iostat .ne. 0 ) pl = 1

     case ( '-prob' )
        i = i + 1
        call get_command_argument( i, arg )
        cnt = scan( arg, ":")
        if ( cnt .gt. 0 ) then
           read( arg(1:cnt-1), *, iostat=iostat ) firstprob
           if ( iostat .ne. 0 ) firstprob = 0

           read( arg(cnt+1:),  *, iostat=iostat ) lastprob
           if ( iostat .ne. 0 ) lastprob = 0
        else
           read( arg, *, iostat=iostat ) firstprob
           if ( iostat .ne. 0 ) firstprob = 0
           lastprob = firstprob
        end if

     case ( '-m' )
        i = i + 1
        call get_command_argument( i, arg )
        read( arg, *, iostat=iostat ) m
        if ( iostat .ne. 0 ) m = 0
        
     case ( '-n' )
        i = i + 1
        call get_command_argument( i, arg )
        read( arg, *, iostat=iostat ) n
        if ( iostat .ne. 0 ) n = 0

     case ( '-J' )
        i = i + 1
        call get_command_argument( i, arg )
        read( arg, *, iostat=iostat ) bigJ
        if ( iostat .ne. 0 ) bigJ = 20

     case ( '-alg21' )
        sigmaini = .false.

     case ( '-track' )
        rmp = .true.
        
     end select

     i = i + 1
  end do

  if ( p .eq. 0 .or. firstprob .eq. 0 ) then
     write ( *, * ) "Command line arguments:"
     write ( *, * ) "    -p      avaliable derivatives ........ {1,2,3},  mandatory"
     write ( *, * ) "    -prob   problem number or interval ... 1:35,     mandatory"
     write ( *, * ) "    -n      dimension .................... >0,       optional, default: problem default"
     write ( *, * ) "    -m      problem equations ............ >0,       optional, default: problem default"
     write ( *, * ) "    -pl     print level .................. 0..2,     optional, default: 1"
     write ( *, * ) "    -J      J value ... .................. >=0,      optional, default: 20"
     write ( *, * ) "    -alg21  runs alg. 2.1 ................ no param, optional, default: runs algorithm 4.1"
     write ( *, * ) "    -track  per iteration information .... no param, optional"
     
     return
  end if
  
  do prob = firstprob, lastprob
     ! Set MGH problem
     call mgh_set_problem( prob, flag )

     ! Set problem dimensions
     if ( n .eq. 0 .and. m .eq. 0 ) then
        call mgh_get_dims( n, m )
     else
        call mgh_set_dims( n, m, flag )
        if ( flag .ne. 0 ) then
           select case ( flag )
           case ( -1 ) ! n is not valid
              call mgh_get_dims( n=n )

           case ( -2 ) ! m is not valid
              call mgh_get_dims( m=m )

           case ( -3 ) ! both n and m are not valid
              write ( *, * ) "Both m and n are not valid. Please choose valid ones."
              stop
           end select
        end if
     end if

     ! Allocate arrays
     allocate( x(n), stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write ( *, * ) "Allocation error in run_mgh."
        return
     end if

     ! Get initial point
     call mgh_get_x0( x )

     ! Initialize default arp parameters
     call arp_load( param )

     ! Set parameters
     param%p          = p
     param%printlevel = pl
     param%uevalf => mgh_evalf
     param%uevalg => mgh_evalg
     param%uevalh => mgh_evalh
     param%uevalt => mgh_evalt

     param%bigJ     = bigJ
     param%sigmaini = sigmaini

     if ( rmp ) then
        write ( filename, '(I0,A)' ) prob, "_rmp.txt"
        open ( 55, file=filename )
        param%trackunit = 55
     else
        param%trackunit = 0
     end if

     ! Optimization
     call arp_solve( x, param, output )

     ! Save results
     write ( filename, '(A,I0)' ) 'results_tab_p', p
     open ( 25, file=filename, position='append' )

     call mgh_get_name( probname )
     write ( 25, 100 ) prob, trim( probname ), n, m, output%info,     &
          output%fval, output%gnorm, output%itcnt, output%fcnt, &
          output%gcnt, output%time, output%sigma, output%nwstep

     close ( 25 )
     if ( rmp ) close ( 55 )

     ! Deallocate array
     deallocate( x, stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write ( *, * ) "Deallocation error in mgh driver."
        stop
     end if
  end do

  ! NONEXECUTABLE STATEMENTS
100 format( 2X, I2, 2X, A42, 2(2X, I4), 1X, 1P, '|  ', I3, 1X, &
            2(D16.8, 1X), 3(I8, 1X), F10.2, 1X, D10.3, 1X, I5, &
            1X, I5, ' |' )
  
end program mgh
