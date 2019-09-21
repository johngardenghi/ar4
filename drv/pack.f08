program pack

  use arp
  
  implicit none

  ! ARP PARAMS
  type(arp_param) :: param

  ! ROUTINE OUTPUT
  type(arp_output) :: output

  ! LOCAL SCALARS
  character(len=3)  :: draw
  character(len=30) :: arg, inivmd, solvmd_arp, solvmd_alg
  character(len=80) :: specfnm, outputfnm
  integer           :: allocstat, geninfo, hnnzmax, i, ind, inform,   &
                       iostat, j, jcnnzmax, m, n, ndim, nite, nnereg, &
                       nreg, nvparam, p, pl
  logical           :: checkder, internal, run_arp, run_algencan
  real(kind=8)      :: cnorm, drand, epsa, efacc, efstain, eoacc,  &
                       eostain, epsr, epsfeas, epsopt, f, nlpsupn, &
                       objrad, r, seed, snorm, theta

  ! LOCAL ARRAYS
  character(len=80),              dimension(10)    :: vparam
  integer,           allocatable, dimension(:)     :: next
  integer,           allocatable, dimension(:,:)   :: nonempty
  integer,           allocatable, dimension(:,:,:) :: start
  logical,                        dimension(0)     :: equatn, linear
  logical,                        dimension(11)    :: coded
  real(kind=8),                   dimension(0)     :: lambda
  real(kind=8),      allocatable, dimension(:)     :: l, u, x, xini

  ! LOCAL SPARSE MATRIX
  type(sparse_matrix) :: hs

  ! EXTERNAL ROUTINES
  external :: drand

  ! Iniitialize standard values
  ndim = 3
  r    = 1.0d0

  ! Expected to get as command line parameter
  nite         = 0
  objrad       = 0.0d0
  draw         = 'pov'
  p            = 3
  pl           = 1
  internal     = .true.
  run_arp      = .true.
  run_algencan = .true.
  
  do i = 1, command_argument_count(), 2
     call get_command_argument( i, arg )

     select case ( trim( arg ) )

     case ( '-N' )
        call get_command_argument( i+1, arg )
        read( arg, *, iostat=iostat ) nite
        if ( iostat .ne. 0 ) nite = 0

     case ( '-p' )
        call get_command_argument( i+1, arg )
        read( arg, *, iostat=iostat ) p
        if ( iostat .ne. 0 ) p = 3
        
     case ( '-r' )
        call get_command_argument( i+1, arg )
        read( arg, *, iostat=iostat ) objrad
        if ( iostat .ne. 0 ) objrad = 0.0d0

     case ( '-d' )
        call get_command_argument( i+1, arg )
        read( arg, *, iostat=iostat ) draw
        if ( iostat .ne. 0 ) draw = 'pov'

     case ( '-i' )
        call get_command_argument( i+1, arg )
        read( arg, *, iostat=iostat ) internal
        if ( iostat .ne. 0 ) internal = .true.

     case ( '-pl' )
        call get_command_argument( i+1, arg )
        read( arg, *, iostat=iostat ) pl
        if ( iostat .ne. 0 ) pl = 1

     case ( '-s' )
        call get_command_argument( i+1, arg )

        if ( trim( arg ) .eq. 'alg' ) then
           run_arp = .false.
        elseif ( trim( arg ) .eq. 'arp' ) then
           run_algencan = .false.
        end if

     end select
  end do

  if ( nite .eq. 0 ) then
     write ( *, * ) "Command line arguments:"
     write ( *, * ) "    -N      number of items ........... >0,        mandatory"
     write ( *, * ) "    -p      available derivatives ..... {2,3},     optional, default: 3"
     write ( *, * ) "    -r      larger sphere radius ...... >0,        optional, default: density 0.3"
     write ( *, * ) "    -d      solution drawing method ... {pov,vmd}, optional, default: pov"
     write ( *, * ) "    -i      initial point internal .... {T,F},     optional, default: T"
     write ( *, * ) "    -s      solver to run ............. {alg,arp}, optional, default: both"
     write ( *, * ) "    -pl     print level ............... 0..2,      optional, default: 1"

     return
  end if

  if ( objrad .eq. 0.0d0 ) then
     objrad = r *  ( nite / 3.0d-01 ) ** ( 1.0d0 / 3.0d0 )
  end if

  n    = ndim * nite
  nreg = ceiling( objrad / r )

  allocate( l(n), u(n), xini(n), x(n), start(0:nreg+1,0:nreg+1,0:nreg+1), &
       next(nite), nonempty(3,nite), stat=allocstat )
  if ( allocstat .ne. 0 ) then
     write (*,*) 'Allocation error in pack'
     stop
  end if

  start(0:nreg+1,0:nreg+1,0:nreg+1) = 0
  nnereg = 0

  ! Randomly generates initial point
  seed = 123456.0d0

  i = 1
  do while ( i .le. nite )
     ind = ndim * ( i - 1 )

     do j = ind+1, ind+ndim
        xini(j) = - objrad + 2.0d0 * objrad * drand( seed )
     end do

     if ( .not. internal .or. ( sum( xini(ind+1:ind+ndim) ** 2 ) - ( objrad - r ) ** 2 ) .le. 0.0d0 ) then
        i = i + 1
     end if
  end do

  write (*,1) "-----------------------------------------"
  write (*,1) "             PACKING PROBLEM             "
  write (*,1) "-----------------------------------------"
  write (*,2) "Dimension ........:", 3
  write (*,2) "Number of items ..:", nite
  write (*,3) "Items radius .....:", 1.0
  write (*,3) "Object radius ....:", objrad
  write (*,4) "Feasible initial? :", internal
  write (*,5) "Draw method ......:", draw
  write (*,1) "-----------------------------------------"
     
  write (inivmd, fmt='(I0,A,A)') nite, '_ini.', draw
  write (solvmd_arp, fmt='(I0,A,I0,A,A)') nite, '_arp_p', p, '.', draw
  write (solvmd_alg, fmt='(I0,A,A)') nite, '_algencan.', draw

  open (15, file=inivmd)
  if ( draw .eq. 'vmd' ) then
     call drawsol( n, xini, 15 )
  else
     call drawsol2( n, xini, 15 )
  end if
  close (15)

  if ( run_arp ) then
     call arp_load( param )

     if ( n .gt. ( huge( 1 ) ** 0.5 ) / 3.0d0 ) then
        param%hnnzmax = 2000000000
     else
        param%hnnzmax = 3 * n ** 2
     end if

     if ( n .gt. ( huge( 1 ) ** ( 1.0d0 / 3.0d0 ) ) / 3.0d0 ) then
        param%tnnzmax = 2000000000
     else
        param%tnnzmax = 3 * n ** 3
     end if

     param%p = p

     param%uevalf  => packing_evalf_arp
     param%uevalg  => packing_evalg_arp
     if ( p .eq. 2 ) then
        param%uevalh => packing_evalh_arp
        param%dense = .true.

        allocate( hs%row(param%hnnzmax), hs%col(param%hnnzmax), &
             hs%val(param%hnnzmax), stat=allocstat )
        if ( allocstat .ne. 0 ) then
           write (*,*) 'Allocation error in pack'
           stop
        end if

     elseif ( p .eq. 3 ) then
        param%uevalhs => packing_evalhs_arp
        param%uevalts => packing_evalts_arp
        param%dense = .false.
     else
        write ( *, * ) "For ARP, p must be 2 or 3."
        stop
     end if

     param%printlevel = pl

     x(1:n) = xini(1:n)
     call arp_solve( x(1:n), param, output )

     write ( *, * ) "ARP output info: ", output%info

     open (20, file=solvmd_arp)
     if ( draw .eq. 'vmd' ) then
        call drawsol( n, x, 20 )
     else
        call drawsol2( n, x, 20 )
     end if

     if ( p .eq. 2 ) then
        deallocate( hs%row, hs%col, hs%val, stat=allocstat )
     end if
     
     close (20)
  end if
  
  if ( run_algencan ) then
     m = 0

     l(1:n) = - huge( 0.0d0 )
     u(1:n) =   huge( 0.0d0 )

     coded(1:3)  = .true.  ! f, g, h
     coded(4:11) = .false. ! c, jac, hc, fc, gjac, gjacp, hl, hlp

     ! Upper bounds on the number of sparse-matrices non-null elements
     jcnnzmax = 0

     if ( n .gt. ( huge( 1 ) ** 0.5 ) ) then
        hnnzmax = 2000000000
     else
        hnnzmax  = n * ( n + 1 )
     end if

     ! Checking derivatives?
     checkder = .false.

     ! Parameters setting
     epsfeas   = 1.0d-08
     epsopt    = 1.0d-08

     efstain   = sqrt( epsfeas )
     eostain   = epsopt ** 1.5d0

     efacc     = sqrt( epsfeas )
     eoacc     = sqrt( epsopt )

     theta  = 1.0d+20

     write (outputfnm,'(I0,A)') n/3, '_algencan_nw.out'
     specfnm   = ''

     ! Optimize
     call system( 'rm .silent' )

     nvparam   = 2
     vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 11'
     vparam(2) = 'NEWTON-LINE-SEARCH-INNER-SOLVER'

     ! Run with Newton Method as inner solver
     x(1:n) = xini(1:n)
     call algencan(packing_evalf, packing_evalg, packing_evalhs,       &
          packing_evalc, packing_evaljac, packing_evalhc,              &
          packing_evalfc, packing_evalgjac, packing_evalgjacp,         &
          packing_evalhl, packing_evalhlp, jcnnzmax, hnnzmax, epsfeas, &
          epsopt, efstain, eostain, efacc, eoacc, outputfnm, specfnm,  &
          nvparam, vparam, theta, 0, n, x, l, u, m, lambda,            &
          equatn, linear, coded, checkder, f, cnorm, snorm, nlpsupn,   &
          geninfo, inform)

     ! Draw solution
     open (20, file=solvmd_alg)
     if ( draw .eq. 'vmd' ) then
        call drawsol( n, x, 20 )
     else
        call drawsol2( n, x, 20 )
     end if
     close (20)

     ! Run with default option for inner solver
     nvparam = 1
     write (outputfnm,'(I0,A)') n/3, '_algencan_df.out'
     x(1:n) = xini(1:n)
     call algencan(packing_evalf, packing_evalg, packing_evalhs,       &
          packing_evalc, packing_evaljac, packing_evalhc,              &
          packing_evalfc, packing_evalgjac, packing_evalgjacp,         &
          packing_evalhl, packing_evalhlp, jcnnzmax, hnnzmax, epsfeas, &
          epsopt, efstain, eostain, efacc, eoacc, outputfnm, specfnm,  &
          nvparam, vparam, theta, 0, n, x, l, u, m, lambda,            &
          equatn, linear, coded, checkder, f, cnorm, snorm, nlpsupn,   &
          geninfo, inform)

     call system( 'touch .silent' )
  end if

  deallocate( l, u, xini, x, start, next, nonempty, stat=allocstat )

  ! NON-EXECUTABLE STATEMENTS
1 format( 1X, A )
2 format( 1X, A, I20 )
3 format( 1X, A, F20.6 )
4 format( 1X, A, 19X, L )
5 format( 1X, A, A20 )

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine packing_evalf_arp(x,f,flag)

    ! SCALAR ARGUMENT
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENT
    real(kind=8), dimension(:), intent(in) :: x

    call packing_evalf( size( x ), x, f, flag )
    
  end subroutine packing_evalf_arp
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine packing_evalf(n,x,f,flag)
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: i,ind1,j,k,l,p

    ! Compute objective function

    flag = 0

    call classify(n,x)

    f = 0.0d0
    do l = 1,nnereg
       i = nonempty(1,l)
       j = nonempty(2,l)
       k = nonempty(3,l)
       p = start(i,j,k)
       do while ( p .ne. 0 )
          f = f + nrdist(n,x,p,next(p))            &
                + nrdist(n,x,p,start(i+1,j-1,k  )) &
                + nrdist(n,x,p,start(i+1,j,  k  )) &
                + nrdist(n,x,p,start(i+1,j+1,k  )) &
                + nrdist(n,x,p,start(i,  j+1,k  )) &
                + nrdist(n,x,p,start(i+1,j-1,k-1)) &
                + nrdist(n,x,p,start(i+1,j,  k-1)) &
                + nrdist(n,x,p,start(i+1,j+1,k-1)) &
                + nrdist(n,x,p,start(i,  j+1,k-1)) &
                + nrdist(n,x,p,start(i+1,j-1,k+1)) &
                + nrdist(n,x,p,start(i+1,j,  k+1)) &
                + nrdist(n,x,p,start(i+1,j+1,k+1)) &
                + nrdist(n,x,p,start(i,  j+1,k+1)) &
                + nrdist(n,x,p,start(i,  j,  k+1))
          p = next(p)
       end do
    end do

    do i = 1, nite
       ind1 = ( i - 1 ) * ndim
       f = f + max( 0.0d0, sum( x(ind1+1:ind1+ndim) ** 2 ) - ( objrad - r ) ** 2 ) ** 4
    end do

  end subroutine packing_evalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalg_arp(x,g,flag)
    ! SCALAR ARGUMENTS
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(:), intent(in)  :: x
    real(kind=8), dimension(:), intent(out) :: g

    call packing_evalg( size( x ), x, g, flag )

  end subroutine packing_evalg_arp

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalg(n,x,g,flag)
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    integer      :: i,ind1,j,k,l,p
    real(kind=8) :: fparc

    ! LOCAL ARRAYS
    real(kind=8) :: val(ndim),xdiff(ndim)

    ! Compute gradient of the objective function

    flag = 0

    call classify(n,x)

    g(1:n) = 0.0d0

    do l = 1,nnereg
       i = nonempty(1,l)
       j = nonempty(2,l)
       k = nonempty(3,l)
       p = start(i,j,k)
       do while ( p .ne. 0 )
          call gnrdist(n,x,g,p,next(p))           
          call gnrdist(n,x,g,p,start(i+1,j-1,k  ))
          call gnrdist(n,x,g,p,start(i+1,j,  k  ))
          call gnrdist(n,x,g,p,start(i+1,j+1,k  ))
          call gnrdist(n,x,g,p,start(i,  j+1,k  ))
          call gnrdist(n,x,g,p,start(i+1,j-1,k-1))
          call gnrdist(n,x,g,p,start(i+1,j,  k-1))
          call gnrdist(n,x,g,p,start(i+1,j+1,k-1))
          call gnrdist(n,x,g,p,start(i,  j+1,k-1))
          call gnrdist(n,x,g,p,start(i+1,j-1,k+1))
          call gnrdist(n,x,g,p,start(i+1,j,  k+1))
          call gnrdist(n,x,g,p,start(i+1,j+1,k+1))
          call gnrdist(n,x,g,p,start(i,  j+1,k+1))
          call gnrdist(n,x,g,p,start(i,  j,  k+1))
          p = next(p)
       end do
    end do

    do i = 1, nite
       ind1 = ( i - 1 ) * ndim
       fparc = sum( x(ind1+1:ind1+ndim) ** 2 ) - ( objrad - r ) ** 2
       if ( fparc .gt. 0.0d0 ) then
          g(ind1+1:ind1+ndim) = g(ind1+1:ind1+ndim) + &
               8.0d0 * fparc ** 3 * x(ind1+1:ind1+ndim)
       end if
    end do
    
  end subroutine packing_evalg

  ! *****************************************************************
  ! *****************************************************************

  subroutine packing_evalh_arp(x,h,flag)
    ! SCALAR ARGUMENT
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(:),   intent(in)  :: x
    real(kind=8), dimension(:,:), intent(out) :: h

    ! This routine compute the dense upper triangle of the Hessian of
    ! the objective function

    ! LOCAL SCALARS
    integer :: i, j, k

    call packing_evalhs_arp( x, hs, flag )

    do j = 1, size( x )
       do i = 1, j
          h(i,j) = 0.0d0
       end do
    end do

    do k = 1, hs%nnz
       i = hs%col(k)
       j = hs%row(k)
       h(i,j) = hs%val(k)
    end do

  end subroutine packing_evalh_arp

  ! *****************************************************************
  ! *****************************************************************

  subroutine packing_evalhs_arp(x,hs,flag)
    ! SCALAR ARGUMENT
    integer, intent(out) :: flag

    ! ARRAY ARGUMENT
    real(kind=8), dimension(:), intent(in) :: x

    ! SPARSE MATRIX ARGUMENT
    type(sparse_matrix), intent(inout) :: hs

    ! LOCAL SCALARS
    logical :: lmem
    integer :: hnnz

    call packing_evalhs( size( x ), x, hs%row, hs%col, hs%val, &
         hs%nnz, size( hs%row ), lmem, flag )

    if ( lmem ) then
       write ( *, * ) "Lack of memory in packing_evalh"
       stop
    end if
    
  end subroutine packing_evalhs_arp

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalhs(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in)  :: lim, n
    integer, intent(out) :: flag, hnnz

    ! ARRAY ARGUMENTS
    integer,      dimension(lim), intent(out) :: hcol, hrow
    real(kind=8), dimension(n),   intent(in)  :: x
    real(kind=8), dimension(lim), intent(out) :: hval

    ! LOCAL SCALARS
    integer      :: i,ind1,j,k,l,p
    real(kind=8) :: fparc

    ! LOCAL ARRAYS
    real(kind=8) :: diagblock(nite,ndim,ndim),xdiff(ndim)

    ! TESTING
    character(len=20) :: filename

    ! Compute (lower triangle of the) Hessian of the objective function

    flag = 0
    lmem = .false.

    call classify(n,x)

    diagblock(1:nite,1:ndim,1:ndim) = 0.0d0

    hnnz = 0

    do l = 1,nnereg
       i = nonempty(1,l)
       j = nonempty(2,l)
       k = nonempty(3,l)
       p = start(i,j,k)
       do while ( p .ne. 0 )
          call hnrdist(n,x,p,next(p)           ,diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j-1,k  ),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j,  k  ),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j+1,k  ),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i,  j+1,k  ),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j-1,k-1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j,  k-1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j+1,k-1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i,  j+1,k-1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j-1,k+1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j,  k+1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i+1,j+1,k+1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i,  j+1,k+1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return                                                           
          call hnrdist(n,x,p,start(i,  j,  k+1),diagblock,hrow,hcol,hval,hnnz,lim,lmem)
          if ( lmem ) return

          p = next(p)
       end do
    end do

    do i = 1, nite
       ind1 = ( i - 1 ) * ndim
       fparc = sum( x(ind1+1:ind1+ndim) ** 2 ) - ( objrad - r ) ** 2
       if ( fparc .gt. 0.0d0 ) then
          do k = 1, ndim
             diagblock(i,k,k) = diagblock(i,k,k) + 8.0d0 * fparc ** 3 + &
                  48.0d0 * fparc ** 2 * x(ind1+k) ** 2
             do l = 1, k-1
                diagblock(i,k,l) = diagblock(i,k,l) + &
                     48.0d0 * fparc ** 2 * x(ind1+k) * x(ind1+l)
             end do
          end do
       end if
    end do

    do i = 1,nite
       do k = 1,ndim
          do l = 1,k
             if ( hnnz + 1 .gt. lim ) then
                lmem = .true.
                return
             end if

             hrow(hnnz+1) = ndim * ( i - 1 ) + k
             hcol(hnnz+1) = ndim * ( i - 1 ) + l
             hval(hnnz+1) = diagblock(i,k,l)
             hnnz = hnnz + 1
          end do
       end do
    end do
    
  end subroutine packing_evalhs

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalts_arp(x,ts,flag)
    ! SCALAR ARGUMENT
    integer, intent(out) :: flag

    ! ARRAY ARGUMENT
    real(kind=8), dimension(:), intent(in) :: x

    ! SPARSE TENSOR ARGUMENT
    type(sparse_tensor), intent(inout) :: ts

    ! LOCAL SCALARS
    logical :: lmem

    call packing_evalts( size( x ), x, ts%id1, ts%id2, ts%id3, ts%val, &
         ts%nnz, size( ts%id1 ), lmem, flag )

    if ( lmem ) then
       write ( *, * ) "Lack of memory in packing_evalt"
       stop
    end if

  end subroutine packing_evalts_arp

  ! *****************************************************************
  ! *****************************************************************

  subroutine packing_evalts(n,x,tid1,tid2,tid3,tval,tnnz,lim,lmem,flag)
    ! SCALAR ARGUMENTS
    integer, intent(in ) :: lim, n
    integer, intent(out) :: flag, tnnz
    logical, intent(out) :: lmem

    ! ARRAY ARGUMENTS
    integer,      dimension(lim), intent(out) :: tid1, tid2, tid3
    real(kind=8), dimension(n),   intent(in)  :: x
    real(kind=8), dimension(lim), intent(out) :: tval

    ! LOCAL SCALARS
    integer      :: i, ind1, j, k, l, p
    real(kind=8) :: fparc, val

    ! LOCAL ARRAYS
    real(kind=8), dimension(nite,ndim,ndim,ndim) :: diagblock
    real(kind=8), dimension(ndim)                :: xdiff

    ! TESTING
    character(len=20) :: filename

    flag = 0
    lmem = .false.

    call classify( n, x )
    
    diagblock(1:nite,1:ndim,1:ndim,1:ndim) = 0.0d0

    tnnz = 0

    do l = 1,nnereg
       i = nonempty(1,l)
       j = nonempty(2,l)
       k = nonempty(3,l)
       p = start(i,j,k)
       do while ( p .ne. 0 )
          call tnrdist(n,x,p,next(p)           ,diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j-1,k  ),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j,  k  ),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j+1,k  ),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i,  j+1,k  ),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j-1,k-1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j,  k-1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j+1,k-1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i,  j+1,k-1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j-1,k+1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j,  k+1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i+1,j+1,k+1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i,  j+1,k+1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return                                                           
          call tnrdist(n,x,p,start(i,  j,  k+1),diagblock,tid1,tid2,tid3,tval,tnnz,lim,lmem)
          if ( lmem ) return

          p = next(p)
       end do
    end do

    do i = 1, nite
       ind1 = ( i - 1 ) * ndim
       fparc = sum( x(ind1+1:ind1+ndim) ** 2 ) - ( objrad - r ) ** 2
       if ( fparc .gt. 0.0d0 ) then
          do k = 1, ndim
             diagblock(i,k,k,k) = diagblock(i,k,k,k) + &
                  48.0d0 * fparc * x(ind1+k) * ( 3.0d0 * fparc + 4.0d0 * x(ind1+k) ** 2 )

             do l = 1, ndim
                if ( l .ne. k ) then
                   val = 48.0d0 * fparc * x(ind1+l) * ( fparc + 4.0d0 * x(ind1+k) ** 2 )

                   if ( l .lt. k ) then
                      diagblock(i,k,k,l) = diagblock(i,k,k,l) + val
                   else
                      diagblock(i,l,k,k) = diagblock(i,l,k,k) + val
                   end if

                   if ( l .lt. k ) then
                      do p = 1, l-1
                         diagblock(i,k,l,p) = diagblock(i,k,l,p) + &
                              192.0d0 * fparc * x(ind1+k) * x(ind1+l) * x(ind1+p)
                      end do
                   end if
                end if
             end do
          end do
       end if
    end do

    do i = 1, nite
       do k = 1, ndim
          do l = 1, ndim
             do p = 1, ndim
                if ( abs( diagblock(i,k,l,p) ) .gt. 0.0d0 ) then
                   if ( tnnz + 1 .gt. lim ) then
                      lmem = .true.
                      return
                   end if

                   tnnz = tnnz + 1
                   tid1(tnnz) = ndim * ( i - 1 ) + k
                   tid2(tnnz) = ndim * ( i - 1 ) + l
                   tid3(tnnz) = ndim * ( i - 1 ) + p
                   tval(tnnz) = diagblock(i,k,l,p)
                end if
             end do
          end do
       end do
    end do

  end subroutine packing_evalts

  ! *****************************************************************
  ! *****************************************************************

  subroutine packing_evalc(n,x,ind,c,flag)
    ! SCALAR ARGUMENTS
    integer, intent(in) :: ind,n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: ini

    ! Compute ind-th constraint

    flag = -1

  end subroutine packing_evalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)
    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: ind,lim,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    ! LOCAL SCALARS
    integer :: i,ini

    ! Compute gradient of the ind-th constraint

    flag = -1
    lmem = .false.

  end subroutine packing_evaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)
    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: ind,lim,n
    integer, intent(out) :: flag,hcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hccol(lim),hcrow(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: hcval(lim)

    ! LOCAL SCALARS
    integer :: i,ini

    ! Compute gradient of the ind-th constraint

    flag = -1
    lmem = .false.

  end subroutine packing_evalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalfc(n,x,f,m,c,flag)
    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m)

    flag = - 1

  end subroutine packing_evalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)
    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: lim,m,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    flag = - 1

  end subroutine packing_evalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalgjacp(n,x,g,m,p,q,work,gotj,flag)
    ! SCALAR ARGUMENTS
    logical, intent(inout) :: gotj
    integer, intent(in) :: m,n
    integer, intent(out) :: flag
    character, intent(in) :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out) :: g(n)

    flag = - 1

  end subroutine packing_evalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)
    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: lim,m,n
    integer, intent(out) :: flag,hlnnz
    real(kind=8), intent(in) :: sf

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hlcol(lim),hlrow(lim)
    real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hlval(lim)

    flag = - 1

  end subroutine packing_evalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine packing_evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)
    ! SCALAR ARGUMENTS
    logical, intent(inout) :: goth
    integer, intent(in) :: m,n
    integer, intent(out) :: flag
    real(kind=8), intent(in) :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    flag = - 1

  end subroutine packing_evalhlp

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine classify(n,x)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in) :: x

    ! LOCAL SCALARS
    integer i,ind,j,k,l

    ! CLEAN-UP THE START STRUCTURE

    do l = 1,nnereg
       i = nonempty(1,l)
       j = nonempty(2,l)
       k = nonempty(3,l)
       start(i,j,k) = 0
    end do

    ! FILL-IN THE START STRUCTURE AGAIN

    nnereg = 0
    do l = 1,nite
       ind = ( l - 1 ) * ndim
       call region(x(ind+1),x(ind+2),x(ind+3),i,j,k)

       if ( start(i,j,k) .eq. 0 ) then
          nnereg = nnereg + 1
          nonempty(1,nnereg) = i
          nonempty(2,nnereg) = j
          nonempty(3,nnereg) = k
       end if

       next(l) = start(i,j,k)
       start(i,j,k) = l
    end do

  end subroutine classify

  ! ******************************************************************
  ! ******************************************************************

  subroutine region(x,y,z,i,j,k)

    ! SCALAR ARGUMENTS
    integer,      intent(out) :: i,j,k
    real(kind=8), intent(in)  :: x,y,z

    i = min( max( 1, 1 + int( (x + nreg * r) / (2.0d0 * r) ) ), nreg )
    j = min( max( 1, 1 + int( (y + nreg * r) / (2.0d0 * r) ) ), nreg )
    k = min( max( 1, 1 + int( (z + nreg * r) / (2.0d0 * r) ) ), nreg )

  end subroutine region

  ! ******************************************************************
  ! ******************************************************************

  real(kind=8) function nrdist(n,x,i,lstart)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i, lstart, n

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in) :: x

    ! LOCAL SCALARS
    integer      :: ind1, ind2, j
    real(kind=8) :: dist

    nrdist = 0.0d0
    j = lstart
    ind1 = ( i - 1 ) * ndim
    do while ( j .ne. 0 )
       ind2 = ( j - 1 ) * ndim
       dist = sum( ( x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim) ) ** 2 )
       nrdist = nrdist + max( 0.0d0, ( 2.0d0 * r ) ** 2 - dist ) ** 4
       j = next(j)
    end do

  end function nrdist

  ! ******************************************************************
  ! ******************************************************************

  subroutine gnrdist(n,x,g,i,lstart)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: i, lstart, n

    ! ARRAY ARGUMENTS
    real(kind=8), dimension(n), intent(in)    :: x
    real(kind=8), dimension(n), intent(inout) :: g

    ! LOCAL SCALARS
    integer      :: ind1, ind2, j
    real(kind=8) :: dist

    ! LOCAL ARRAYS
    real(kind=8), dimension(ndim) :: val, xdiff
  
    j = lstart

    ind1 = ( i - 1 ) * ndim

    do while ( j .ne. 0 )
       ind2 = ( j - 1 ) * ndim
       xdiff(1:ndim) = x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim)
       dist = ( 2.0d0 * r ) ** 2 - sum( xdiff(1:ndim) ** 2 )

       if ( dist .gt. 0.d0 ) then
          val(1:ndim) = 8.0d0 * dist ** 3 * xdiff(1:ndim)
          g(ind1+1:ind1+ndim) = g(ind1+1:ind1+ndim) - val(1:ndim)
          g(ind2+1:ind2+ndim) = g(ind2+1:ind2+ndim) + val(1:ndim)
       end if

       j = next(j)
    end do
    
  end subroutine gnrdist

  ! ******************************************************************
  ! ******************************************************************

  subroutine hnrdist(n,x,i,lstart,diagblock,hrow,hcol,hval,hnnz,lim,lmem)

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: i,lim,lstart,n
    logical, intent(inout) :: lmem
    integer, intent(inout) :: hnnz

    ! ARRAY ARGUMENTS
    integer,      dimension(lim),            intent(out)   :: hcol,hrow
    real(kind=8), dimension(n),              intent(in)    :: x
    real(kind=8), dimension(nite,ndim,ndim), intent(inout) :: diagblock
    real(kind=8), dimension(lim),            intent(out)   :: hval

    ! LOCAL SCALARS
    integer      :: col, ind1, ind2, j, k, l, row
    real(kind=8) :: fparc, val

    ! LOCAL ARRAYS
    real(kind=8), dimension(ndim) :: xdiff

    j = lstart
    do while ( j .ne. 0 )
       col = min( i, j )
       row = max( i, j )

       ind1 = ndim * ( col - 1 )
       ind2 = ndim * ( row - 1 )

       xdiff(1:ndim) = x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim)
       fparc = ( 2.0d0 * r ) ** 2 - sum( xdiff(1:ndim) ** 2 )

       if ( fparc .gt. 0.0d0 ) then
          do k = 1,ndim
             val = 8.0d0 * fparc ** 2 * ( 6.0d0 * xdiff(k) ** 2 - fparc )
             diagblock(i,k,k) = diagblock(i,k,k) + val
             diagblock(j,k,k) = diagblock(j,k,k) + val

             if ( hnnz + 1 .gt. lim ) then
                lmem = .true.
                return
             end if

             hrow(hnnz+1) = ind2 + k
             hcol(hnnz+1) = ind1 + k
             hval(hnnz+1) = - val
             hnnz = hnnz + 1

             do l = 1,k - 1
                val = 48.0d0 * fparc ** 2 * xdiff(k) * xdiff(l)
                diagblock(i,k,l) = diagblock(i,k,l) + val
                diagblock(j,k,l) = diagblock(j,k,l) + val

                if ( hnnz + 1 .gt. lim ) then
                   lmem = .true.
                   return
                end if

                hrow(hnnz+1) = ind2 + k
                hcol(hnnz+1) = ind1 + l
                hval(hnnz+1) = - val
                hnnz = hnnz + 1
             end do

             do l = k + 1,ndim
                val = 48.0d0 * fparc ** 2 * xdiff(k) * xdiff(l)

                if ( hnnz + 1 .gt. lim ) then
                   lmem = .true.
                   return
                end if

                hrow(hnnz+1) = ind2 + k
                hcol(hnnz+1) = ind1 + l
                hval(hnnz+1) = - val
                hnnz = hnnz + 1
             end do
          end do
       end if

       j = next(j)
    end do
    
  end subroutine hnrdist

  ! ******************************************************************
  ! ******************************************************************

  subroutine tnrdist(n,x,i,lstart,diagblock,tid1,tid2,tid3,tval,tnnz,&
       lim,lmem)

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: i,lim,lstart,n
    logical, intent(inout) :: lmem
    integer, intent(inout) :: tnnz

    ! ARRAY ARGUMENTS
    integer,      dimension(lim),                 intent(out)   :: tid1, tid2, tid3
    real(kind=8), dimension(n),                   intent(in)    :: x
    real(kind=8), dimension(nite,ndim,ndim,ndim), intent(inout) :: diagblock
    real(kind=8), dimension(lim),                 intent(out)   :: tval
    
    ! LOCAL SCALARS
    integer      :: col, ind1, ind2, j, k, l, p, row
    real(kind=8) :: fparc, val

    ! LOCAL ARRAYS
    real(kind=8), dimension(ndim) :: xdiff

    lmem = .false.
    
    j = lstart
    do while ( j .ne. 0 )
       col = min( i, j )
       row = max( i, j )

       ind1 = ndim * ( col - 1 )
       ind2 = ndim * ( row - 1 )

       xdiff(1:ndim) = x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim)
       fparc = ( 2.0d0 * r ) ** 2 - sum( xdiff(1:ndim) ** 2 )

       if ( fparc .gt. 0.0d0 ) then
          do k = 1, ndim
             val = 48.0d0 * fparc * xdiff(k) * ( 3.0d0 * fparc - 4.0d0 * xdiff(k) ** 2 )
             diagblock(i,k,k,k) = diagblock(i,k,k,k) + val
             diagblock(j,k,k,k) = diagblock(j,k,k,k) - val

             if ( tnnz + 2 .gt. lim ) then
                lmem = .true.
                return
             end if

             tnnz = tnnz + 1
             tid1(tnnz) = ind2 + k
             tid2(tnnz) = ind1 + k
             tid3(tnnz) = ind1 + k
             tval(tnnz) = - val

             tnnz = tnnz + 1
             tid1(tnnz) = ind2 + k
             tid2(tnnz) = ind2 + k
             tid3(tnnz) = ind1 + k
             tval(tnnz) = val

             do l = 1, ndim
                if ( l .ne. k ) then
                   if ( tnnz + 8 .gt. lim ) then
                      lmem = .true.
                      return
                   end if

                   val = 48.0d0 * fparc * xdiff(l) * ( fparc - 4.0d0 * xdiff(k) ** 2 )

                   if ( l .lt. k ) then
                      diagblock(i,k,k,l) = diagblock(i,k,k,l) + val
                      diagblock(j,k,k,l) = diagblock(j,k,k,l) - val
                   else
                      diagblock(i,l,k,k) = diagblock(i,l,k,k) + val
                      diagblock(j,l,k,k) = diagblock(j,l,k,k) - val
                   end if

                   tnnz = tnnz + 1
                   tid1(tnnz) = ind2 + l
                   tid2(tnnz) = ind1 + k
                   tid3(tnnz) = ind1 + k
                   tval(tnnz) = - val
                   
                   tnnz = tnnz + 1
                   tid1(tnnz) = ind2 + k
                   tid2(tnnz) = max( ind1 + k, ind1 + l )
                   tid3(tnnz) = min( ind1 + k, ind1 + l )
                   tval(tnnz) = - val

                   tnnz = tnnz + 1
                   tid1(tnnz) = max( ind2 + k, ind2 + l )
                   tid2(tnnz) = min( ind2 + k, ind2 + l )
                   tid3(tnnz) = ind1 + k
                   tval(tnnz) = val

                   tnnz = tnnz + 1
                   tid1(tnnz) = ind2 + k
                   tid2(tnnz) = ind2 + k
                   tid3(tnnz) = ind1 + l
                   tval(tnnz) = val

                   if ( l .lt. k ) then
                      do p = 1, l-1
                         if ( tnnz + 6 .gt. lim ) then
                            lmem = .true.
                            return
                         end if
                         
                         val = 192.0d0 * fparc * xdiff(k) * xdiff(l) * xdiff(p)
                         
                         diagblock(i,k,l,p) = diagblock(i,k,l,p) - val
                         diagblock(j,k,l,p) = diagblock(j,k,l,p) + val
                         
                         tnnz = tnnz + 1
                         tid1(tnnz) = ind2 + p
                         tid2(tnnz) = ind1 + k
                         tid3(tnnz) = ind1 + l
                         tval(tnnz) = val

                         tnnz = tnnz + 1
                         tid1(tnnz) = ind2 + l
                         tid2(tnnz) = ind1 + k
                         tid3(tnnz) = ind1 + p
                         tval(tnnz) = val

                         tnnz = tnnz + 1
                         tid1(tnnz) = ind2 + l
                         tid2(tnnz) = ind2 + p
                         tid3(tnnz) = ind1 + k
                         tval(tnnz) = - val

                         tnnz = tnnz + 1
                         tid1(tnnz) = ind2 + k
                         tid2(tnnz) = ind1 + l
                         tid3(tnnz) = ind1 + p
                         tval(tnnz) = val

                         tnnz = tnnz + 1
                         tid1(tnnz) = ind2 + k
                         tid2(tnnz) = ind2 + p
                         tid3(tnnz) = ind1 + l
                         tval(tnnz) = - val

                         tnnz = tnnz + 1
                         tid1(tnnz) = ind2 + k
                         tid2(tnnz) = ind2 + l
                         tid3(tnnz) = ind1 + p
                         tval(tnnz) = - val
                      end do
                   end if
                end if
             end do
          end do
       end if

       j = next(j)
    end do

  end subroutine tnrdist

  ! *****************************************************************
  ! *****************************************************************
  
  subroutine drawsol(n,x,unit)
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,unit

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: i, ind

    if ( ndim .ne. 3 ) then
       return
    end if

    write (unit,'(A)') 'mol new'
    write (unit,'(A)') 'draw material Opaque'

    do i = 1, nite
       ind = ( i - 1 ) * ndim
       write (unit,10) 'draw color ', modulo(i,17)
       write (unit,20) 'draw sphere {', x(ind+1:ind+ndim), '} radius', r, 'resolution 20'
    end do

    write (unit,'(A)') 'mol new'
    write (unit,'(A)') 'draw material Transparent'
    write (unit,'(A)') 'draw color 800'
    write (unit,20)    'draw sphere {', [0.0, 0.0, 0.0], '} radius', objrad, 'resolution 20'
    write (unit,'(A)') 'display resetview'
    write (unit,'(A)') 'light 0 off'
    write (unit,'(A)') 'light 1 on'
    write (unit,'(A)') 'light 2 on'
    write (unit,'(A)') 'light 3 off'
    write (unit,'(A)') 'axes location off'
    write (unit,'(A)') 'display projection orthographic'

! NONEXECUTABLE STATEMENTS
10 format( A, I6 )
20 format( A, 3(1X, F20.10), 1X, A, 1X, F20.10, 1X, A )
  end subroutine drawsol

  ! *****************************************************************
  ! *****************************************************************
  
  subroutine drawsol2(n,x,unit)
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n, unit

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: i, ind

    ! LOCAL ARRAY
    character(len=10), dimension(8) :: colors

    real(kind=8) :: drand, seed

    seed = 123456.0d0
    
    colors(1) = 'Red'
    colors(2) = 'Green'
    colors(3) = 'Blue'
    colors(4) = 'Yellow'
    colors(5) = 'Orange'
    colors(6) = 'Brown'
    colors(7) = 'Navy'

    write (unit,10) '#include "/Users/john/Software/POVRAY-3.7/include/colors.inc"'
    write (unit,*)
    write (unit,10) 'background { color White }'
    write (unit,*)
    write (unit,10) 'camera {'
    write (unit,20) '  location <0, 0, -', ceiling(2.5*objrad), '>'
    write (unit,10) '  look_at <0, 0, 0>'
    write (unit,10) '}'
    write (unit,*)

    write (unit,20) 'light_source { <0, 0, -', ceiling(2.5*objrad), '> color White}'
    write (unit,*)

    write (unit,10) 'sphere {'
    write (unit,50) '  <0, 0, 0>, ', objrad
    write (unit,10) '  texture {'
    write (unit,40) '    pigment { color Gray05 transmit 0.8 }'
    write (unit,10) '  }'
    write (unit,10) '}'
    write (unit,*)
       
    do i = 1, nite
       ind = ndim * ( i - 1 )
       write (unit,10) 'sphere {'
       write (unit,30) '  <', x(ind+1), ', ', x(ind+2), ', ', x(ind+3), '>, 1'
       write (unit,10) '  texture {'
       write (unit,40) '    pigment { color ', trim(colors(ceiling(6*drand(seed)))), ' }'
       write (unit,10) '  }'
       write (unit,10) '}'
       write (unit,*)
    end do

    ! NONEXECUTABLE STATEMENTS
10  format( A )
20  format( A, I0, A )
30  format( A, 3(F0.10, A) )
40  format( 3(A) )
50  format( A, F0.10 )
  end subroutine drawsol2

  
end program pack
