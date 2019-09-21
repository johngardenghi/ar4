module arp_model

  use arp_types
  
  implicit none

  ! GLOBAL SCALARS
  integer,       pointer, private :: p
  logical,       pointer, private :: dense
  real(kind=dp), pointer, private :: f, sigma

  ! GLOBAL ARRAYS
  real(kind=dp),       private, pointer :: g(:), h(:,:), t(:,:,:)
  type(sparse_matrix), private, pointer :: hs
  type(sparse_tensor), private, pointer :: ts

contains

  subroutine set_model_pointers( p_in, sigma_in, dense_in, f_in, g_in, &
       h_in, t_in, hs_in, ts_in )
    
    ! SCALAR ARGUMENTS
    integer,       intent(in), target :: p_in
    logical,       intent(in), target :: dense_in
    real(kind=dp), intent(in), target :: f_in, sigma_in

    ! ARRAY ARGUMENTS
    real(kind=dp),       intent(in), target           :: g_in(:)
    real(kind=dp),       intent(in), target, optional :: h_in(:,:), t_in(:,:,:)
    type(sparse_matrix), intent(in), target, optional :: hs_in
    type(sparse_tensor), intent(in), target, optional :: ts_in

    nullify( dense, f, g, h, p, sigma, t )

    p => p_in

    f => f_in
    g => g_in
 
    if ( present( h_in ) ) h => h_in
    if ( present( t_in ) ) t => t_in

    if ( present( hs_in ) ) hs => hs_in
    if ( present( ts_in ) ) ts => ts_in

    sigma => sigma_in
    dense => dense_in

  end subroutine set_model_pointers

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalf( n, s, m, flag )
    
    ! SCALAR ARGUMENT
    integer,       intent(in)  :: n
    integer,       intent(out) :: flag
    real(kind=dp), intent(out) :: m

    ! ARRAY ARGUMENT
    real(kind=dp), dimension(n), intent(in) :: s

    ! ==============================================================
    ! This subroutine computes the model value at a specified
    ! direction s
    ! ==============================================================

    ! LOCAL SCALAR
    integer :: i, j, k, l

    ! LOCAL ARRAY
    ! real(kind=dp), dimension(n,n) :: w

    m = f + dot_product( g, s )
    ! m = f
    if ( dense ) then
       ! m = m + 0.5_dp * dot_product( matmul( h, s ), s )
       do j = 1, n
          do i = 1, j-1
             ! v(i) = v(i) + h(i,j)*s(j)
             ! v(j) = v(j) + h(i,j)*s(i)
             m = m + s(i)*h(i,j)*s(j)
          end do
          m = m + 0.5_dp*h(j,j)*s(j)**2
       end do
    else
       do k = 1, hs%nnz
          i = hs%row(k)
          j = hs%col(k)

          if ( i .eq. j ) then
             m = m + 0.5_dp * s(i) * hs%val(k) * s(j)
          else
             m = m + s(i) * hs%val(k) * s(j)
          end if
       end do
    end if

    if ( p .eq. 2 ) then
       m = m + ( sigma / 3.0_dp ) * norm2( s ) ** 3
    elseif ( p .eq. 3 ) then
       if ( dense ) then
          do k = 1, n
             do j = 1, k
                do i = 1, j
                   if ( i .eq. j .and. j .eq. k ) then
                      m = m + ( 1.0_dp / 6.0_dp ) * t(i,j,k)*s(i)*s(j)*s(k)
                   elseif ( i .ne. j .and. i .ne. k .and. j .ne. k ) then
                      m = m + t(i,j,k)*s(i)*s(j)*s(k)
                   else
                      m = m + 0.5_dp*t(i,j,k)*s(i)*s(j)*s(k)
                   end if
                end do
             end do
          end do
          
          ! w = s(1) * t(:,:,1)
          ! do k = 2, n
          !    w = w + s(k) * t(:,:,k)
          ! end do

          ! m = m + ( 1.0_dp / 6.0_dp ) * dot_product( matmul( w, s ), s ) + &
          !      0.25_dp * sigma * dot_product( s, s ) ** 2

          m = m + 0.25_dp*sigma*dot_product( s, s )**2
       else
          do l = 1, ts%nnz
             i = ts%id1(l)
             j = ts%id2(l)
             k = ts%id3(l)

             if ( i .eq. j .and. j .eq. k ) then
                m = m + ( 1.0_dp / 6.0_dp ) * ts%val(l) * s(i) * s(j) * s(k)
             elseif ( i .ne. j .and. i .ne. k .and. j .ne. k ) then
                m = m + ts%val(l) * s(i) * s(j) * s(k)
             else
                m = m + 0.5_dp * ts%val(l) * s(i) * s(j) * s(k)
             end if
          end do

          m = m + 0.25_dp * sigma * sum( s ** 4 )
       end if
    end if

    flag = 0

  end subroutine model_evalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalg( n, s, mg, flag )

    ! SCALAR ARGUMENT
    integer, intent(in)  :: n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENT
    real(kind=dp), dimension(n), intent(in)  :: s
    real(kind=dp), dimension(n), intent(out) :: mg

    ! ==============================================================
    ! This subroutine computes the model gradient at a specified
    ! direction s
    ! ==============================================================

    ! LOCAL SCALAR
    integer :: i, j, k, l

    ! LOCAL ARRAY
    ! real(kind=dp), dimension(n,n) :: w

    mg = g

    if ( dense ) then
       do j = 1, n
          do i = 1, j-1
             mg(i) = mg(i) + h(i,j)*s(j)
             mg(j) = mg(j) + h(i,j)*s(i)
          end do
          mg(j) = mg(j) + h(j,j)*s(j)
       end do
       ! mg = mg + matmul( h, s )
    else
       do k = 1, hs%nnz
          i = hs%row(k)
          j = hs%col(k)
       
          mg(i) = mg(i) + hs%val(k) * s(j)

          if ( i .ne. j ) mg(j) = mg(j) + hs%val(k) * s(i)
       end do
    end if

    if ( p .eq. 2 ) then
       mg = mg + sigma * norm2( s ) * s
    elseif ( p .eq. 3 ) then
       if ( dense ) then
          ! w = s(1) * t(:,:,1)
          ! do k = 2, n
          !    w = w + s(k) * t(:,:,k)
          ! end do

          ! mg = mg + 0.5_dp * matmul( w, s ) + sigma * dot_product( s, s ) * s

          do k = 1, n
             do j = 1, k
                do i = 1, j
                   if ( i .ne. j .and. j .ne. k ) then
                      mg(i) = mg(i) + t(i,j,k)*s(j)*s(k)
                      mg(j) = mg(j) + t(i,j,k)*s(i)*s(k)
                      mg(k) = mg(k) + t(i,j,k)*s(i)*s(j)

                   elseif ( i .eq. j .and. j .eq. k ) then
                      mg(i) = mg(i) + 0.5_dp*t(i,j,k)*s(j)*s(k)

                   elseif ( i .eq. j ) then
                      mg(i) = mg(i) +        t(i,j,k)*s(j)*s(k)
                      mg(k) = mg(k) + 0.5_dp*t(i,j,k)*s(i)*s(j)
             
                   elseif ( j .eq. k ) then
                      mg(i) = mg(i) + 0.5*t(i,j,k)*s(j)*s(k)
                      mg(j) = mg(j) +     t(i,j,k)*s(i)*s(k)

                   end if
                end do
             end do
          end do

          mg = mg + sigma*dot_product( s, s )*s
       else
          do l = 1, ts%nnz
             i = ts%id1(l)
             j = ts%id2(l)
             k = ts%id3(l)

             if ( i .ne. j .and. j .ne. k ) then
                mg(i) = mg(i) + ts%val(l) * s(j) * s(k)
                mg(j) = mg(j) + ts%val(l) * s(i) * s(k)
                mg(k) = mg(k) + ts%val(l) * s(i) * s(j)

             elseif ( i .eq. j .and. j .eq. k ) then
                mg(i) = mg(i) + 0.5_dp * ts%val(l) * s(j) * s(k)

             elseif ( i .eq. j ) then
                mg(i) = mg(i) +          ts%val(l) * s(j) * s(k)
                mg(k) = mg(k) + 0.5_dp * ts%val(l) * s(i) * s(j)
             
             elseif ( j .eq. k ) then
                mg(i) = mg(i) + 0.5 * ts%val(l) * s(j) * s(k)
                mg(j) = mg(j) +       ts%val(l) * s(i) * s(k)
                   
             end if
          end do
       
          mg = mg + sigma * ( s ** 3 )
       end if
    end if

    flag = 0

  end subroutine model_evalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalh(n, s, hrow, hcol, hval, hnnz, lim, lmem, flag)

    ! SCALAR ARGUMENT
    integer, intent(in)  :: lim, n
    integer, intent(out) :: flag, hnnz
    logical, intent(out) :: lmem

    ! ARRAY ARGUMENT
    integer,       dimension(lim), intent(out) :: hcol, hrow
    real(kind=dp), dimension(n),   intent(in)  :: s
    real(kind=dp), dimension(lim), intent(out) :: hval

    ! ==============================================================
    ! This subroutine computes the model Hessian at a specified
    ! direction s
    ! ==============================================================

    ! LOCAL SCALAR
    integer       :: i, j, k, l
    real(kind=dp) :: snorm, snorm2

    ! LOCAL ARRAY
    real(kind=dp), dimension(n,n) :: w

    if ( dense ) then
       w(1:n,1:n) = h(1:n,1:n)

    else
       hnnz = 0
       if ( hs%nnz .gt. lim ) then
          print *, "Lack of memory in model_evalh", hs%nnz, lim
          lmem = .true.
          return
       end if
    
       do k = 1, hs%nnz
          hnnz = hnnz + 1
          hrow(hnnz) = hs%row(k)
          hcol(hnnz) = hs%col(k)
          hval(hnnz) = hs%val(k)
       end do
    end if


    if ( p .eq. 2 ) then
       snorm = norm2( s )

       do j = 1, n
          do i = 1, j-1
             w(i,j) = w(i,j) + sigma / snorm * s(i) * s(j)
          end do
          w(j,j) = w(j,j) + sigma * ( 1.0_dp / snorm * s(j) ** 2 + snorm )
       end do

    elseif ( p .eq. 3 ) then
       if ( dense ) then
          do k = 1, n
             do j = 1, k
                do i = 1, j
                   w(i,j) = w(i,j) + s(k)*t(i,j,k)
                end do
             end do
          end do
          
          ! w(1:n,1:n) = w(1:n,1:n) + s(1) * t(1:n,1:n,1)
          ! do k = 2, n
          !    w = w + s(k) * t(1:n,1:n,k)
          ! end do

          snorm2 = dot_product( s, s )

          do j = 1, n
             do i = 1, j-1
                w(i,j) = w(i,j) + 2.0_dp * sigma * s(i) * s(j)
             end do
             w(j,j) = w(j,j) + sigma * ( 2.0_dp * s(j) * s(j) + snorm2 )
          end do

       else
          do l = 1, ts%nnz
             if ( hnnz + 1 .gt. lim ) then
                print *, "Lack of memory in model_evalh", hnnz, lim
                lmem = .true.
                return
             end if

             i = ts%id1(l)
             j = ts%id2(l)
             k = ts%id3(l)

             hnnz = hnnz + 1
             hrow(hnnz) = i
             hcol(hnnz) = j
             hval(hnnz) = ts%val(l) * s(k)

             if ( j .ne. k ) then
                if ( hnnz + 1 .gt. lim ) then
                   print *, "Lack of memory in model_evalh", hnnz, lim
                   lmem = .true.
                   return
                end if
                
                hnnz = hnnz + 1
                hrow(hnnz) = i
                hcol(hnnz) = k
                hval(hnnz) = ts%val(l) * s(j)
             end if

             if ( i .ne. j ) then
                if ( hnnz + 1 .gt. lim ) then
                   print *, "Lack of memory in model_evalh", hnnz, lim
                   lmem = .true.
                   return
                end if

                hnnz = hnnz + 1
                hrow(hnnz) = j
                hcol(hnnz) = k
                hval(hnnz) = ts%val(l) * s(i)
             end if
          end do

          if ( hnnz + n .gt. lim ) then
             print *, "Lack of memory in model_evalh", hnnz, lim
             lmem = .true.
             return
          end if
          
          do i = 1, n
             hnnz = hnnz + 1
             hrow(hnnz) = i
             hcol(hnnz) = i
             hval(hnnz) = 3.0_dp * sigma * s(i) ** 2
          end do

       end if
    end if

    if ( dense ) then
       hnnz = 0
       do j = 1, n
          do i = 1, j
             if ( abs( w(i,j) ) .gt. macheps ) then
                hnnz = hnnz + 1
                if ( hnnz .gt. lim ) then
                   lmem = .true.
                   flag = - 1
                   return
                end if
                hrow(hnnz) = j
                hcol(hnnz) = i
                hval(hnnz) = w(i,j)
             end if
          end do
       end do
    end if

    flag = 0
    lmem = .false.

  end subroutine model_evalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalc(n,x,ind,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: ind,n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)

    flag = - 1

  end subroutine model_evalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in)  :: ind,lim,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    flag = - 1

  end subroutine model_evaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in)  :: ind,lim,n
    integer, intent(out) :: flag,hcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hccol(lim),hcrow(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hcval(lim)

    flag = - 1

  end subroutine model_evalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalfc(n,x,f,m,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: m,n
    integer,      intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    flag = - 1

  end subroutine model_evalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,m,n
    integer,      intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    flag = - 1

  end subroutine model_evalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalgjacp(n,x,g,m,p,q,work,gotj,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m,n
    integer,   intent(out)   :: flag
    character, intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)    :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out)   :: g(n)

    flag = - 1

  end subroutine model_evalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz, &
       lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(out) :: lmem
    integer,      intent(in)  :: lim,m,n
    integer,      intent(out) :: flag,hlnnz
    real(kind=8), intent(in)  :: sf

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hlcol(lim),hlrow(lim)
    real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hlval(lim)

    flag = - 1

  end subroutine model_evalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine model_evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(inout) :: goth
    integer,      intent(in)    :: m,n
    integer,      intent(out)   :: flag
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    flag = - 1

  end subroutine model_evalhlp

end module arp_model
