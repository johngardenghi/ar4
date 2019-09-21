module arp_types

  implicit none
  
  ! PRECISION
  integer, parameter :: dp      = kind( 0.0d0 )
  integer, parameter :: macheps = epsilon( 0.0d0 )

  ! ------------------------------------------------------------------

  ! Sparse Matrix
  type sparse_matrix
     ! SCALAR FIELDS
     integer :: nnz

     ! ARRAY FIELDS
     integer,       allocatable, dimension(:) :: col
     integer,       allocatable, dimension(:) :: row
     real(kind=dp), allocatable, dimension(:) :: val
  end type sparse_matrix

  ! Sparse Tensor
  type sparse_tensor
     ! SCALAR FIELDS
     integer :: nnz

     ! ARRAY FIELDS
     integer,       allocatable, dimension(:) :: id1
     integer,       allocatable, dimension(:) :: id2
     integer,       allocatable, dimension(:) :: id3
     real(kind=dp), allocatable, dimension(:) :: val
  end type sparse_tensor

  ! ------------------------------------------------------------------

  ! Interfaces to functional subroutines
  interface
     subroutine evalf(x,f,flag)
       import :: dp
       ! SCALAR ARGUMENT
       integer,       intent(out) :: flag
       real(kind=dp), intent(out) :: f
       ! ARRAY ARGUMENT
       real(kind=dp), dimension(:), intent(in) :: x
     end subroutine evalf
     
     subroutine evalg(x,g,flag)
       import :: dp
       ! SCALAR ARGUMENT
       integer, intent(out) :: flag
       ! ARRAY ARGUMENTS
       real(kind=dp), dimension(:), intent(in)  :: x
       real(kind=dp), dimension(:), intent(out) :: g
     end subroutine evalg
     
     subroutine evalh(x,h,flag)
       import :: dp
       ! SCALAR ARGUMENT
       integer, intent(out) :: flag
       ! ARRAY ARGUMENTS
       real(kind=dp), dimension(:),   intent(in)  :: x
       real(kind=dp), dimension(:,:), intent(out) :: h
     end subroutine evalh

     subroutine evalhs(x,hs,flag)
       import :: dp, sparse_matrix
       ! SCALAR ARGUMENT
       integer, intent(out) :: flag
       ! ARRAY ARGUMENT
       real(kind=dp), dimension(:), intent(in) :: x
       ! SPARSE MATRIX ARGUMENT
       type(sparse_matrix), intent(inout) :: hs
     end subroutine evalhs

     subroutine evalt(x,t,flag)
       import :: dp 
       ! SCALAR ARGUMENT
       integer, intent(out) :: flag
       ! ARRAY ARGUMENTS
       real(kind=dp), dimension(:),     intent(in)  :: x
       real(kind=dp), dimension(:,:,:), intent(out) :: t
     end subroutine evalt

     subroutine evalts(x,ts,flag)
       import :: dp, sparse_tensor
       ! SCALAR ARGUMENT
       integer, intent(out) :: flag
       ! ARRAY ARGUMENT
       real(kind=dp), dimension(:), intent(in) :: x
       ! SPARSE TENSOR ARGUMENT
       type(sparse_tensor), intent(inout) :: ts       
     end subroutine evalts
  end interface

  ! ------------------------------------------------------------------

  type arp_param

     ! ALGORITHMIC PARAMETERS
     integer       :: p           ! Order of available derivatives
     integer       :: bigJ        ! Step control iterations limit
     real(kind=dp) :: alpha       ! Sufficient decrease scale factor
     real(kind=dp) :: eta1        ! Step control - function related
     real(kind=dp) :: eta2        ! Step control - step related
     real(kind=dp) :: gamma1      ! Scale factor for decreasing sigma
     real(kind=dp) :: gamma2      ! Scale factor for increasing sigma
     real(kind=dp) :: sigmalow    ! Lower bound for sigma(ini)
     real(kind=dp) :: theta       ! Subproblem stopping criterion

     ! IMPLEMENTATION PARAMETERS
     integer       :: hnnzmax     ! Maximum number of non-null elements in the Hessian
     integer       :: kmax        ! Maximum number of iterations
     integer       :: printlevel  ! Level of printing
     integer       :: tnnzmax     ! Maximum number of non-null elements in the Tensor
     integer       :: trackunit   ! Output unit to print tracking info
     logical       :: dense       ! Dense routines if true, sparse otherwise
     logical       :: sigmaini    ! Increase sigma from sigmaini?
     real(kind=dp) :: epsilon     ! Optimality precision
     real(kind=dp) :: epsilonrel  ! Optimality precision

     ! Routines
     procedure(evalf),  pointer, nopass :: uevalf
     procedure(evalg),  pointer, nopass :: uevalg
     procedure(evalh),  pointer, nopass :: uevalh
     procedure(evalhs), pointer, nopass :: uevalhs
     procedure(evalt),  pointer, nopass :: uevalt
     procedure(evalts), pointer, nopass :: uevalts

  end type arp_param

  ! ------------------------------------------------------------------
  
  type arp_output

     integer       :: fcnt       ! Total function evaluations
     integer       :: gcnt       ! Total gradient evaluations
     integer       :: hcnt       ! Total Hessian evaluations
     integer       :: info       ! Main program returning flag
     integer       :: itcnt      ! Total outer iterations
     integer       :: nwstep     ! Total iterations in which sigma=0
     integer       :: tcnt       ! Total tensor evaluations
     real(kind=4)  :: time       ! Total execution time
     real(kind=dp) :: fval       ! Objective function at the final point
     real(kind=dp) :: gnorm      ! Gradient norm at the final point
     real(kind=dp) :: sigma      ! Final regularization parameter value

  end type arp_output

end module arp_types
