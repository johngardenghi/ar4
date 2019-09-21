subroutine allocatearrays(n,m,jcnnzmax,hnnzmax,allocerr)

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: hnnzmax,jcnnzmax,m,n
  integer, intent(out) :: allocerr

  call allocatesvdgrad(n,m,jcnnzmax,allocerr)
  if ( allocerr .ne. 0 ) return

  call allocatesvdhess(hnnzmax,allocerr)
  if ( allocerr .ne. 0 ) return

  call allocateminsq(n,m,jcnnzmax,allocerr)
  if ( allocerr .ne. 0 ) return

end subroutine allocatearrays

! *****************************************************************
! *****************************************************************

subroutine deallocatearrays(deallocerr)

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  call deallocatesvdgrad(deallocerr)
  if ( deallocerr .ne. 0 ) return

  call deallocatesvdhess(deallocerr)
  if ( deallocerr .ne. 0 ) return

  call deallocateminsq(deallocerr)
  if ( deallocerr .ne. 0 ) return

end subroutine deallocatearrays

! *****************************************************************
! *****************************************************************

subroutine allocategencanarrays(n,allocerr)

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n
  integer, intent(out) :: allocerr

  call allocaterspace(n,allocerr)
  if ( allocerr .ne. 0 ) return

  call allocatesydat(n,allocerr)
  if ( allocerr .ne. 0 ) return

  call allocatehappdat(n,allocerr)
  if ( allocerr .ne. 0 ) return

  call allocatehpredat(n,allocerr)
  if ( allocerr .ne. 0 ) return

end subroutine allocategencanarrays

! *****************************************************************
! *****************************************************************

subroutine deallocategencanarrays(deallocerr)

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  call deallocaterspace(deallocerr)
  if ( deallocerr .ne. 0 ) return

  call deallocatesydat(deallocerr)
  if ( deallocerr .ne. 0 ) return

  call deallocatehappdat(deallocerr)
  if ( deallocerr .ne. 0 ) return

  call deallocatehpredat(deallocerr)
  if ( deallocerr .ne. 0 ) return

end subroutine deallocategencanarrays

! *****************************************************************
! *****************************************************************

subroutine allocaterspace(n,allocerr)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n
  integer, intent(out) :: allocerr

  ! Allocate global arrays in module rspace

  allocate(ind(n),xfull(n),stat=allocerr)

end subroutine allocaterspace

! *****************************************************************
! *****************************************************************

subroutine deallocaterspace(deallocerr)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module rspace

  deallocate(ind,xfull,stat=deallocerr)

end subroutine deallocaterspace

! *****************************************************************
! *****************************************************************

subroutine allocatesydat(n,allocerr)

  use modsydat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n
  integer, intent(out) :: allocerr

  ! Allocate global arrays in module sysdat

  allocate(s(n),y(n),stat=allocerr)

end subroutine allocatesydat

! *****************************************************************
! *****************************************************************

subroutine deallocatesydat(deallocerr)

  use modsydat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module sysdat

  deallocate(s,y,stat=deallocerr)

end subroutine deallocatesydat

! *****************************************************************
! *****************************************************************

subroutine allocatesvdgrad(n,m,jcnnzmax,allocerr)

  use modsvdgrad

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: jcnnzmax,m,n
  integer, intent(out) :: allocerr

  ! Allocate global arrays in module svdgrad

  allocate(svdjcvar(jcnnzmax),svdjcsta(m),svdjclen(m),svdc(m), &
       svddpdc(m),svdg(n),svdgparc(n),svdjcval(jcnnzmax),stat=allocerr)

end subroutine allocatesvdgrad

! *****************************************************************
! *****************************************************************

subroutine deallocatesvdgrad(deallocerr)

  use modsvdgrad

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module svdgrad

  deallocate(svdjcvar,svdjcsta,svdjclen,svdc,svddpdc,svdg,svdgparc, &
       svdjcval,stat=deallocerr)

end subroutine deallocatesvdgrad

! *****************************************************************
! *****************************************************************

subroutine allocatesvdhess(hnnzmax,allocerr)

  use modsvdhess

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: hnnzmax
  integer, intent(out) :: allocerr

  ! Allocate global arrays in module svdhess

  allocate(svdhcol(hnnzmax),svdhrow(hnnzmax),svdhval(hnnzmax),stat=allocerr)

end subroutine allocatesvdhess

! *****************************************************************
! *****************************************************************

subroutine deallocatesvdhess(deallocerr)

  use modsvdhess

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module svdhess

  deallocate(svdhcol,svdhrow,svdhval,stat=deallocerr)

end subroutine deallocatesvdhess

! *****************************************************************
! *****************************************************************

subroutine allocatehappdat(n,allocerr)

  use modhappdat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n
  integer, intent(out) :: allocerr

  ! Allocate global arrays in module happdat

  allocate(hds(n),stat=allocerr)

end subroutine allocatehappdat

! *****************************************************************
! *****************************************************************

subroutine deallocatehappdat(deallocerr)

  use modhappdat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module happdat

  deallocate(hds,stat=deallocerr)

end subroutine deallocatehappdat

! *****************************************************************
! *****************************************************************

subroutine allocatehpredat(n,allocerr)

  use modhpredat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: n
  integer, intent(out) :: allocerr

  ! Allocate global arrays in module hpredat

  allocate(pdiag(n),psmdy(n),stat=allocerr)

end subroutine allocatehpredat

! *****************************************************************
! *****************************************************************

subroutine deallocatehpredat(deallocerr)

  use modhpredat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module hpredat

  deallocate(pdiag,psmdy,stat=deallocerr)
  
end subroutine deallocatehpredat

! *****************************************************************
! *****************************************************************

subroutine allocateminsq(n,m,jcnnzmax,allocerr)

  use modminsq

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: jcnnzmax,m,n
  integer, intent(out) :: allocerr

  ! LOCAL SCALARS
  integer annzmax

  ! Allocate global arrays in module minsq

  ! Matrix A will hold the some data related to a min squares 
  ! problem that may be solved to correct the multipliers' signal 
  ! after running the acceleration process

  ! Matrix A will be formed considering the Jacobian of the constraints
  ! plus entrances related to slack variables for the inequality 
  ! constraints plus entrances related to the bound constraints on x.
  ! This is the reason for its maximum size below.

  annzmax = jcnnzmax + m + 2 * n

  allocate(b(n),arow(annzmax),acol(annzmax),aval(annzmax),stat=allocerr)

end subroutine allocateminsq

! *****************************************************************
! *****************************************************************

subroutine deallocateminsq(deallocerr)

  use modminsq

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(out) :: deallocerr

  ! Deallocate global arrays in module minsq

  deallocate(b,arow,acol,aval,stat=deallocerr)
  
end subroutine deallocateminsq

