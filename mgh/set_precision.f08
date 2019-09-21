! ********************************************************************
! This module contains definitions about real kind precision used in
! every module of this package
! ********************************************************************

module set_precision
  
  implicit none

  integer, parameter :: kind_s  = kind(0.0e0)             ! Single precision
  integer, parameter :: kind_d  = kind(0.0d0)             ! Double precision
  integer, parameter :: kind_q  = selected_real_kind(30)  ! Quad precision

  ! CHANGE THIS VARIABLE TO DEFINE THE PRECISION OF THE MODULE
  integer, parameter :: rk = kind_d

end module set_precision
