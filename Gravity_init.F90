!!****if* source/physics/Gravity/GravityMain/BinaryNonInertial/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes the gravitational physics unit for BinaryNonInertial.
!!  Assumes that the orbital axis is "z" and the the binary is aligned along the "x" axis.
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Gravity_init()

  use Gravity_data
  use Grid_Data, ONLY: gr_globalMe
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"

  real :: pi

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshMe)

  call PhysicalConstants_get("newton", grv_newton)
  call PhysicalConstants_get("Pi", pi)

  call RuntimeParameters_get("totmass", grv_totmass)
  call RuntimeParameters_get("q", grv_q)
  call RuntimeParameters_get("sma", grv_sma)
  call RuntimeParameters_get("useGravity", useGravity)

!==============================================================================

  ! We do the following normalization because this module is intended for Bondi
  ! simulations. Change this (you could do if Bondi == true or something) if you 
  ! do not want to use Bondi units. 
  ! Set sim_totmass = 1/G such that at cinf = 1, bondi radius = 1
  grv_totmass = 1.0/grv_newton

  ! Set binary parameters from runtime parameters; (1) is primary, (2) is secondary

  grv_ptmass(1) = grv_totmass/(1.0 + grv_q)
  grv_ptmass(2) = grv_totmass*grv_q/(1.0 + grv_q)

  ! Semimajor axis along x-axis. NK: Add runtime parameter to change this!
  grv_ptxpos(1) = -grv_sma*grv_ptmass(2)/grv_totmass
  grv_ptxpos(2) = grv_sma*grv_ptmass(1)/grv_totmass
  grv_ptypos(1) = 0.0
  grv_ptypos(2) = 0.0

  ! Orbital angular momentum in z direction; keep in z plane
  grv_ptzpos(1) = 0.0
  grv_ptzpos(2) = 0.0

  ! Gravitational potential for each component
  grv_factor(1) = -grv_newton * grv_ptmass(1)
  grv_factor(2) = -grv_newton * grv_ptmass(2)

  ! Other orbital parameters
  grv_period = 2.0*pi*sqrt((grv_sma**3.0)/(grv_newton*grv_totmass))
  grv_w      = 2.0*pi/grv_period
  grv_w2     = grv_w*grv_w

!==============================================================================

  if (gr_globalMe  == MASTER_PE) then


     print*,"######## Checking Binary parameters  ##################"
     print*,"grv_totmass (1/G)", grv_newton*grv_totmass
     print*,"grv_q", grv_q
     print*,"grv_sma", grv_sma
     print*,"Primary mass (1/G)", grv_newton*grv_ptmass(1)
     print*,"Secondary mass (1/G)", grv_newton*grv_ptmass(2)
     print*,"Primary xpos", grv_ptxpos(1)
     print*,"Primary ypos", grv_ptypos(1)
     print*,"Primary zpos", grv_ptzpos(1)
     print*,"Secondary xpos", grv_ptxpos(2)
     print*,"Secondary ypos", grv_ptypos(2)
     print*,"Secondary zpos", grv_ptzpos(2)
     print*,"Binary period",  grv_period
     print*,"##############################################"

  endif


  return
end subroutine Gravity_init
