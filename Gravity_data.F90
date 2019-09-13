!!****if* source/physics/Gravity/GravityMain/BinaryNonInertial/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for binary gravity. 
!!
!! PARAMETERS
!!   grv_totmass, grv_q, grv_sma
!!
!!***

module Gravity_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: grv_totmass, grv_q, grv_sma

  !! *** Binary components  *** !!
  !! (1) ~ primary, (2) ~ secondary

  real, dimension(2), save :: grv_factor, grv_ptmass, grv_ptxpos, grv_ptypos, grv_ptzpos
  real, save :: grv_period, grv_w, grv_w2

  !! *** Physical Constants *** !!

  real, save :: grv_newton


  integer, save :: grv_meshMe, grv_meshNumProcs
  logical, save :: useGravity

end module Gravity_data
