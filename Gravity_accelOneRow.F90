!!****if* source/physics/Gravity/GravityMain/BinaryNonInertial/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in pointmass gravity
!!  extraAccelVars :  optional, ignored in this implementation
!! 
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!=======================================================================

  use Gravity_data, ONLY: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_factor, &
       useGravity, grv_w, grv_w2
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer,intent(IN),OPTIONAL :: extraAccelVars(MDIM)

!==========================================================================


#ifdef FIXEDBLOCKSIZE
  real,dimension(2,GRID_KHI_GC) :: zCenter
  real,dimension(2,GRID_JHI_GC) :: yCenter
  real,dimension(2,GRID_IHI_GC) :: xCenter
  real,dimension(GRID_KHI_GC) :: xPos
  real,dimension(GRID_JHI_GC) :: yPos
  real,dimension(GRID_IHI_GC) :: zPos
#else
  real,allocatable,dimension(:,:) ::xCenter,yCenter,zCenter
  real,allocatable,dimension(:) ::xPos,yPos,zPos
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
#endif
  real,dimension(2) :: dr32, tmpdr32
  real              :: tmpvel

  integer :: sizeX,sizeY,sizez

  real, pointer :: solnData(:,:,:,:)

  integer :: ii,j,k
  logical :: gcell = .true.

!==============================================================================

  if (.NOT.useGravity) return

  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(2,sizeX))
  allocate(yCenter(2,sizeY))
  allocate(zCenter(2,sizeZ))
  allocate(xPos(sizeX))
  allocate(yPos(sizeY))
  allocate(zPos(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zPos, sizeZ)
     zCenter(1,:) = zPos(:) - grv_ptzpos(1)
     zCenter(2,:) = zPos(:) - grv_ptzpos(2)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yPos, sizeY)
     yCenter(1,:) = yPos(:) - grv_ptypos(1)
     yCenter(2,:) = yPos(:) - grv_ptypos(2)
  endif
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xPos, sizeX)
  xCenter(1,:) = xPos(:) - grv_ptxpos(1)
  xCenter(2,:) = xPos(:) - grv_ptxpos(2)
 
  ! need solnData to get velocities for coriolis terms
  call Grid_getBlkPtr(blockID, solnData)

  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     tmpdr32(1) = yCenter(1,j)*yCenter(1,j) + zCenter(1,k)*zCenter(1,k)
     tmpdr32(2) = yCenter(2,j)*yCenter(2,j) + zCenter(2,k)*zCenter(2,k)

     do ii = 1, numCells

        dr32(1) = sqrt(xCenter(1,ii)*xCenter(1,ii) + tmpdr32(1))
        dr32(2) = sqrt(xCenter(2,ii)*xCenter(2,ii) + tmpdr32(2))
        dr32(1) = dr32(1)*dr32(1)*dr32(1)
        dr32(2) = dr32(2)*dr32(2)*dr32(2)

        tmpvel = solnData(VELY_VAR,ii,j,k) ! need y vel for coriolis x acceleration

        grav(ii) = grv_factor(1)*xCenter(1,ii)/dr32(1) + grv_factor(2)*xCenter(2,ii)/dr32(2) ! gravitational terms
        grav(ii) = grav(ii) + grv_w2*xPos(ii) ! centrifugal term in x direction
        grav(ii) = grav(ii) + 2.0*grv_w*tmpvel ! coriolis term in x direction
     enddo


  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32(1) = xCenter(1,j)*xCenter(1,j) + zCenter(1,k)*zCenter(1,k) 
     tmpdr32(2) = xCenter(2,j)*xCenter(2,j) + zCenter(2,k)*zCenter(2,k) 

     do ii = 1, numCells
        
        dr32(1) = sqrt(yCenter(1,ii)*yCenter(1,ii) + tmpdr32(1))
        dr32(2) = sqrt(yCenter(2,ii)*yCenter(2,ii) + tmpdr32(2))
        dr32(1) = dr32(1)*dr32(1)*dr32(1)
        dr32(2) = dr32(2)*dr32(2)*dr32(2)

        tmpvel = solnData(VELX_VAR,ii,j,k) ! need x vel for coriolis y acceleration

        grav(ii) = grv_factor(1)*yCenter(1,ii)/dr32(1) + grv_factor(2)*yCenter(2,ii)/dr32(2)
        grav(ii) = grav(ii) + grv_w2*yPos(ii)  ! centrifugal term in y direction 
        grav(ii) = grav(ii) - 2.0*grv_w*tmpvel ! coriolis term in y direction 
     enddo

  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32(1) = xCenter(1,j)*xCenter(1,j) + yCenter(1,k)*yCenter(1,k) 
     tmpdr32(2) = xCenter(2,j)*xCenter(2,j) + yCenter(2,k)*yCenter(2,k) 

     do ii = 1, numCells
        
        dr32(1) = sqrt(zCenter(1,ii)*zCenter(1,ii) + tmpdr32(1))
        dr32(2) = sqrt(zCenter(2,ii)*zCenter(2,ii) + tmpdr32(2))
        dr32(1) = dr32(1)*dr32(1)*dr32(1)
        dr32(2) = dr32(2)*dr32(2)*dr32(2)
        
        grav(ii) = grv_factor(1)*zCenter(1,ii)/dr32(1) + grv_factor(2)*zCenter(2,ii)/dr32(2)
     enddo

  endif

  call Grid_releaseBlkPtr(blockID, solnData)

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  deallocate(xPos)
  deallocate(yPos)
  deallocate(zPos)
#endif

  return

end subroutine Gravity_accelOneRow
