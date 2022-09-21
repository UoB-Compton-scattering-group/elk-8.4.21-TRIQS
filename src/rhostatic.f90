
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhostatic
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer is,ias,np,i
! allocatable arrays
real(8), allocatable :: jrmt0(:,:,:),jrir0(:,:)
! external functions
real(8), external :: rfmtint,rfint
! store original parameters
maxscl0=maxscl
afieldc0(:)=afieldc(:)
tforce0=tforce
tjr0=tjr
! initialise global variables
call init0
! only one self-consistent loop required
maxscl=1
! read potential from STATE.OUT
trdstate=.true.
! calculate current density with zero A-field
tjr=.true.
afieldc(:)=0.d0
! switch off forces
tforce=.false.
! run the ground-state calculation
call gndstate
! store the current density
allocate(jrmt0(npmtmax,natmtot,3),jrir0(ngtot,3))
do i=1,3
  call rfcopy(jrmt(:,:,i),jrir(:,i),jrmt0(:,:,i),jrir0(:,i))
end do
! allocate static density and charge global arrays
if (allocated(rhosmt)) deallocate(rhosmt)
allocate(rhosmt(npmtmax,natmtot,3))
if (allocated(rhosir)) deallocate(rhosir)
allocate(rhosir(ngtot,3))
if (allocated(chgsmt)) deallocate(chgsmt)
allocate(chgsmt(natmtot,3))
! loop over three directions
do i=1,3
  afieldc(:)=0.d0
  afieldc(i)=1.d0
! run the ground-state calculation
  call gndstate
! muffin-tin static density
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    rhosmt(1:np,ias,i)=rhomt(1:np,ias) &
     -solsc*(jrmt(1:np,ias,i)-jrmt0(1:np,ias,i))
! compute the muffin-tin static charge
    chgsmt(ias,i)=rfmtint(nrmt(is),nrmti(is),wrmt(:,is),rhosmt(:,ias,i))
  end do
! interstitial static density
  rhosir(:,i)=rhoir(:)-solsc*(jrir(:,i)-jrir0(:,i))
! compute the static charge
  chgstot(i)=rfint(rhosmt(:,:,i),rhosir(:,i))
end do
! write static density and charge to file
if (mp_mpi) then
  open(100,file='RHOSTAT.OUT',form='UNFORMATTED',action='WRITE')
  write(100) natmtot
  write(100) npmtmax
  write(100) ngtot
  write(100) rhosmt,rhosir
  write(100) chgsmt
  write(100) chgstot
  close(100)
  open(50,file='CHGSTAT.OUT',form='FORMATTED',action='WRITE')
  write(50,'(3G18.10)') chgstot(:)
  close(50)
end if
deallocate(jrmt0,jrir0)
! restore original input parameters
maxscl=maxscl0
afieldc(:)=afieldc0(:)
tforce=tforce0
tjr=tjr0
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

