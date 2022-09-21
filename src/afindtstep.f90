
! Copyright (C) 2020 Peter Elliott, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine afindtstep
use modmain
use modtddft
use modmpi
implicit none
! local variables
real(8) dt,t1
! time step length
dt=times(itimes+1)-times(itimes)
! add to the time derivative of the induced A-field
t1=fourpi*solsc*afindscf*dt/omega
afindt(:,1)=afindt(:,1)+t1*jtot(:)
! add to the induced A-field
afindt(:,0)=afindt(:,0)+afindt(:,1)*dt
! add to the total A-field
afieldt(:,itimes+1)=afieldt(:,itimes+1)+afindt(:,0)
! write the induced A-field and its time derivative to file
if (mp_mpi) then
  open(50,file='AFINDT.OUT',form='FORMATTED')
  write(50,'(6G18.10)') afindt(:,:)
  close(50)
end if
end subroutine

