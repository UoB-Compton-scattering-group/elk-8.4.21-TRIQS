
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfuv(evecu,evecv,wfmt,wfir,wfumt,wfuir,wfvmt,wfvir)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: evecu(nstsv,nstsv),evecv(nstsv,nstsv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngtc,nspinor,nstsv)
complex(8), intent(out) :: wfumt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfuir(ngtc,nspinor,nstsv)
complex(8), intent(out) :: wfvmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfvir(ngtc,nspinor,nstsv)
! local variables
integer ist,jst,ispn
integer is,ias,npc,nthd
complex(8) z1
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ist,z1,ispn) &
!$OMP PRIVATE(ias,is,npc) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do jst=1,nstsv
  wfumt(:,:,:,jst)=0.d0
  wfuir(:,:,jst)=0.d0
  do ist=1,nstsv
    z1=evecu(ist,jst)
    if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
      do ispn=1,nspinor
        do ias=1,natmtot
          is=idxis(ias)
          npc=npcmt(is)
          call zaxpy(npc,z1,wfmt(:,ias,ispn,ist),1,wfumt(:,ias,ispn,jst),1)
        end do
        call zaxpy(ngtc,z1,wfir(:,ispn,ist),1,wfuir(:,ispn,jst),1)
      end do
    end if
  end do
end do
!$OMP END DO NOWAIT
!$OMP DO
do jst=1,nstsv
  wfvmt(:,:,:,jst)=0.d0
  wfvir(:,:,jst)=0.d0
  do ist=1,nstsv
    z1=evecv(ist,jst)
    if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
      do ispn=1,nspinor
        do ias=1,natmtot
          is=idxis(ias)
          npc=npcmt(is)
          call zaxpy(npc,z1,wfmt(:,ias,ispn,ist),1,wfvmt(:,ias,ispn,jst),1)
        end do
        call zaxpy(ngtc,z1,wfir(:,ispn,ist),1,wfvir(:,ispn,jst),1)
      end do
    end if
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

