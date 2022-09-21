
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vclcore(wfmt,vmat)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist1,ist2,ist3
integer is,ia,ias,m,nthd
integer nrc,nrci,npc
! automatic arrays
complex(8) wfcr(npcmtmax,2),zfmt(npcmtmax),v(nstsv,nstsv)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:)
! external functions
complex(8), external :: zdotc
allocate(zrhomt(npcmtmax,nstsv))
v(:,:)=0.d0
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m,npcmtmax,wfcr)
          call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt) &
!$OMP NUM_THREADS(nthd)
          do ist1=1,nstsv
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call zrho2(npc,wfcr,wfcr(:,2),wfmt(:,ias,1,ist1), &
               wfmt(:,ias,2,ist1),zfmt)
            else
              call zrho1(npc,wfcr,wfmt(:,ias,1,ist1),zfmt)
            end if
            call zfsht(nrc,nrci,zfmt,zrhomt(:,ist1))
          end do
!$OMP END PARALLEL DO
          call freethd(nthd)
          call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,ist1) &
!$OMP NUM_THREADS(nthd)
          do ist2=1,nstsv
            call zpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
             zrhomt(:,ist2),zfmt)
            call zfmtwr(nrc,nrci,wrcmt(:,is),zfmt)
            do ist1=1,ist2
              v(ist1,ist2)=v(ist1,ist2)-zdotc(npc,zrhomt(:,ist1),1,zfmt,1)
            end do
          end do
!$OMP END PARALLEL DO
          call freethd(nthd)
        end do
      end if
    end do
  end do
end do
! set the lower triangular part of the matrix
do ist1=1,nstsv
  do ist2=1,ist1-1
    v(ist1,ist2)=conjg(v(ist2,ist1))
  end do
end do
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) v(:,:)=hybridc*v(:,:)
! add to input matrix
vmat(:,:)=vmat(:,:)+v(:,:)
deallocate(zrhomt)
return

contains

pure subroutine zrho1(n,x,y,z)
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: x(n),y(n)
complex(8), intent(out) :: z(n)
z(:)=conjg(x(:))*y(:)
end subroutine

pure subroutine zrho2(n,x1,x2,y1,y2,z)
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: x1(n),x2(n),y1(n),y2(n)
complex(8), intent(out) :: z(n)
z(:)=conjg(x1(:))*y1(:)+conjg(x2(:))*y2(:)
end subroutine

end subroutine

