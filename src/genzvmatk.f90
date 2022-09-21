
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvmatk(zvmt,zvir,ngp,igpig,wfmt,wfir,wfgp,vmat)
use modmain
use modomp
implicit none
! arguments
! the potential is multiplied by the radial integration weights in the
! muffin-tin and by the characteristic function in the interstitial region
complex(8), intent(in) :: zvmt(npcmtmax,natmtot),zvir(ngtot)
integer, intent(in) :: ngp,igpig(ngp)
! note that the muffin-tin wavefunction is single-precision
complex(4), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
! note that wfir does not have a 1/sqrt(omega) prefactor
complex(8), intent(in) :: wfir(ngtot,nspinor,nstsv)
complex(8), intent(in) :: wfgp(ngp,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn
integer is,ias,npc,nthd
! automatic arrays
complex(4) wfmt1(npcmtmax)
complex(8) z(ngp)
! allocatable arrays
complex(8), allocatable :: wfir1(:)
! external functions
complex(4), external :: cdotc
complex(8), external :: zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ispn,ias) &
!$OMP PRIVATE(is,npc,ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nstsv
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
! apply complex potential to wavefunction
      wfmt1(1:npc)=zvmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
! compute the inner products
      do ist=1,nstsv
        vmat(ist,jst)=vmat(ist,jst)+cdotc(npc,wfmt(:,ias,ispn,ist),1,wfmt1,1)
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
!---------------------------!
!     interstitial part     !
!---------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,z,ispn,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir1(ngtot))
!$OMP DO
do jst=1,nstsv
  do ispn=1,nspinor
! apply potential to wavefunction
    wfir1(1:ngtot)=zvir(1:ngtot)*wfir(1:ngtot,ispn,jst)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir1)
    z(1:ngp)=wfir1(igfft(igpig(1:ngp)))
    do ist=1,nstsv
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+zdotc(ngp,wfgp(:,ispn,ist),1,z,1)
    end do
  end do
end do
!$OMP END DO
deallocate(wfir1)
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

