
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,ngp,igpig,wfmt,ld,wfgp,vmat)
use modmain
use moddftu
use modomp
implicit none
! arguments
! the potential is multiplied by the radial integration weights in the
! muffin-tin and by the characteristic function in the interstitial region
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci,nrco
integer npc,npc2,ipco,igp,nthd
! automatic arrays
complex(4) wfmt2(npcmtmax)
complex(8) wfmt3(npcmtmax),z(ngkmax)
! allocatable arrays
complex(4), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfir(:)
! external functions
real(4), external :: sdot
real(8), external :: ddot
complex(4), external :: cdotc
complex(8), external :: zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(npcmtmax,nstsv))
call holdthd(nstsv,nthd)
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    nrco=nrc-nrci
    npc=npcmt(is)
    npc2=npc*2
    ipco=npcmti(is)+1
! make single-precision copy of wavefunction
    wfmt1(1:npc,:)=wfmt(1:npc,ias,ispn,:)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt2,wfmt3,ist) &
!$OMP NUM_THREADS(nthd)
    do jst=1,nstsv
      wfmt2(1:npc)=vmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
! apply muffin-tin DFT+U potential matrix if required
      if (dftu.ne.0) then
        if (any(tvmmt(0:lmaxdm,ias))) then
          call zgemm('N','N',lmmaxi,nrci,lmmaxi,zone,vmatmti(:,:,1,1,ias), &
           lmmaxi,wfmt(:,ias,ispn,jst),lmmaxi,zzero,wfmt3,lmmaxi)
          call zgemm('N','N',lmmaxo,nrco,lmmaxo,zone,vmatmto(:,:,1,1,ias), &
           lmmaxo,wfmt(ipco,ias,ispn,jst),lmmaxo,zzero,wfmt3(ipco),lmmaxo)
          call zfcmtwr(nrc,nrci,wrcmt(:,is),wfmt3)
          wfmt2(1:npc)=wfmt2(1:npc)+wfmt3(1:npc)
        end if
      end if
! compute the inner products
      do ist=1,jst-1
        vmat(ist,jst)=vmat(ist,jst)+cdotc(npc,wfmt1(:,ist),1,wfmt2,1)
      end do
      vmat(jst,jst)=vmat(jst,jst)+sdot(npc2,wfmt1(:,jst),1,wfmt2,1)
    end do
!$OMP END PARALLEL DO
  end do
end do
call freethd(nthd)
deallocate(wfmt1)
!---------------------------!
!     interstitial part     !
!---------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,z,ispn,jspn) &
!$OMP PRIVATE(igp,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir(ngtot))
!$OMP DO
do jst=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform wavefunction to real-space
    wfir(:)=0.d0
    do igp=1,ngp(jspn)
      wfir(igfft(igpig(igp,jspn)))=wfgp(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir)
! apply potential to wavefunction
    wfir(1:ngtot)=vir(1:ngtot)*wfir(1:ngtot)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir)
    do igp=1,ngp(jspn)
      z(igp)=wfir(igfft(igpig(igp,jspn)))
    end do
    do ist=1,jst-1
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+zdotc(ngp(jspn),wfgp(:,ispn,ist),1,z,1)
    end do
    vmat(jst,jst)=vmat(jst,jst)+ddot(ngp(jspn)*2,wfgp(:,ispn,jst),1,z,1)
  end do
end do
!$OMP END DO
deallocate(wfir)
!$OMP END PARALLEL
call freethd(nthd)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
end subroutine

