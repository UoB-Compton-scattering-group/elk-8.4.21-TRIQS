
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dwfmtfv(lrstp,ias,ngp,ngpq,apwalmq,dapwalm,evecfv,devecfv,dwfmt)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: lrstp,ias,ngp,ngpq
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax),devecfv(nmatmax)
complex(8), intent(out) :: dwfmt(*)
! local variables
integer is,io,ilo
integer nri,nro,iro
integer l,lm,np,npi,i
complex(8) z1
! external functions
complex(8), external :: zdotu
is=idxis(ias)
iro=nrmti(is)+lrstp
if (lrstp.eq.1) then
  nri=nrmti(is)
  nro=nrmt(is)-nri
  np=npmt(is)
  npi=npmti(is)
else
  nri=nrcmti(is)
  nro=nrcmt(is)-nri
  np=npcmt(is)
  npi=npcmti(is)
end if
! zero the wavefunction derivative
dwfmt(1:np)=0.d0
!-----------------------!
!     APW functions     !
!-----------------------!
do l=0,lmaxo
  do lm=l**2+1,(l+1)**2
    i=npi+lm
    do io=1,apword(l,is)
      z1=zdotu(ngpq,devecfv,1,apwalmq(:,io,lm),1)
      if (ias.eq.iasph) then
        z1=z1+zdotu(ngp,evecfv,1,dapwalm(:,io,lm),1)
      end if
      if (l.le.lmaxi) then
        call zfzrf(nri,z1,lrstp,apwfr(1,1,io,l,ias),lmmaxi,dwfmt(lm))
      end if
      call zfzrf(nro,z1,lrstp,apwfr(iro,1,io,l,ias),lmmaxo,dwfmt(i))
    end do
  end do
end do
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    i=npi+lm
    z1=devecfv(ngpq+idxlo(lm,ilo,ias))
    if (l.le.lmaxi) then
      call zfzrf(nri,z1,lrstp,lofr(1,1,ilo,ias),lmmaxi,dwfmt(lm))
    end if
    call zfzrf(nro,z1,lrstp,lofr(iro,1,ilo,ias),lmmaxo,dwfmt(i))
  end do
end do
return

contains

pure subroutine zfzrf(n,z,ld1,rf,ld2,zf)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: z
integer, intent(in) :: ld1
real(8), intent(in) :: rf(ld1,n)
integer, intent(in) :: ld2
complex(8), intent(inout) :: zf(ld2,n)
zf(1,:)=zf(1,:)+z*rf(1,:)
end subroutine

end subroutine

