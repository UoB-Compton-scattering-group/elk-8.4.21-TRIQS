
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wfmtfv
! !INTERFACE:
subroutine wfmtfv(lrstp,ias,ngp,apwalm,evecfv,wfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   wfmt   : complex muffin-tin wavefunction passed in as real array
!            (out,real(2,*))
! !DESCRIPTION:
!   Calculates the first-variational wavefunction in the muffin-tin in terms of
!   a spherical harmonic expansion. For atom $\alpha$ and a particular $k$-point
!   ${\bf p}$, the $r$-dependent $(l,m)$-coefficients of the wavefunction for
!   the $i$th state are given by
!   $$ \Phi^{i{\bf p}}_{\alpha lm}(r)=\sum_{\bf G}b^{i{\bf p}}_{\bf G}
!    \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})u^{\alpha}_{jl}(r)
!    +\sum_{j=1}^{N^{\alpha}}b^{i{\bf p}}_{(\alpha,j,m)}v^{\alpha}_j(r)
!    \delta_{l,l_j}, $$
!   where $b^{i{\bf p}}$ is the $i$th eigenvector returned from routine
!   {\tt eveqn}; $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient;
!   $M^{\alpha}_l$ is the order of the APW; $u^{\alpha}_{jl}$ is the APW radial
!   function; $N^{\alpha}$ is the number of local-orbitals; $v^{\alpha}_j$ is
!   the $j$th local-orbital radial function; and $(\alpha,j,m)$ is a compound
!   index for the location of the local-orbital in the eigenvector. See routines
!   {\tt genapwfr}, {\tt genlofr}, {\tt match} and {\tt eveqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed description, October 2004 (C. Brouder)
!   Removed argument ist, November 2006 (JKD)
!   Changed arguments and optimised, December 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp,ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: evecfv(nmatmax)
complex(8), intent(out) :: wfmt(*)
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
! zero the wavefunction
wfmt(1:np)=0.d0
!-----------------------!
!     APW functions     !
!-----------------------!
do l=0,lmaxo
  do lm=l**2+1,(l+1)**2
    i=npi+lm
    do io=1,apword(l,is)
      z1=zdotu(ngp,evecfv,1,apwalm(:,io,lm),1)
      if (l.le.lmaxi) then
        call zfzrf(nri,z1,lrstp,apwfr(1,1,io,l,ias),lmmaxi,wfmt(lm))
      end if
      call zfzrf(nro,z1,lrstp,apwfr(iro,1,io,l,ias),lmmaxo,wfmt(i))
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
    z1=evecfv(ngp+idxlo(lm,ilo,ias))
    if (l.le.lmaxi) then
      call zfzrf(nri,z1,lrstp,lofr(1,1,ilo,ias),lmmaxi,wfmt(lm))
    end if
    call zfzrf(nro,z1,lrstp,lofr(iro,1,ilo,ias),lmmaxo,wfmt(i))
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
!EOC

