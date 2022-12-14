
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fyukawa0
! !INTERFACE:
real(8) function fyukawa0(is,l,k)
! !USES:
use modmain
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   is : species type (in,integer)
!   l  : an angular momentum (in,integer)
!   k  : order of Slater parameter (in,integer)
! !DESCRIPTION:
!   Calculates the Slater parameters in the unscreened case. See {\it Phys. Rev.
!   B} {\bf 52}, 1421 (1995) and {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (LN)
!   Modified and tested July 2008 (LN and FC)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is,l,k
! local variables
integer, parameter :: nstart=1
integer ias,nr,ir
integer ir1,ir2,nr1,nr2
real(8) r2,x
! automatic arrays
real(8) bint(nrmtmax),cint(nrmtmax),fint(nrmtmax)
real(8) blow(nrmtmax),bhigh(nrmtmax)
real(8) clow(nrmtmax),chigh(nrmtmax)
real(8) ak(nrmtmax),bk(nrmtmax)
ias=idxas(1,is)
nr=nrmt(is)
ak(1:nr)=0.d0
bk(1:nr)=0.d0
! calculate unscreened Slater parameters
do ir=1,nr
  r2=rlmt(ir,2,is)
  bint(ir)=fdufr(ir,l,ias)*fdufr(ir,l,ias)*r2
  x=rsp(ir,is)**k
  ak(ir)=x
  bk(ir)=1.d0/(x*rsp(ir,is))
end do
do ir=nstart,nr
  nr1=ir-nstart+1
  nr2=nr-ir+1
  do ir1=1,nr1
    ir2=ir1+nstart-1
    blow(ir1)=bint(ir2)*ak(ir2)
  end do
  call fderiv(-1,nr1,rsp(nstart,is),blow,clow)
  do ir1=1,nr2
    ir2=ir1+ir-1
    bhigh(ir1)=bint(ir2)*bk(ir2)
  end do
  call fderiv(-1,nr2,rsp(ir,is),bhigh,chigh)
  cint(ir-nstart+1)=bint(ir)*(bk(ir)*clow(nr1)+ak(ir)*chigh(nr2))
end do
nr1=nr-nstart+1
call fderiv(-1,nr1,rsp(nstart,is),cint,fint)
fyukawa0=fint(nr1)
end function
!EOC

