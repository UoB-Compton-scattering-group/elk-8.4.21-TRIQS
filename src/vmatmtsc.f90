
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vmatmtsc
use modmain
use moddftu
implicit none
! local variables
integer ispn,jspn,ias,lm
! automatic arrays
complex(8) a(lmmaxo,lmmaxo),b(lmmaxo,lmmaxo)
lm=min(lmmaxi,lmmaxdm)
do ias=1,natmtot
  if (any(tvmmt(0:lmaxdm,ias))) then
    do jspn=1,nspinor
      do ispn=1,nspinor
! transform the muffin-tin potential matrix elements from the spherical harmonic
! basis to the spherical coordinate basis
        a(1:lmmaxi,1:lmmaxi)=0.d0
        a(1:lm,1:lm)=vmatmt(1:lm,ispn,1:lm,jspn,ias)
        call zgemm('N','N',lmmaxi,lmmaxi,lmmaxi,zone,a,lmmaxo,zfshti,lmmaxi, &
         zzero,b,lmmaxo)
        call zgemm('N','N',lmmaxi,lmmaxi,lmmaxi,zone,zbshti,lmmaxi,b,lmmaxo, &
         zzero,vmatmti(:,:,ispn,jspn,ias),lmmaxi)
        a(:,:)=0.d0
        a(1:lmmaxdm,1:lmmaxdm)=vmatmt(1:lmmaxdm,ispn,1:lmmaxdm,jspn,ias)
        call zgemm('N','N',lmmaxo,lmmaxo,lmmaxo,zone,a,lmmaxo,zfshto,lmmaxo, &
         zzero,b,lmmaxo)
        call zgemm('N','N',lmmaxo,lmmaxo,lmmaxo,zone,zbshto,lmmaxo,b,lmmaxo, &
         zzero,vmatmto(:,:,ispn,jspn,ias),lmmaxo)
      end do
    end do
  end if
end do
end subroutine

