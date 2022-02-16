!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_input
!**********************************************************************
!                                                                     *
!     Module for data input                                           *
!     (Here will be defined uniform decomposition and padding with    *
!      zeros of control variables, required by the filter)            *
!                                                                     *
!**********************************************************************

public input_2d
public input_spec1_2d
public input_3d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine input_2d                             &
!***********************************************************************
!                                                                      !
!   Define some function for testing redecomposition                   !
!   (for analysis grid)                                                !
!                                                                      !
!***********************************************************************
(V,imin,jmin,imax,jmax,imax0,ampl)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
use mg_mppstuff, only: nx,my
implicit none
integer(i_kind),intent(in):: imax,jmax
integer(i_kind),intent(in):: imin,jmin
integer(i_kind),intent(in):: imax0
integer(i_kind),intent(in):: ampl
real(r_kind),dimension(imin:imax,jmin:jmax),intent(out):: V
real(i_kind):: ng,mg,L,m,n
!-----------------------------------------------------------------------


     do m=imin,jmax
       mg = (my-1)*jmax+m
     do n=jmin,imax
       ng = (nx-1)*imax+n
         V(n,m)=ampl*(mg*imax0+ng)
!         V(n,m)=0.
     enddo
     enddo

!-----------------------------------------------------------------------
                        endsubroutine input_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine input_spec1_2d                       &
!***********************************************************************
!                                                                      !
!   Define some function for testing redecomposition                   !
!   (for analysis grid)                                                !
!                                                                      !
!***********************************************************************
(V,nx0,my0,flag)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
use mg_mppstuff, only: nx,my
use mg_parameter, only: nm,mm
implicit none
integer(i_kind),intent(in):: nx0,my0
real(r_kind),dimension(0:nm,0:mm),intent(out):: V
character(len=2), intent(in):: flag
integer(r_kind):: v0=1.
!-----------------------------------------------------------------------

    V(:,:)=0.

if(flag=='md') then
 if(nx==nx0.and.my==my0) then
    V(nm/2,mm/2)=v0
 endif
else &
if(flag=='rt') then
 if(nx==nx0.and.my==my0) then
    V(nm,mm)=v0
 endif
 if(nx==nx0+1.and.my==my0) then
    V(0,mm)=v0
 endif
 if(nx==nx0.and.my==my0+1) then
    V(nm,0)=v0
 endif
 if(nx==nx0+1.and.my==my0+1) then
    V(0,0)=v0
 endif
endif

!-----------------------------------------------------------------------
                        endsubroutine input_spec1_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine input_3d                             &
!***********************************************************************
!                                                                      !
!   Define some function for testing redecomposition                   !
!   (for analysis grid)                                                !
!                                                                      !
!***********************************************************************
(V,imin,jmin,lmin,imax,jmax,lmax,imax0,ampl,incrm)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
use mg_mppstuff, only: nx,my
implicit none
integer(i_kind),intent(in):: imin,jmin,lmin
integer(i_kind),intent(in):: imax,jmax,lmax
integer(i_kind),intent(in):: imax0
integer(i_kind),intent(in):: ampl,incrm
real(r_kind),dimension(lmin:lmax,imin:imax,jmin:jmax),intent(out):: V
real(i_kind):: ng,mg,L,m,n
!-----------------------------------------------------------------------


   do l=lmin,lmax
     do m=imin,jmax
       mg = (my-1)*jmax+m
     do n=jmin,imax
       ng = (nx-1)*imax+n
         V(l,n,m)=ampl*(mg*imax0+ng) +(l-1)*incrm
!         V(l,n,m)=0.
     enddo
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine input_3d


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_input
