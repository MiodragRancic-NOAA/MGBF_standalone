!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_output
!**********************************************************************
!                                                                     *
!     Module for data output                                          *
!                                                                     *
!**********************************************************************

public output_spec1_2d
public output_spec2_2d
public output_vertical_2d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_spec1_2d                      &
!***********************************************************************
!                                                                      !
!   Outpyt a 2D array                                                  !
!                                                                      !
!***********************************************************************
(V)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
use mg_parameter, only: nm,mm
use mg_mppstuff, only: mype
implicit none
real(r_kind),dimension(0:nm,0:mm),intent(in):: V
real(i_kind):: m,n
!-----------------------------------------------------------------------

    do m=0,mm
      write(100+mype,'(f9.3)') (V(n,m),n=0,nm)
    enddo

!-----------------------------------------------------------------------
                        endsubroutine output_spec1_2d


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_vertical_2d                   &
!***********************************************************************
!                                                                      !
!   Output of verical a 2D array                                       !
!                                                                      !
!***********************************************************************
(V,my0)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
use mg_parameter, only: nm,mm,Lm,mym
use mg_mppstuff, only: mype,nx,my,barrierMPI
implicit none
real(r_kind),dimension(0:nm,1:LM),intent(in):: V
integer(i_kind):: my0
real(i_kind):: n,l
!-----------------------------------------------------------------------

if(my==my0) then
    do L=1,Lm
      write(500+nx,'(f9.3)') (V(n,L),n=0,nm)
    enddo
endif
call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine output_vertical_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_spec2_2d                      &
!***********************************************************************
!                                                                      !
!   Outlyt a 2D array                                                  !
!                                                                      !
!***********************************************************************
(V)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
use mg_parameter, only: imL,jmL
use mg_mppstuff, only: mype
implicit none
real(r_kind),dimension(-1:imL+1,-1:jmL+1),intent(in):: V
real(i_kind):: j,i
!-----------------------------------------------------------------------

    do j=jmL+1,-1,-1
      write(100+mype,'(f9.3)') (V(i,j),i=-1,imL+1)
    enddo

!-----------------------------------------------------------------------
                        endsubroutine output_spec2_2d


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_output
