!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_interpolate
!***********************************************************************
!                                                                      !
!    general mapping between 2d arrays using linerly squared           !
!    interpolations                                                    !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use kinds
use mg_parameter, only: x0,y0,dxa,dxf,dya,dyf                           &
                        ,nm,mm,km,km2,km3,lm,lm_all                     &
                        ,im,jm,ib,jb                                    
use mg_intstate, only: iref,jref                                        &
                      ,cx0,cx1,cx2,cx3                                  &
                      ,cy0,cy1,cy2,cy3                                  &
                      ,cf00,cf01,cf02,cf03                              &
                      ,cf10,cf11,cf12,cf13                              &
                      ,cf20,cf21,cf22,cf23                              &
                      ,cf30,cf31,cf32,cf33
!use mpimod, only: mype
use mg_mppstuff, only: mype
use mg_mppstuff, only: finishMPI
implicit none

public lsqr_mg_coef
public lwq_vertical_coef
public lwq_vertical_direct
public lwq_vertical_adjoint


interface lsqr_forward_xk
  module procedure lsqr_forward_x_k2d
  module procedure lsqr_forward_x_k3d
endinterface

interface lsqr_adjoint_xk
  module procedure lsqr_adjoint_x_k2d
  module procedure lsqr_adjoint_x_k3d
endinterface

interface lsqr_forward_xyk
  module procedure lsqr_forward_xy_k2d
  module procedure lsqr_forward_xy_k3d
endinterface

interface lsqr_adjoint_xyk
  module procedure lsqr_adjoint_xy_k2d
  module procedure lsqr_adjoint_xy_k3d
endinterface

interface lsqr_forward
  module procedure lsqr_forward_k2d
  module procedure lsqr_forward_k3d
endinterface

interface lsqr_adjoint
  module procedure lsqr_adjoint_k2d
  module procedure lsqr_adjoint_k3d
endinterface

public lsqr_forward_all
public lsqr_adjoint_all

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_mg_coef                            
!***********************************************************************
!                                                                      !
!   Prepare coeficients for mapping between:                           !
!        filter grid on analysis decomposition: W(-ib:im+ib,-jb:jm+jb) !
!        and analysis grid:                     V(0:nm,0:mm)           !  
!                                                                      !
!              (  im <= nm  and  jm < mm   )                           !
!                                                                      !
!***********************************************************************
implicit none
real(r_kind), dimension(0:nm):: xa
real(r_kind), dimension(-ib:im+ib):: xf
real(r_kind), dimension(0:mm):: ya
real(r_kind), dimension(-jb:jm+jb):: yf
integer(i_kind):: i,j,n,m
real(r_kind) x1,x2,x3,x4,x
real(r_kind) x1x,x2x,x3x,x4x
real(r_kind) rx2x1,rx3x1,rx4x1,rx3x2,rx4x2,rx4x3
real(r_kind) y1,y2,y3,y4,y
real(r_kind) y1y,y2y,y3y,y4y
real(r_kind) ry2y1,ry3y1,ry4y1,ry3y2,ry4y2,ry4y3
real(r_kind) cfl1,cfl2,cfl3,cll
real(r_kind) cfr1,cfr2,cfr3,crr
!-----------------------------------------------------------------------
!
! Initialize
!
 
   do n=0,nm
     xa(n)=x0+n*dxa
   enddo

   do i=-ib,im+ib
     xf(i)=x0+i*dxf
   enddo

   do m=0,mm
     ya(m)=y0+m*dya
   enddo

   do j=-jb,jm+jb
     yf(j)=y0+j*dyf
   enddo

!
! Find iref and jref
!
   do n=0,nm
     do i=-ib,im+ib-1
       if( xa(n)< xf(i)) then
         iref(n)=i-2
         exit
       endif
     enddo
   enddo

   do m=0,mm
     do j=-jb,jm+jb-1
       if(ya(m) < yf(j)) then
         jref(m)=j-2
         exit
       endif
     enddo
   enddo


   do n=0,nm
     i=iref(n)
     x1=xf(i)
     x2=xf(i+1)
     x3=xf(i+2)
     x4=xf(i+3)
     x = xa(n)
       x1x = x1-x   
       x2x = x2-x   
       x3x = x3-x   
       x4x = x4-x   
       rx2x1 = 1./(x2-x1)
       rx3x1 = 1./(x3-x1)
       rx4x1 = 1./(x4-x1)
       rx3x2 = 1./(x3-x2)
       rx4x2 = 1./(x4-x2)
       rx4x3 = 1./(x4-x3)
     CFL1 = x2x*x3x*rx2x1*rx3x1
     CFL2 =-x1x*x3x*rx2x1*rx3x2
     CFL3 = x1x*x2x*rx3x1*rx3x2
     CLL = x3x*rx3x2
     CFR1 = x3x*x4x*rx3x2*rx4x2
     CFR2 =-x2x*x4x*rx3x2*rx4x3
     CFR3 = x2x*x3x*rx4x2*rx4x3
     CRR =-x2x*rx3x2
       cx0(n)=CFL1*CLL
       cx1(n)=CFL2*CLL+CFR1*CRR
       cx2(n)=CFL3*CLL+CFR2*CRR
       cx3(n)=CFR3*CRR
   enddo

   do m=0,mm
     j=jref(m)
     y1=yf(j)
     y2=yf(j+1)
     y3=yf(j+2)
     y4=yf(j+3)
     y = ya(m)
       y1y = y1-y   
       y2y = y2-y   
       y3y = y3-y   
       y4y = y4-y   
       ry2y1 = 1./(y2-y1)
       ry3y1 = 1./(y3-y1)
       ry4y1 = 1./(y4-y1)
       ry3y2 = 1./(y3-y2)
       ry4y2 = 1./(y4-y2)
       ry4y3 = 1./(y4-y3)
     CFL1 = y2y*y3y*ry2y1*ry3y1
     CFL2 =-y1y*y3y*ry2y1*ry3y2
     CFL3 = y1y*y2y*ry3y1*ry3y2
     CLL = y3y*ry3y2
     CFR1 = y3y*y4y*ry3y2*ry4y2
     CFR2 =-y2y*y4y*ry3y2*ry4y3
     CFR3 = y2y*y3y*ry4y2*ry4y3
     CRR =-y2y*ry3y2
       cy0(m)=CFL1*CLL
       cy1(m)=CFL2*CLL+CFR1*CRR
       cy2(m)=CFL3*CLL+CFR2*CRR
       cy3(m)=CFR3*CRR
   enddo

   do m=0,mm
   do n=0,nm
     cf00(n,m)=cx0(n)*cy0(m)
     cf01(n,m)=cx0(n)*cy1(m)
     cf02(n,m)=cx0(n)*cy2(m)
     cf03(n,m)=cx0(n)*cy3(m)
     cf10(n,m)=cx1(n)*cy0(m)
     cf11(n,m)=cx1(n)*cy1(m)
     cf12(n,m)=cx1(n)*cy2(m)
     cf13(n,m)=cx1(n)*cy3(m)
     cf20(n,m)=cx2(n)*cy0(m)
     cf21(n,m)=cx2(n)*cy1(m)
     cf22(n,m)=cx2(n)*cy2(m)
     cf23(n,m)=cx2(n)*cy3(m)
     cf30(n,m)=cx3(n)*cy0(m)
     cf31(n,m)=cx3(n)*cy1(m)
     cf32(n,m)=cx3(n)*cy2(m)
     cf33(n,m)=cx3(n)*cy3(m)
   enddo
   enddo
 
!-----------------------------------------------------------------------
                        endsubroutine lsqr_mg_coef

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_coef                    &
!***********************************************************************
!                                                                      !
!  Prepare coeficients for vetical mapping between:                    !
!  analysis grid vertical resolution (nm) and                          !
!  generation one of filter grid vertical resoluition (im)             !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(nm,im,c1,c2,c3,c4,iref)
use mg_mppstuff, only: mype
implicit none

integer(i_kind), intent(in):: nm,im
real(r_kind), dimension(1:nm), intent(out):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(out):: iref

real(r_kind), dimension(1:nm):: y
real(r_kind), dimension(0:im+1):: x
real(r_kind):: dy,x1,x2,x3,x4,dx1,dx2,dx3,dx4 
real(r_kind):: dx13,dx23,dx24

integer(i_kind):: i,n
!-----------------------------------------------------------------------

   do i=0,im+1
     x(i)=(i-1)*1.
   enddo

    dy = 1.*(im-1)/(nm-1)
  do n=1,nm
    y(n)=(n-1)*dy
  enddo
    y(nm)=x(im)
 
  do n=2,nm-1
    i = y(n)+1
      x1 = x(i-1)
      x2 = x(i)
      x3 = x(i+1)
      x4 = x(i+2)
    iref(n)=i
      dx1 = y(n)-x1
      dx2 = y(n)-x2
      dx3 = y(n)-x3
      dx4 = y(n)-x4
      dx13 = dx1*dx3
      dx23 = 0.5*dx2*dx3
      dx24 = dx2*dx4
    c1(n) = -dx23*dx3
    c2(n) =  (    dx13+0.5*dx24)*dx3
    c3(n) = -(0.5*dx13+    dx24)*dx2
    c4(n) = dx23*dx2

    if(iref(n)==1) then
      c3(n)=c3(n)+c1(n)
      c1(n)=0.
    endif
    if(iref(n)==im-1) then
      c2(n)=c2(n)+c4(n)
      c4(n)=0.
    endif
  enddo
     iref(1)=1; c1(1)=0.; c2(1)=1.; c3(1)=0.; c4(1)=0.
     iref(nm)=im; c1(nm)=0.; c2(nm)=1.; c3(nm)=0.; c4(n)=0.
     
!TEST
!if(mype==0) then 
!     print *, '---------------------'
!     print *,'nm,im=',nm,im
!     print *,'y(n),iref(n),n'
!  do n=1,nm
!     print *,y(n),iref(n),n
!  enddo 
!     print *,' '
!     print *,'x(i),i'
!  do i=0,im+1
!    print *,x(i),i
!  enddo
!     print *, '---------------------'
!     print *,'nm,im=',nm,im
!     print *,'iref(n),c1(n),c2(n),c3(n)'
!  do n=1,nm
!     print *,iref(n),c1(n),c2(n),c3(n)
!  enddo 
!endif
!TEST

       
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_coef                            

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_adjoint                 &
!***********************************************************************
!                                                                      !
!  Direct linerly weighted quadratic adjoint interpolation in vertical !
!  from reslution nm to resolution im                                  !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(nm,km,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,w,f)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: nm,km,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:nm,imin:imax,jmin:jmax), intent(in):: w
real(r_kind), dimension(1:km,imin:imax,jmin:jmax), intent(out):: f
integer(i_kind):: k,n
!-----------------------------------------------------------------------
  f = 0.
do n=2,nm-1
  k = kref(n)
  if( k==1 ) then
    f(1,:,:) = f(1,:,:)+c2(n)*w(n,:,:)
    f(2,:,:) = f(2,:,:)+c3(n)*w(n,:,:)
    f(3,:,:) = f(3,:,:)+c4(n)*w(n,:,:)
  elseif &
    ( k==km-1) then
    f(km-2,:,:) = f(km-2,:,:)+c1(n)*w(n,:,:)
    f(km-1,:,:) = f(km-1,:,:)+c2(n)*w(n,:,:)
    f(km  ,:,:) = f(km  ,:,:)+c3(n)*w(n,:,:)
  elseif( k==km) then
    f(k  ,:,:) = f(k  ,:,:)+c2(n)*w(n,:,:)
  else
    f(k-1,:,:) = f(k-1,:,:)+c1(n)*w(n,:,:)
    f(k  ,:,:) = f(k  ,:,:)+c2(n)*w(n,:,:)
    f(k+1,:,:) = f(k+1,:,:)+c3(n)*w(n,:,:)
    f(k+2,:,:) = f(k+2,:,:)+c4(n)*w(n,:,:)
  endif
enddo
    f(1,:,:)=f(1,:,:)+w(1,:,:)
    f(km,:,:)=f(km,:,:)+w(nm,:,:)

!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_adjoint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_direct                  &
!***********************************************************************
!                                                                      !
!  Linerly weighted direct quadratic interpolation in vertical         !
!  from reslouion im to resolution nm                                  !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(km,nm,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,f,w)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,nm,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:km,imin:imax,jmin:jmax), intent(in):: f
real(r_kind), dimension(1:nm,imin:imax,jmin:jmax), intent(out):: w
integer(i_kind):: k,n
!-----------------------------------------------------------------------
do n=2,nm-1
  k = kref(n)
  if( k==1 ) then
    w(n,:,:) =             c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  elseif &
    ( k==km-1) then
    w(n,:,:) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)
  elseif &
    ( k==km)   then
    w(n,:,:) =                 c2(n)*f(k,:,:)
  else
    w(n,:,:) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,: )+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  endif
enddo
    w(1,:,:)=f(1,:,:)
    w(nm,:,:)=f(km,:,:)
    
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_direct

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_k2d                     &
!***********************************************************************
!                                                                      !
! Given a source array  V(-ib:im+ib,-jb:jm+jb,km2) perform             !
! forward  interpolations to get target array W(0:nm,0,mm,km2)         ! 
!                                                                      !
!***********************************************************************
(V,W)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-ib:im+ib,-ib:im+ib,km2), intent(in):: V
real(r_kind), dimension(0:nm,0:mm,km2),intent(out):: W  
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km2):: v00,v01,v02,v03                            &
                             ,v10,v11,v12,v13                            &
                             ,v20,v21,v22,v23                            &
                             ,v30,v31,v32,v33                                            
!-----------------------------------------------------------------------

   do m=0,mm
     j = jref(m)
   do n=0,nm
     i = iref(n)
     v00(:)=V(i  ,j  ,:)
     v10(:)=V(i+1,j  ,:)
     v20(:)=V(i+2,j  ,:)
     v30(:)=V(i+3,j  ,:)
     v01(:)=V(i  ,j+1,:)
     v11(:)=V(i+1,j+1,:)
     v21(:)=V(i+2,j+1,:)
     v31(:)=V(i+3,j+1,:)
     v02(:)=V(i  ,j+2,:)
     v12(:)=V(i+1,j+2,:)
     v22(:)=V(i+2,j+2,:)
     v32(:)=V(i+3,j+2,:)
     v03(:)=V(i  ,j+3,:)
     v13(:)=V(i+1,j+3,:)
     v23(:)=V(i+2,j+3,:)
     v33(:)=V(i+3,j+3,:)
     W(n,m,:) = cf00(n,m)*v00(:)+cf10(n,m)*v10(:)+cf20(n,m)*v20(:)+cf30(n,m)*v30(:)   &
               +cf01(n,m)*v01(:)+cf11(n,m)*v11(:)+cf21(n,m)*v21(:)+cf31(n,m)*v31(:)   &
               +cf02(n,m)*v02(:)+cf12(n,m)*v12(:)+cf22(n,m)*v22(:)+cf32(n,m)*v32(:)   &
               +cf03(n,m)*v03(:)+cf13(n,m)*v13(:)+cf23(n,m)*v23(:)+cf33(n,m)*v33(:) 
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_k3d                     &
!***********************************************************************
!                                                                      !
! Given a source array  V(-ib:im+ib,-jb:jm+jb,lm,km3) perform forward  !
! interpolations to get target array W(0:nm,0,mm,lm,km3)               ! 
!                                                                      !
!***********************************************************************
(V,W)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,lm,km3), intent(in):: V
real(r_kind), dimension(0:nm,0:mm,lm,km3),intent(out):: W  
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km3):: v00,v01,v02,v03                             &
                             ,v10,v11,v12,v13                             &
                             ,v20,v21,v22,v23                             &
                             ,v30,v31,v32,v33                                            
!-----------------------------------------------------------------------

   do l=1,lm
   do m=0,mm
     j = jref(m)
   do n=0,nm
     i = iref(n)
     v00(:)=V(i  ,j  ,l,:)
     v10(:)=V(i+1,j  ,l,:)
     v20(:)=V(i+2,j  ,l,:)
     v30(:)=V(i+3,j  ,l,:)
     v01(:)=V(i  ,j+1,l,:)
     v11(:)=V(i+1,j+1,l,:)
     v21(:)=V(i+2,j+1,l,:)
     v31(:)=V(i+3,j+1,l,:)
     v02(:)=V(i  ,j+2,l,:)
     v12(:)=V(i+1,j+2,l,:)
     v22(:)=V(i+2,j+2,l,:)
     v32(:)=V(i+3,j+2,l,:)
     v03(:)=V(i  ,j+3,l,:)
     v13(:)=V(i+1,j+3,l,:)
     v23(:)=V(i+2,j+3,l,:)
     v33(:)=V(i+3,j+3,l,:)
     W(n,m,l,:) = cf00(n,m)*v00(:)+cf10(n,m)*v10(:)+cf20(n,m)*v20(:)+cf30(n,m)*v30(:)   &
                 +cf01(n,m)*v01(:)+cf11(n,m)*v11(:)+cf21(n,m)*v21(:)+cf31(n,m)*v31(:)   &
                 +cf02(n,m)*v02(:)+cf12(n,m)*v12(:)+cf22(n,m)*v22(:)+cf32(n,m)*v32(:)   &
                 +cf03(n,m)*v03(:)+cf13(n,m)*v13(:)+cf23(n,m)*v23(:)+cf33(n,m)*v33(:) 
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_k2d                      &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,km2) perform adjoint interpolations !
! in order to get source array V(-ib:im+ib,-jb:jm+jb,km2)              ! 
!                                                                      !
!***********************************************************************
(W,V)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:nm,0:mm,km2),intent(in):: W  
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,km2), intent(out):: V
integer(i_kind):: i,j,n,m
integer(i_kind):: ip1,ip2,ip3,jp1,jp2,jp3
!-----------------------------------------------------------------------
     
    V(:,:,:) = 0.

!$OMP PARALLEL DO 
   do m=0,mm
       j = jref(m)
       jp1=j+1
       jp2=j+2
       jp3=j+3
   do n=0,nm
       i = iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
     V(i  ,j  ,:) = V(i  ,j  ,:)+W(n,m,:)*cf00(n,m)
     V(ip1,j  ,:) = V(ip1,j  ,:)+W(n,m,:)*cf10(n,m)
     V(ip2,j  ,:) = V(ip2,j  ,:)+W(n,m,:)*cf20(n,m)
     V(ip3,j  ,:) = V(ip3,j  ,:)+W(n,m,:)*cf30(n,m)
     V(i  ,jp1,:) = V(i  ,jp1,:)+W(n,m,:)*cf01(n,m) 
     V(ip1,jp1,:) = V(ip1,jp1,:)+W(n,m,:)*cf11(n,m)
     V(ip2,jp1,:) = V(ip2,jp1,:)+W(n,m,:)*cf21(n,m)
     V(ip3,jp1,:) = V(ip3,jp1,:)+W(n,m,:)*cf31(n,m)
     V(i  ,jp2,:) = V(i  ,jp2,:)+W(n,m,:)*cf02(n,m)
     V(ip1,jp2,:) = V(ip1,jp2,:)+W(n,m,:)*cf12(n,m)
     V(ip2,jp2,:) = V(ip2,jp2,:)+W(n,m,:)*cf22(n,m)
     V(ip3,jp2,:) = V(ip3,jp2,:)+W(n,m,:)*cf32(n,m)
     V(i  ,jp3,:) = V(i  ,jp3,:)+W(n,m,:)*cf03(n,m)
     V(ip1,jp3,:) = V(ip1,jp3,:)+W(n,m,:)*cf13(n,m)
     V(ip2,jp3,:) = V(ip2,jp3,:)+W(n,m,:)*cf23(n,m)
     V(ip3,jp3,:) = V(ip3,jp3,:)+W(n,m,:)*cf33(n,m)
   enddo
   enddo

!$OMP END PARALLEL DO 
!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_k2d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_k3d                      &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,lm,km3) perform adjoint             !
! interpolations in order to get source V(-ib:im+ib,-jb:jm+jb,lm,km3)  ! 
!                                                                      !
!***********************************************************************
(W,V)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:nm,0:mm,lm,km3),intent(in):: W  
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,lm,km3), intent(out):: V
integer(i_kind):: i,j,n,m,l
integer(i_kind):: ip1,ip2,ip3,jp1,jp2,jp3
!-----------------------------------------------------------------------
     
    V(:,:,:,:) = 0.

   do l=1,lm
   do m=0,mm
       j = jref(m)
       jp1=j+1
       jp2=j+2
       jp3=j+3
   do n=0,nm
       i = iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
     V(i  ,j  ,l,:) = V(i  ,j  ,l,:)+W(n,m,l,:)*cf00(n,m)
     V(ip1,j  ,l,:) = V(ip1,j  ,l,:)+W(n,m,l,:)*cf10(n,m)
     V(ip2,j  ,l,:) = V(ip2,j  ,l,:)+W(n,m,l,:)*cf20(n,m)
     V(ip3,j  ,l,:) = V(ip3,j  ,l,:)+W(n,m,l,:)*cf30(n,m)
     V(i  ,jp1,l,:) = V(i  ,jp1,l,:)+W(n,m,l,:)*cf01(n,m) 
     V(ip1,jp1,l,:) = V(ip1,jp1,l,:)+W(n,m,l,:)*cf11(n,m)
     V(ip2,jp1,l,:) = V(ip2,jp1,l,:)+W(n,m,l,:)*cf21(n,m)
     V(ip3,jp1,l,:) = V(ip3,jp1,l,:)+W(n,m,l,:)*cf31(n,m)
     V(i  ,jp2,l,:) = V(i  ,jp2,l,:)+W(n,m,l,:)*cf02(n,m)
     V(ip1,jp2,l,:) = V(ip1,jp2,l,:)+W(n,m,l,:)*cf12(n,m)
     V(ip2,jp2,l,:) = V(ip2,jp2,l,:)+W(n,m,l,:)*cf22(n,m)
     V(ip3,jp2,l,:) = V(ip3,jp2,l,:)+W(n,m,l,:)*cf32(n,m)
     V(i  ,jp3,l,:) = V(i  ,jp3,l,:)+W(n,m,l,:)*cf03(n,m)
     V(ip1,jp3,l,:) = V(ip1,jp3,l,:)+W(n,m,l,:)*cf13(n,m)
     V(ip2,jp3,l,:) = V(ip2,jp3,l,:)+W(n,m,l,:)*cf23(n,m)
     V(ip3,jp3,l,:) = V(ip3,jp3,l,:)+W(n,m,l,:)*cf33(n,m)
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_k3d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_x_k2d                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(-ib:im+ib,0:mm,km2) perform forward          !
! interpolations to get target array W(0:nm,0:mm,km2)                  ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(V,W)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-ib:im+ib,0:mm,km2), intent(in):: V
real(r_kind), dimension(0:nm,0:mm,km2),intent(out):: W  
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km2):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

   do m=0,mm
       j = m
   do n=0,nm
       i = iref(n)
     v0(:)=V(i  ,j  ,:)
     v1(:)=V(i+1,j  ,:)
     v2(:)=V(i+2,j  ,:)
     v3(:)=V(i+3,j  ,:)
     W(n,m,:) = cx0(n)*v0(:)+cx1(n)*v1(:)+cx2(n)*v2(:)+cx3(n)*v3(:)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_x_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_x_k3d                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(-ib:im+ib,-jb:jm+jb,lm,km3) perform forward  !
! interpolations to get target array W(0:nm,0:mm,lm,km3)               ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(V,W)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-ib:im+ib,0:mm,lm,km3), intent(in):: V
real(r_kind), dimension(0:nm,0:mm,lm,km3),intent(out):: W  
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km3):: v0,v1,v2,v3 
!-----------------------------------------------------------------------

   do l=1,lm
   do m=0,mm
       j = m
   do n=0,nm
       i = iref(n)
     v0(:)=V(i  ,j  ,l,:)
     v1(:)=V(i+1,j  ,l,:)
     v2(:)=V(i+2,j  ,l,:)
     v3(:)=V(i+3,j  ,l,:)
     W(n,m,l,:) = cx0(n)*v0(:)+cx1(n)*v1(:)+cx2(n)*v2(:)+cx3(n)*v3(:)  
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_x_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_x_k2d                   &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,km2) perform adjoint                !
! interpolations in order to get source V(-ib:im+ib,0:mm,km2)          ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(W,V)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:nm,0:mm,km2),intent(in):: W  
real(r_kind), dimension(-ib:im+ib,0:mm,km2), intent(out):: V
integer(i_kind):: i,j,n,m
integer(i_kind):: ip1,ip2,ip3
!-----------------------------------------------------------------------
     
    V(:,:,:) = 0.

   do m=0,mm
       j = m
   do n=0,nm
       i = iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
       j = m
     V(i  ,j  ,:) = V(i  ,j  ,:)+W(n,m,:)*cx0(n)
     V(ip1,j  ,:) = V(ip1,j  ,:)+W(n,m,:)*cx1(n)
     V(ip2,j  ,:) = V(ip2,j  ,:)+W(n,m,:)*cx2(n)
     V(ip3,j  ,:) = V(ip3,j  ,:)+W(n,m,:)*cx3(n)
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_x_k2d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_x_k3d                   &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,lm,km3) perform adjoint             !
! interpolations in order to get source V(-ib:im+ib,0:mm,lm,km3)       ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(W,V)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:nm,0:mm,lm,km3),intent(in):: W  
real(r_kind), dimension(-ib:im+ib,0:mm,lm,km3), intent(out):: V
integer(i_kind):: i,j,n,m,l
integer(i_kind):: ip1,ip2,ip3
!-----------------------------------------------------------------------
     
    V(:,:,:,:) = 0.

   do l=1,lm
   do m=0,mm
       j = m
   do n=0,nm
       i = iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
       j = m
     V(i  ,j  ,l,:) = V(i  ,j  ,l,:)+W(n,m,l,:)*cx0(n)
     V(ip1,j  ,l,:) = V(ip1,j  ,l,:)+W(n,m,l,:)*cx1(n)
     V(ip2,j  ,l,:) = V(ip2,j  ,l,:)+W(n,m,l,:)*cx2(n)
     V(ip3,j  ,l,:) = V(ip3,j  ,l,:)+W(n,m,l,:)*cx3(n)
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_x_k3d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_xy_k2d                  &
!***********************************************************************
!                                                                      !
!  Given a source array  V(-ib:im+ib,-jb:jm+jb,km2) perform forward    !
!  interpolations to get target array W(0:nm,0:mm,km2)                 ! 
!  using two passes of 1d interpolator                                 !
!                                                                      !
!***********************************************************************
(V,W)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,km2), intent(in):: V
real(r_kind), dimension(0:nm,0:mm,km2),intent(out):: W  

real(r_kind), dimension(0:nm,-jb:jm+jb,km2):: VX
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km2):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

   do j=-jb,jm+jb
   do n=0,nm
       i = iref(n)
     v0(:)=V(i  ,j  ,:)
     v1(:)=V(i+1,j  ,:)
     v2(:)=V(i+2,j  ,:)
     v3(:)=V(i+3,j  ,:)
     VX(n,j,:) = cx0(n)*v0(:)+cx1(n)*v1(:)+cx2(n)*v2(:)+cx3(n)*v3(:)
   enddo
   enddo

   do m=0,mm
     j = jref(m)
   do n=0,nm
     v0(:)=VX(n,j  ,:) 
     v1(:)=VX(n,j+1,:) 
     v2(:)=VX(n,j+2,:) 
     v3(:)=VX(n,j+3,:) 
     W(n,m,:) =  cy0(m)*v0(:)+cy1(m)*v1(:)+cy2(m)*v2(:)+cy3(m)*v3(:)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_xy_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_xy_k3d                  &
!***********************************************************************
!                                                                      !
! Given a source array  V(-ib:im+ib,-jb:jm+jb,lm,km3) perform forward  !
! interpolations to get target array W(0:nm,0:mm,lm,km3)               ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(V,W)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,lm,km3), intent(in):: V
real(r_kind), dimension(0:nm,0:mm,lm,km3),intent(out):: W  

real(r_kind), dimension(0:nm,-jb:jm+jb,km3):: VX
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km3):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

 do l=1,lm

   do j=-jb,jm+jb
   do n=0,nm
       i = iref(n)
     v0(:)=V(i  ,j  ,l,:)
     v1(:)=V(i+1,j  ,l,:)
     v2(:)=V(i+2,j  ,l,:)
     v3(:)=V(i+3,j  ,l,:)
     VX(n,j,:) = cx0(n)*v0(:)+cx1(n)*v1(:)+cx2(n)*v2(:)+cx3(n)*v3(:)
   enddo
   enddo

   do m=0,mm
     j = jref(m)
   do n=0,nm
     v0(:)=VX(n,j  ,:) 
     v1(:)=VX(n,j+1,:) 
     v2(:)=VX(n,j+2,:) 
     v3(:)=VX(n,j+3,:) 
     W(n,m,l,:) =  cy0(m)*v0(:)+cy1(m)*v1(:)+cy2(m)*v2(:)+cy3(m)*v3(:)
   enddo
   enddo

 enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_xy_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_all                     &
!***********************************************************************
!                                                                      !
! Given a source array  V(-ib:im+ib,-jb:jm+jb,lm_all) perform forward  !
! interpolations to get target array W(0:nm,0:mm,lm_all)               ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(V,W,km_all)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: km_all
real(r_kind), dimension(km_all,-ib:im+ib,-jb:jm+jb), intent(in):: V
real(r_kind), dimension(km_all,0:nm,0:mm),intent(out):: W  

real(r_kind), dimension(km_all,0:nm,-jb:jm+jb):: VX
integer(i_kind):: i,j,n,m,l
real(r_kind),dimension(km_all):: v0,v1,v2,v3     
!-----------------------------------------------------------------------


   do j=-jb,jm+jb
   do n=0,nm
       i = iref(n)
     v0(:)=V(:,i  ,j)
     v1(:)=V(:,i+1,j)
     v2(:)=V(:,i+2,j)
     v3(:)=V(:,i+3,j)
     VX(:,n,j) = cx0(n)*v0(:)+cx1(n)*v1(:)+cx2(n)*v2(:)+cx3(n)*v3(:)
   enddo
   enddo

   do m=0,mm
     j = jref(m)
   do n=0,nm
     v0(:)=VX(:,n,j  ) 
     v1(:)=VX(:,n,j+1) 
     v2(:)=VX(:,n,j+2) 
     v3(:)=VX(:,n,j+3) 
     W(:,n,m) =  cy0(m)*v0(:)+cy1(m)*v1(:)+cy2(m)*v2(:)+cy3(m)*v3(:)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_xy_k2d                  &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,km2) perform adjoint                !
! interpolations to get source array V(-ib:im+ib,-jb:jm+jb,km2)        ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(W,V)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:nm,0:mm,km2),intent(in):: W  
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,km2), intent(out):: V
real(r_kind), dimension(0:nm,-jb:jm+jb,km2):: VX
integer(i_kind):: i,j,n,m,l
integer(i_kind):: ip1,ip2,ip3
integer(i_kind):: jp1,jp2,jp3
!-----------------------------------------------------------------------
   
   VX(:,:,:)=0.

   do m=0,mm
     j = jref(m)
     jp1=j+1
     jp2=j+2
     jp3=j+3
   do n=0,nm
     VX(n,j  ,:) = VX(n,j  ,:)+W(n,m,:)*cy0(m)
     VX(n,jp1,:) = VX(n,jp1,:)+W(n,m,:)*cy1(m)
     VX(n,jp2,:) = VX(n,jp2,:)+W(n,m,:)*cy2(m)
     VX(n,jp3,:) = VX(n,jp3,:)+W(n,m,:)*cy3(m)
   enddo
   enddo
 
   V(:,:,:) = 0.

   do j=-jb,jm+jb
   do n=0,nm
     i = iref(n)
     ip1=i+1
     ip2=i+2
     ip3=i+3

     V(i  ,j,:) = V(i  ,j,:)+VX(n,j,:)*cx0(n)
     V(ip1,j,:) = V(ip1,j,:)+VX(n,j,:)*cx1(n)
     V(ip2,j,:) = V(ip2,j,:)+VX(n,j,:)*cx2(n)
     V(ip3,j,:) = V(ip3,j,:)+VX(n,j,:)*cx3(n)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_xy_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_xy_k3d                  &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,lm,km3) perform adjoint             !
! interpolations to get source array V(-ib:im+ib,-jb:jm+jb,lm,km3)     ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(W,V)
!-----------------------------------------------------------------------
implicit none
real(r_kind), dimension(0:nm,0:mm,lm,km3),intent(in):: W  
real(r_kind), dimension(-ib:im+ib,-jb:jm+jb,lm,km3), intent(out):: V
real(r_kind), dimension(0:nm,-jb:jm+jb,km3):: VX
integer(i_kind):: i,j,n,m,l,k
integer(i_kind):: ip1,ip2,ip3
integer(i_kind):: jp1,jp2,jp3
!-----------------------------------------------------------------------

   V(:,:,:,:) = 0.
   
 do l=1,lm

   VX(:,:,:)=0.

 do k=1,km3
   do m=0,mm
     j = jref(m)
     jp1=j+1
     jp2=j+2
     jp3=j+3
   do n=0,nm
     VX(n,j  ,k) = VX(n,j  ,k)+W(n,m,l,k)*cy0(m)
     VX(n,jp1,k) = VX(n,jp1,k)+W(n,m,l,k)*cy1(m)
     VX(n,jp2,k) = VX(n,jp2,k)+W(n,m,l,k)*cy2(m)
     VX(n,jp3,k) = VX(n,jp3,k)+W(n,m,l,k)*cy3(m)
   enddo
   enddo
 end do
 

 do k=1,km3
   do j=-jb,jm+jb
   do n=0,nm
     i = iref(n)
     ip1=i+1
     ip2=i+2
     ip3=i+3

     V(i  ,j,l,k) = V(i  ,j,l,k)+VX(n,j,k)*cx0(n)
     V(ip1,j,l,k) = V(ip1,j,l,k)+VX(n,j,k)*cx1(n)
     V(ip2,j,l,k) = V(ip2,j,l,k)+VX(n,j,k)*cx2(n)
     V(ip3,j,l,k) = V(ip3,j,l,k)+VX(n,j,k)*cx3(n)
   enddo
   enddo
 enddo

 enddo
!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_xy_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_all                     &
!***********************************************************************
!                                                                      !
! Given a target array W(0:nm,0:mm,lm,km3) perform adjoint             !
! interpolations to get source array V(-ib:im+ib,-jb:jm+jb,lm,km3)     ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(W,V,km_all)
!-----------------------------------------------------------------------
implicit none
integer(i_kind):: km_all
real(r_kind), dimension(km_all,0:nm,0:mm),intent(in):: W  
real(r_kind), dimension(km_all,-ib:im+ib,-jb:jm+jb), intent(out):: V
real(r_kind), dimension(km_all,0:nm,-jb:jm+jb):: VX
integer(i_kind):: i,j,n,m,l,k
integer(i_kind):: ip1,ip2,ip3
integer(i_kind):: jp1,jp2,jp3
!-----------------------------------------------------------------------

   V(:,:,:) = 0.

   VX(:,:,:)=0.

   do m=0,mm
     j = jref(m)
     jp1=j+1
     jp2=j+2
     jp3=j+3
   do n=0,nm
     VX(:,n,j  ) = VX(:,n,j  )+W(:,n,m)*cy0(m)
     VX(:,n,jp1) = VX(:,n,jp1)+W(:,n,m)*cy1(m)
     VX(:,n,jp2) = VX(:,n,jp2)+W(:,n,m)*cy2(m)
     VX(:,n,jp3) = VX(:,n,jp3)+W(:,n,m)*cy3(m)
   enddo
   enddo
 

   do j=-jb,jm+jb
   do n=0,nm
     i = iref(n)
     ip1=i+1
     ip2=i+2
     ip3=i+3

     V(:,i  ,j) = V(:,i  ,j)+VX(:,n,j)*cx0(n)
     V(:,ip1,j) = V(:,ip1,j)+VX(:,n,j)*cx1(n)
     V(:,ip2,j) = V(:,ip2,j)+VX(:,n,j)*cx2(n)
     V(:,ip3,j) = V(:,ip3,j)+VX(:,n,j)*cx3(n)
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_all  


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_interpolate
