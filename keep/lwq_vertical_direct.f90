!
! From analysis to generation one:
!
!  call lwq_vertical_direct(lm,lmf,cvf1,cvf2,cvf3,cvf4,lref ,farr,warr)
!
! From generation one to generation two
!
!  call lwq_vertical_coef(lmf,lmh,cvh1,cvh2,cvh3,cvh4,lref_h,farr,warr)
!
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
(km,nm,imin,imax,jmin,jmax,c1,c1,c3,c4,kref,f,w)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,nm
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:km,imin:imax,jmin:jmax), intent(in):: f
real(r_kind), dimension(1:nm,imin:imax,jmin:jmax), intent(out):: w
integer(i_kind):: k,n
!-----------------------------------------------------------------------
do n=1,nm
  k = kref(n)
  if( k==1 ) then
    w(n,:,:) =             c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  elseif &
    ( k==km-1) then
    w(n) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)
  else
    w(n) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,: )+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  endif
enddo
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_direct
