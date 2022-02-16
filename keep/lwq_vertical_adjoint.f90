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
(nm,km,imin,imax,jmin,jmax,c1,c1,c3,c4,iref,w,f)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: nm,km,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:nm,imin:imax,jmin:jmax), intent(in):: w
real(r_kind), dimension(1:km,imin:imax,jmin:jmax), intent(out):: f
integer(i_kind):: k,n
!-----------------------------------------------------------------------
do n=1,nm
  k = kref(n)
  if( k==1 ) then
    f(1,:,:) = f(1,:,:)+cf2(n)*w(n,:,:)
    f(2,:,:) = f(2,:,:)+cf3(n)*w(n,:,:)
    f(3,:,:) = f(3,:,:)+cf4(n)*w(n,:,:)
  elseif &
    ( k==km-1) then
    f(km-2,:,:) = f(km-2,:,:)+cf1(n)*w(n,:,:)
    f(km-1,:,:) = f(km-1,:,:)+cf2(n)*w(n,:,:)
    f(km  ,:,:) = f(km  ,:,:)+cf3(n)*w(n,:,:)
  else
    f(k-2,:,:) = f(k-2,:,:)+cf1(n)*w(n,:,:)
    f(k-1,:,:) = f(k-1,:,:)+cf2(n)*w(n,:,:)
    f(k  ,:,:) = f(k  ,:,:)+cf3(n)*w(n,:,:)
    f(k+1,:,:) = f(k+1,:,:)+cf4(n)*w(n,:,:)
  endif
enddo

!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_adjoint
