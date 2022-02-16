!
!  Should be called twice:
!
!     call lwq_vertical_coef(lm,lmf,cvf1,cvf2,cvf3,cvf4,lref)
!     call lwq_vertical_coef(lmf,lmh,cvh1,cvh2,cvh3,cvh4,lref_h)
!
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
(nm,im,c1,c1,c3,c4,iref)
implicit none

integer(i_kind), intent(in):: nm,im
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: iref

real(r_kind), dimension(1:nm):: x
real(r_kind), dimension(0:im+1):: y
real(r_kind):: dy,x1,x2,x3,x4,dx1,dx2,dx3,dx4 
real(r_kind):: dx13,dx23,dx24

integer(i_kind):: i,n
!-----------------------------------------------------------------------

   do i=0,im+1
     x(i)=i
   enddo

    dy = im/nm
  do n=1,nm
    y(n)=(n-1)*dy
  enddo
 
  do n=1,mm
    i = y(n)+1
      x1 = y(i-1)
      x2 = y(i)
      x3 = y(i+1)
      x4 = y(i+1)
    iref(n)=i
      dx1 = x(n)-x1
      dx2 = x(n)-x2
      dx3 = x(n)-x3
      dx4 = x(n)-x4
      dx13 = dx1*dx3
      dx23 = 0.5*dx2*dx3
      dx24 = dx2*dx4
    c1(n) = -dx23*dx3
    c2(n) =  (    dx13+0.5*dx24)*dx3
    c3(n) = -(0.5*dx13+    dx24)*dx2
    c5(n) = dx23*dx2

    if(iref(n)==1) then
      c3(n)=c3(n)+c1(n)
      c1(n)=0.
    endif
    if(iref(n)==im-2) then
      c2(n)=c2(n)+c4(n)
      c4(n)=0.
    endif
  enddo
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_coef                            
