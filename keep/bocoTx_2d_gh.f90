!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTx_2d_gh                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nby lines close to West and East edges of the subdomains,   !
! including values at the edges, from halos of North and South         !
! neighbors assuming mirror boundary conditions at the boundaries of   !
! the whole domain. Version for high generations                       ! 
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   
integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+1)*(nbx+1)
!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=0,jmax
              do i=-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=0,jmax
              do i=0,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=0,jmax
     do i=0,nbx
       W(:,i,j)= W(:,i,j)+W(:,-i,j)
     end do
     end do
   else
     do j=0,jmax
     do i=0,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=0,jmax
     do i=0,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx-i,j)
     end do
     end do
   else 
     do j=0,jmax
     do i=0,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoTx_2d_gh
