!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTx_2d_g1                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nbx lines of halos in x directions, including               !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions. Just for generation 1                                    !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=1
!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im    
          jmax = jm


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
                              mpi_comm_world, sHandle(1), isend)

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
                              mpi_comm_world, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


         allocate( rBuf_W(1:km,0:nbx,0:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_world, rHandle(1), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if

!
! Assign received halos from WEST and EAST to interrior of domains
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
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if

!-----------------------------------------------------------------------
                        endsubroutine bocoTx_2d_g1
