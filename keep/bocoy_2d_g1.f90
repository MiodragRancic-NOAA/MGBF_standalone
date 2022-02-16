!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoy_2d_g1                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in y direction, including values at the  !
! edges of the subdomains and assuming mirror boundary conditions at   !
! end of domain. Version for generation 1                              !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Flsouth,Flnorth                     

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout)::      &
                                  Warray
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby


!
!  SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_S(:,i,j) = Warray(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,0:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=0,imax
                    sBuf_N(:,i,j)=Warray(:,i,jmax-nby-1+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if

!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=Warray(:,i,jmax-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From SOUTH

   if(lsouth) then

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=Warray(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=0,imax
       Warray(:,i,-nby-1+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif



!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoy_2d_g1
