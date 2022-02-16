!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoy_2d_gh                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in u direction, including values at the  !
! edges of the subdomains and assuming mirror boundary conditions at   !
! end of domain. Version for high generations                          !
!                                                                      !
!**********************************************************************!
(Warray,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout)::      &
                                  Warray
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S              &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
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
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


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


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby

!
!  SEND boundaries to SOUTH and NORTH
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
                              mpi_comm_work, sHandle(3), isend)
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
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,0:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
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

! From south

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

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoy_2d_gh
