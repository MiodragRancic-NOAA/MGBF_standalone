!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTy_2d_gh                         &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies nby lines close to North and South edges of the subdomains, !
! including values at the edges, from halos of North and South         !
! neighbors assuming mirror boundary conditions at the boundaries of   !
! the whole domain. Version for high generations                       !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fleast,Flsouth,Flnorth                             &
                    ,Fitarg_n,Fitarg_s
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,-nbx:im+nbx,-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::  sBuf_N,sBuf_S            &
                                              ,rBuf_N,rBuf_S
integer(i_kind) itarg_n,itarg_s,itarg_e,imax,jmax
logical least,lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
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
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

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


!----------------------------------------------------------------------
!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,0:imax,0:nby), stat = iaerr )

              do j=-nby,0
              do i=0,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,0:imax,0:nby), stat = iaerr )

              do j=0,nby
              do i=0,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,0:imax,0:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+W(:,i,-j)
     end do
     end do
   else
     do j=0,nby
     do i=0,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby-j)
     enddo
     enddo
   else
     do j=0,nby
     do i=0,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoTy_2d_gh
