    Program visco
      use mpi
      use readsnap
      use fields
      implicit none
      
      ! Our programming variables
      Integer :: N_sample = 20000
      Integer :: downsample_size = 200000
      Real*4, allocatable :: Random_array(:,:), Random_sample_position(:,:)
      Real*4, allocatable :: xds(:,:), vds(:,:)
      Real*4, allocatable :: x(:,:), v(:,:)
      Real*4  :: Box_size = 60 ! in Mpc/h
      Real*8  :: a
      Real*8  :: Lambda, radius, m,pi,G,c,h

      ! Varialbles for the manager worker algorithm
      integer :: ierr
      integer :: rank,processors
      integer :: id
      integer :: rootprocess = 0
      Integer :: p_ids, p_idr ! particles per processor sent and received
      Integer :: part_id
      Integer :: prc_done, prc_fin = 0
      integer :: status(MPI_STATUS_SIZE)
      
      ! for writing to the file
      character(len=8) :: fmt, x1
      character(len=50) :: filename ! format descriptor
      Integer :: i

      COMMON /VARS/    Lambda,radius,m,pi,G,c,h

      h      = 0.67
      Lambda = (h/3.0) ! in Mpc
      radius = 7.0/(Lambda)
      m      = (128**3/downsample_size)*0.7146691620   ! in the units of 10^10 Msun
      pi     = 3.14159265358979
      G      = 43.0208   ! (Mpc/M_10sun)*(km/s^2)
      c      =  3*10**5 !9.72*1.0d-15   ! Mpc/sec

    
      allocate(Random_array(N_sample,3))
      allocate(Random_sample_position(N_sample,3))

      call MPI_INIT ( ierr )
      call MPI_comm_size(MPI_comm_World, processors, ierr)
      call MPI_comm_rank(MPI_comm_World, rank, ierr)
      
      
      if (rank == rootprocess) then
          
          ! write the position file
          call Random_number(Random_array)
          Random_sample_position = Box_size*Random_array
          open(10000,file='./data/positions.dat')
          do i=1,N_sample
              write(10000,*) Random_sample_position(i,:)
          end do
          close(10000)

          do id=1,processors-1
            p_ids = id
            print *, 'sending particle id',p_ids, 'to processor', id
            call MPI_SEND(p_ids,1, MPI_INT, id, 0, MPI_comm_World, ierr)
          end do

          part_id = processors-1
          do
            call MPI_RECV(prc_done,1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_comm_World, status, ierr)
            if (prc_done == -1) then
              prc_fin = prc_fin +1
            else
              if (part_id<N_sample) then
                part_id = part_id + 1
              else
                part_id = N_sample+1
              end if
              call MPI_SEND(part_id,1, MPI_INT, prc_done, 0, MPI_comm_World, ierr)
            end if

            if (prc_fin == processors-1) then
              print *, 'Master node is exiting'
            end if

            if (prc_fin == processors-1) Exit
          end do

      else
!           call MPI_comm_rank(MPI_comm_World, rank, ierr)
!           print *, 'processor', rank, 'received', p_idr
          call Random_number(Random_array)
          Random_sample_position = Box_size*Random_array

          call readposvel(x,v,a)
          call downsample(x,downsample_size,xds)
          call downsample(v,downsample_size,vds)
          xds = xds/1000.0

          do 
            call MPI_RECV(p_idr,1, MPI_INT, rootprocess, MPI_ANY_TAG, MPI_comm_World, status, ierr)

            !opening file for writing data
            fmt = '(I3.3)' ! an integer of width 5 with zeros at the left
            write (x1,fmt) rank ! converting integer to string using a 'internal file'
!             print *, 'opening file on slave nodes', rank,'....'
            filename = './data/tf'//trim(x1)//'.dat'
            open(unit = rank+10,file=filename)

            if (p_idr <= N_sample) then
              print *, 'Processor', rank, 'Calculating', p_idr
              call calculate_fields(Random_sample_position(p_idr,:),xds,vds,a,rank)
            else
              print *, 'processor', rank, 'exiting in this loop'
              rank = -1
              close(rank+10) 
            end if

!             print *, 'Communicating back'
            call MPI_SEND(rank,1, MPI_INT, rootprocess, 0, MPI_comm_World, ierr)
            
            if (p_idr > N_sample) Exit
          end do


      end if
      call MPI_FINALIZE ( ierr )
      stop
    end Program visco

