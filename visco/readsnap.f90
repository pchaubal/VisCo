module readsnap 
implicit none 

   Type snapshotdata
      character*200       ::  filename
      real*4,allocatable  ::  positions(:,:)
      real*4,allocatable  ::  velocities(:,:)
      real*8              ::  a
      real*8              ::  mass
      integer*4           ::  n_particles
   end type snapshotdata


contains      

subroutine readposvel(snapdata)
   integer*4          :: flag_sfr,flag_feedback
   integer*4          :: npart(0:5), nall(0:5)
   real*8             :: massarr(0:5)
   real*8             :: a
   real*8             :: redshift
   integer*4          :: unused(34)
   integer            :: N
   real*4,allocatable :: pos(:,:)
   real*4,allocatable :: vel(:,:)

   type(snapshotdata) :: snapdata

   
   open (1, file=snapdata%filename, form='unformatted')
   read (1) npart, massarr ,a, redshift, flag_sfr,flag_feedback, nall, unused


   N=sum(npart)
   allocate(pos(1:3,1:N))
   allocate(vel(1:3,1:N))

   read (1) pos
   read (1) vel
   close (1)
   
   
   print *,'reading snapshot from:  ', snapdata%filename

   snapdata%a           = a
   snapdata%mass        = massarr(1)
   snapdata%n_particles = int(N)


   allocate(snapdata%positions(1:N,1:3))
   snapdata%positions = Transpose(pos)

   allocate(snapdata%velocities(1:N,1:3))
   snapdata%velocities = Transpose(vel)
   
   deallocate(pos)
   deallocate(vel)
   
   print *,'Done with reading snapshot'
 
   ! At this point, the coordinates of all particles of type 1 will
   ! be stored in the snapdata  
end subroutine readposvel


end module readsnap

