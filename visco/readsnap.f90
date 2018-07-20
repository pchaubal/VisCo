module readsnap 
implicit none 

   Type snapshotdata
      character*200       ::  filename
      real*4,allocatable  ::  positions(:,:)
      real*4,allocatable  ::  velocities(:,:)
      real*8              ::  a
      real*8              ::  mass
      integer             ::  n_particles
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

   print *,'a,m:', a,snapdata%mass


   allocate(snapdata%positions(1:N,1:3))
   snapdata%positions = Transpose(pos)/1000.0 ! kpc to Mpc 

   allocate(snapdata%velocities(1:N,1:3))
   snapdata%velocities = Transpose(vel)
   
   deallocate(pos)
   deallocate(vel)
   
   print *,'Done with reading snapshot'
 
   ! At this point, the coordinates of all particles of type 1 will
   ! be stored in the snapdata  
end subroutine readposvel



subroutine read_illustris(snapdata)
   character*200        :: filepos,filevel
   real*8, allocatable  :: pos(:,:)
   real*8, allocatable  :: vel(:,:)
   type(snapshotdata)   :: snapdata

   filepos= '/home/prakrut/codes/Illustris-3/binarysnaps/positions.bin'
   filevel= '/home/prakrut/codes/Illustris-3/binarysnaps/velocities.bin'

   snapdata%mass = 1.0
   snapdata%a    = 1.0

   open (1, file=filepos, form='unformatted', access='stream')

   allocate(pos(1,1))
   pos = 0.0
   read (1) pos
   snapdata%n_particles = int(pos(1,1) +1)
   deallocate(pos)

   allocate(pos(snapdata%n_particles,3))
   pos = 0.0
   rewind(1)
   read (1) pos

   allocate(snapdata%positions(snapdata%n_particles-1,3))
   snapdata%positions(1:snapdata%n_particles-1,:) = pos(2:snapdata%n_particles,:)/1000.0 !Mpc
   deallocate(pos)

   close (1)


   open (12, file=filevel, form='unformatted', access='stream')
   allocate(vel(snapdata%n_particles-1,3))
   read (12) vel

   allocate(snapdata%velocities(snapdata%n_particles-1,3))
   snapdata%velocities = vel/sqrt(snapdata%a)

   deallocate(vel)

   close (12)


end subroutine read_illustris

end module readsnap

