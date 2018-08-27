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

subroutine downsample_snap(snapdata,downsample_size)
   implicit none
   Integer, intent(in) :: downsample_size
   type(snapshotdata)  :: snapdata
   Real*4, allocatable :: xds(:,:), vds(:,:)


   call downsample(snapdata%positions,downsample_size,xds)
   call downsample(snapdata%velocities,downsample_size,vds)
!       xds = xds/1000.0 ! this is only because the snapdata is in kpc and output data in Mpc
   deallocate(snapdata%positions)
   deallocate(snapdata%velocities)
   allocate(snapdata%positions(downsample_size,3))
   allocate(snapdata%velocities(downsample_size,3))

   snapdata%positions  = xds
   snapdata%velocities = vds
   snapdata%mass       = snapdata%mass*(snapdata%n_particles/downsample_size)
   snapdata%n_particles = downsample_size
end subroutine downsample_snap


subroutine downsample(array,downsample_size,downsampled)
   ! USE IFPORT
   implicit none
   real :: array(:,:)
   real, allocatable :: downsampled(:,:)
   integer, intent(in) :: downsample_size
   integer :: i,random_index
   integer :: array_size

   allocate(downsampled(downsample_size,3))

   array_size = size(array,1)
   
   do i=1,downsample_size
     random_index =  int(rand(0)*(array_size+1-1))+1 
     downsampled(i,:) = array(random_index,:)
   end do

end subroutine downsample



subroutine readposvel(snapdata)
   use constants

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
   call downsample_snap(snapdata,downsample_size)

end subroutine readposvel



subroutine read_illustris(snapdata)
   use constants
!    use utils
   character*200        :: filepos,filevel
   real*8, allocatable  :: pos(:,:)
   real*8, allocatable  :: vel(:,:)
   type(snapshotdata)   :: snapdata

   filepos= '/home/prakrut/codes/Illustris-3/binarysnaps/positions.bin'
   filevel= '/home/prakrut/codes/Illustris-3/binarysnaps/velocities.bin'

   snapdata%mass = (455.0**3.0)/downsample_size
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

   call downsample_snap(snapdata,downsample_size)

!    This prints the position of particles with indices from 10 to 15 
   print *, 'The positions of snap from binary:', snapdata%positions(10:15,:)



end subroutine read_illustris

end module readsnap

