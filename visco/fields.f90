MODULE fields
  IMPLICIT NONE

  TYPE primaryfields
     REAL*8                    ::  rho
     REAL*8, DIMENSION(3)      ::  d_rho
     REAL*8, DIMENSION(3,3)    ::  dd_rho
     REAL*8, DIMENSION(3,3)    ::  d2d_rho
     REAL*8, DIMENSION(3)      ::  xi
     REAL*8, DIMENSION(3,3)    ::  dd_sigma
     REAL*8, DIMENSION(3)      ::  p
     REAL*8, DIMENSION(3,3)    ::  d_pi
     REAL*8, DIMENSION(3,3,3)  ::  dd_pi
     REAL*8, DIMENSION(3,3,3)  ::  d2d_pi
     REAL*8, DIMENSION(3,3)    ::  d_rho_d_phi
  END TYPE primaryfields

  TYPE secondaryfields
     REAL*8, DIMENSION(3,3)    ::  dd_phi
     REAL*8, DIMENSION(3)      ::  v
     REAL*8, DIMENSION(3,3)    ::  d_v
     REAL*8, DIMENSION(3,3,3)  ::  dd_v
     REAL*8, DIMENSION(3,3)    ::  d2d_v
     REAL*8, DIMENSION(3,3)    ::  dd_kappa
  END TYPE secondaryfields

  TYPE tertiaryfields
     REAL*8  ::  delta
     REAL*8  ::  theta
     REAL*8  ::  As
     REAL*8  ::  d2delta
     REAL*8  ::  d2theta
  END TYPE tertiaryfields

CONTAINS
  !==============================================
  ! Calculate distance between two points
  !==============================================

  SUBROUTINE euclidean_dist(r,rn,d)
    implicit none
    REAL*4, DIMENSION(3), INTENT(in) ::  r, rn
    REAL*4, INTENT(out):: d

    d = ( dot_product(r-rn,r-rn) )**0.5
    !   print *, d
  END SUBROUTINE euclidean_dist







  SUBROUTINE calculate_xi(snapdata,xi_rshp)
    use mpi
    use readsnap
    use constants
    implicit none 
    Type(snapshotdata)  :: snapdata
    REAL*4              :: d                    ! A parameter for euclidean dist
    REAL*4, allocatable :: x(:,:)
    Real*4, allocatable :: xi(:), xi_sub(:,:), xi_rshp(:,:)

    ! Varibles for calculation of the xi field
    REAL*8              :: mod_rnrm,Wrnrm,xi_kernel
    Integer             :: prtcl_id, eta_m_ind,i,j
    Real                :: Box_size = 60
    ! variables for parallelization
    integer             :: ierr
    integer             :: rank,processors,start_index,stop_index, size_x, ppp
    

    
    size_x = Size(snapdata%positions,1)
    !Tell each processor how many particles to calculate
    call MPI_comm_rank(MPI_comm_World, rank, ierr)
    call MPI_comm_size(MPI_comm_World, processors, ierr)
    ppp = size_x/(processors)


    ! Each processor at this point has received the positions data.
    ! Reshape the data to its original form
    allocate(x(size_x,3))  
    x = snapdata%positions

    ! print *, 'x=',x(1,:)

    start_index = int((rank*ppp)+1)
    stop_index  = int((rank+1)*ppp)
    ! print *, 'rank start/stop index:', rank, start_index, stop_index



    ! print *, 'x', x
    allocate(xi_sub(ppp,3)) 
     
    DO prtcl_id = start_index, stop_index
      xi_sub = 0.0


      ! print *,'xisub1=',xi_sub(:,:)
      do eta_m_ind=1, SIZE(x,1)
        
        !Periodic bondary conditions 
         if (x(prtcl_id,1) < Lambda) then
            if (x(eta_m_ind,1) > Box_size - Lambda) then
              x(eta_m_ind,1) = x(eta_m_ind,1) - Box_size
            end if
         end if

         if (x(prtcl_id,1) > Box_size - Lambda) then
            if (x(eta_m_ind,1) < Lambda) then
              x(eta_m_ind,1) = x(eta_m_ind,1) + Box_size
            end if
         end if

         if (x(prtcl_id,2) < Lambda) then
            if (x(eta_m_ind,2) > Box_size - Lambda) then
              x(eta_m_ind,2) = x(eta_m_ind,2) - Box_size
            end if
         end if

         if (x(prtcl_id,2) > Box_size - Lambda) then
            if (x(eta_m_ind,2) < Lambda) then
              x(eta_m_ind,2) = x(eta_m_ind,2) + Box_size
            end if
         end if

         if (x(prtcl_id,3) < Lambda) then
            if (x(eta_m_ind,3) > Box_size - Lambda) then
              x(eta_m_ind,3) = x(eta_m_ind,3) - Box_size
            end if
         end if

         if (x(prtcl_id,3) > Box_size - Lambda) then
            if (x(eta_m_ind,3) < Lambda) then
              x(eta_m_ind,3) = x(eta_m_ind,3) + Box_size
            end if
         end if
        

         Call euclidean_dist(x(prtcl_id,:),x(eta_m_ind,:),d)
         If(d .le. radius .and. d > 0.1) then
            mod_rnrm = (dot_product( x(prtcl_id,:)-x(eta_m_ind,:) , x(prtcl_id,:)-x(eta_m_ind,:) ))**0.5
            Wrnrm = EXP( -0.5*Lambda**2*dot_product( x(prtcl_id,:)-x(eta_m_ind,:) , x(prtcl_id,:)-x(eta_m_ind,:) ) )
            xi_kernel =( (1.0/mod_rnrm) * ERFC(Lambda*mod_rnrm/sqrt(2.0)) & 
              + (2.0/pi)**0.5*Lambda*Wrnrm )*((1/mod_rnrm**2)*exp(-6*(mod_rnrm/radius)**6))!*GAMMA(mod_rnrm*Lambda)
            do i=1,3
                ! xi_sub(prtcl_id-rank*ppp,i) = xi_sub(prtcl_id-rank*ppp,i) + (x(prtcl_id-rank*ppp,i)-x(eta_m_ind,i))*xi_kernel
              xi_sub(prtcl_id-rank*ppp,i) = i 
            end do
         end if

      end do
      ! print *,'xisub=',xi_sub(prtcl_id,:)

    END DO
    ! do j =1,ppp
    !   print *, 'xi',xi_sub(j,:)
    ! end do
    ! print *, ''
    ! print *,'xisub=',xi_sub
    call MPI_Barrier(MPI_comm_World,ierr)

    ! At this point each processor has calculated its chunk of xi field. Combine these chunks
    allocate(xi(3*size_x))
    call MPI_Allgather(xi_sub, ppp*3, Mpi_Real, xi, ppp*3, Mpi_Real, MPI_COMM_WORLD, ierr)
    ! print *, 'Calculating the xi field'


    ! Reshape this combined 1-D array to get back the calculated xi field
    allocate(xi_rshp(size_x,3))
    xi_rshp = Reshape(xi,(/size_x,3/))
    ! print *, SIZE(xi,1)
    ! call MPI_FINALIZE ( ierr )

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! print *, "xi field final", xi_rshp
    


21 END SUBROUTINE calculate_xi









  !==========================================================================================
  !                Calculate primary fields
  !               ++++++++++++++++++++++++++
  ! This subroutine takes the value of the processor number (pro_id), a position in the cube (r),
  ! the positions of all particles read from snapshot (x), the velocities (v), the model parameters
  ! (Lambda) and (radius) and total number of sampled points (N_sample) and then calculates the
  ! W fields at the given position r
  !==========================================================================================

  SUBROUTINE primary_fields(r,snapdata,pf)
    use constants
    use readsnap

    IMPLICIT NONE
    ! Integer, intent(in) :: N_sample ! The number of sampled points in the simulation cube
    ! Real*8, intent(in)  :: a
    REAL*4, allocatable   :: x(:,:)
    REAL*4, allocatable   :: v(:,:)
    REAL*4, INTENT(in)  :: r(3)

    REAL*4  :: d                    ! A parameter for euclidean dist

    TYPE(primaryfields) :: pf
    Type(snapshotdata)  :: snapdata

    ! We now define the quantities required for calculation of W
    REAL*8  :: Wkernel
    REAL*8  :: d_Wkernel(3)
    REAL*8  :: dd_Wkernel(3,3)
    REAL*8  :: d2d_Wkernel(3,3)
    REAL*8  :: mod_rnrm,Wrnrm,xi_kernel
    INTEGER :: i,j,k
    INTEGER :: eta_n_ind, eta_m_ind
    REAL*8  :: norm1,norm
    Real*8  :: dp

    ! COMMON /VARS/    m

    pf%rho         = 0.0
    pf%d_rho       = 0.0
    pf%dd_rho      = 0.0
    pf%d2d_rho     = 0.0
    pf%p           = 0.0
    pf%d_pi        = 0.0
    pf%dd_pi       = 0.0
    pf%d2d_pi      = 0.0
    pf%dd_sigma    = 0.0
    pf%d_rho_d_phi = 0.0
    pf%xi          = 0.0


      ! print *, 'size = ' , size(x,1)
    x = snapdata%positions
    v = snapdata%velocities

    ! print *, "calculating xi field"
    DO eta_n_ind=1,SIZE(x,1)
       ! print *, 'r and x', r,x(eta_n_ind,:)
       CALL euclidean_dist(r,x(eta_n_ind,:),d)
             ! print *, d
       IF (d<radius) THEN

          !Sum the W first
          Wkernel = EXP( -0.5*Lambda**2*dot_product(r-x(eta_n_ind,:),r-x(eta_n_ind,:)) )
          pf%rho = pf%rho + Wkernel
          dp = dot_product(r-x(eta_n_ind,:),r-x(eta_n_ind,:))
          ! PRINT *, "d,radius", d,radius

          ! print *, pf%rho

          !Now sum the d_W
          DO i=1,3
             d_Wkernel(i) = -Lambda**2*r(i)*Wkernel
             pf%d_rho(i) = pf%d_rho(i) + d_Wkernel(i)
          END DO

          !Now sum the second derivative of W
          DO i=1,3
             DO j=1,3
                IF (i==j) THEN
                   dd_Wkernel(i,i) = Lambda**2*Wkernel*(Lambda**2*(r(i)-x(eta_n_ind,i))**2 - 1)
                ELSE
                   dd_Wkernel(i,j) = Lambda**2*Wkernel*(Lambda**2*(r(i)-x(eta_n_ind,i))*(r(j)-x(eta_n_ind,j)))
                END IF
                pf%dd_rho(i,j) = pf%dd_rho(i,j) + dd_Wkernel(i,j)
             END DO
          END DO

          !Now sum the third derivatives
          DO i=1,3
             DO j=1,3
                IF (i==j) THEN
                   d2d_Wkernel(i,i) = Lambda**2*( d_Wkernel(i)*(Lambda**2*(r(i)-x(eta_n_ind,i))**2 - 1) &
                        + Wkernel*Lambda**2*2*(r(i)-x(eta_n_ind,i)) )
                ELSE
                   d2d_Wkernel(i,j) = Lambda**2*( d_Wkernel(i)*(Lambda**2*(r(i)-x(eta_n_ind,i))*(r(j)-x(eta_n_ind,j))) &
                        + Wkernel*Lambda**2*(r(j)-x(eta_n_ind,j)) )
                END IF
                pf%d2d_rho(i,j) = pf%d2d_rho(i,j) + d2d_Wkernel(i,j)
             END DO
          END DO
          ! Calculate the xi field now by constructing another spehere of around x(i,:)
          DO eta_m_ind=1, SIZE(x,1)
             CALL euclidean_dist(x(eta_n_ind,:),x(eta_m_ind,:),d)
!             print *, d, radius

             IF(d .le. radius .and. d > 10.0) THEN
             	! print*, "calculating xi"
                ! calculate xi
                mod_rnrm = (dot_product( x(eta_n_ind,:)-x(eta_m_ind,:) , x(eta_n_ind,:)-x(eta_m_ind,:) ))**0.5
               ! print *, "mod_rnrm =", mod_rnrm, eta_m_ind, eta_n_ind
                Wrnrm = EXP( -0.5*Lambda**2*dot_product( x(eta_n_ind,:)-x(eta_m_ind,:) , x(eta_n_ind,:)-x(eta_m_ind,:) ) )
!                print *, "Wrnrm = ", Wrnrm
                xi_kernel =( (1.0/mod_rnrm) * ERFC(Lambda*mod_rnrm/sqrt(2.0)) & 
                	+ (2.0/pi)**0.5*Lambda*Wrnrm )*((1/mod_rnrm**2)*exp(-6*(mod_rnrm/radius)**6))!*GAMMA(mod_rnrm*Lambda)
               ! print *, "mod_rnrm = ", mod_rnrm,"xi_kernel = ", xi_kernel
                DO i=1,3
                   pf%xi(i) = pf%xi(i) + (x(eta_n_ind,1)-x(eta_m_ind,1))*xi_kernel
                  ! print *, pf%xi(i)
                END DO
             END IF
          END DO

          ! Calculate the pi fields now
          DO i = 1,3
             pf%p(i) = pf%p(i) + pf%p(i)*pf%rho
          END DO
          ! The first derivative of pi field
          DO i = 1,3
             DO j = 1,3
                pf%d_pi(i,j) = pf%d_pi(i,j) +  d_Wkernel(i)*v(eta_n_ind,j)
             END DO
          END DO
          ! The second derivative of the pi field
          DO i=1,3
             DO j=1,3
                DO k=1,3
                   pf%dd_pi(i,j,k) = pf%dd_pi(i,j,k) + dd_Wkernel(i,j)*v(eta_n_ind,k)
                END DO
             END DO
          END DO
          ! The third derivative
          DO i=1,3
             DO j=1,3
                DO k=1,3
                   pf%d2d_pi(i,j,k) = pf%d2d_pi(i,j,k) + d2d_Wkernel(i,j)*v(eta_n_ind,k)
                END DO
             END DO
          END DO
          ! The sigma terms
          DO i=1,3
             DO j=1,3
                pf%dd_sigma(i,j) = pf%dd_sigma(i,j) + dd_Wkernel(i,j)*v(eta_n_ind,i)*v(eta_n_ind,j)
             END DO
          END DO

          ! print *, "here1"
          DO i=1,3
             DO j=1,3
                pf%d_rho_d_phi(i,j) = pf%d_rho_d_phi(i,j) + d_Wkernel(i)*pf%xi(j)
               ! print *,  d_Wkernel(i), pf%xi(j)
             END DO
             ! print *, pf%d_rho_d_phi
          END DO
          ! print *, eta_n_ind
       ! ELSE
       ! print *, "pf%rho =", pf%rho
       END IF
    END DO

    ! Normalize the fields now
    norm  = (snapdata%mass*Lambda**3)/(snapdata%a**3*(2*pi)**(3/2))
    pf%rho     = norm*pf%rho
    ! print *,"Sum W", pf%rho, norm
    ! print *
    pf%d_rho   = norm*pf%d_rho
    pf%dd_rho  = norm*pf%dd_rho
    pf%d2d_rho = norm*pf%d2d_rho

    norm1 = (snapdata%mass**2*G)/(snapdata%a**4*c**2) *(Lambda/(sqrt(2*pi)))**3
    ! norm1 = (m**2*G)/(a**4) *(Lambda/(sqrt(2*pi)))**5

    pf%d_rho_d_phi = norm1*pf%d_rho_d_phi

    !velocity field normalization
    pf%p        = (1/c)*norm*pf%p
    pf%d_pi     = (1/c)*norm*pf%d_pi
    pf%dd_pi    = (1/c)*norm*pf%dd_pi
    pf%d2d_pi   = (1/c)*norm*pf%d2d_pi
    pf%dd_sigma = (1/c**2)*norm*pf%dd_sigma

    !    print *, a, Lambda, radius, m, norm, x(1,:)
    ! Print *, 'normalized rho = ', norm,m,a,Lambda
    !   Print *, 'normalized d_rho = ', norm*pf%d_rho    

  END SUBROUTINE primary_fields

  !****************************************************************************************************************

  SUBROUTINE secondary_fields(pf,a,sf)
    use constants
    implicit none
    TYPE(primaryfields), INTENT(in) :: pf
    TYPE(secondaryfields) :: sf
    INTEGER :: i,j,k
    REAL*8, INTENT(in) :: a
    real*8  :: norm2
    ! COMMON /VARS/    m

    norm2 = (4.0*pi*G*a**2)/(2*Lambda**2*c**2)
    ! norm2 = (4.0*pi*G*a**2)/(2*Lambda**2)

    !    print *, a, Lambda,norm2
    ! print *, "calculating secondary field"
    DO i =1,3
       DO j=1,3
          sf%dd_phi(i,j) = pf%d_rho_d_phi(i,j) + norm2*pf%d_rho(i)*pf%d_rho(j) + pf%rho*pf%dd_rho(i,j)
          !          print *, pf%d_rho_d_phi(i,j), sf%dd_phi(i,j)
       END DO
    END DO

    DO i=1,3
       sf%v(i) = pf%p(i)/(pf%rho)
    END DO

    DO i=1,3
       DO j=1,3
          sf%d_v(i,j) = pf%d_pi(i,j)/pf%rho - pf%p(j)*pf%d_rho(i)/pf%rho**2
          ! print *,"rho, d_pi",pf%rho,pf%d_pi(i,j)
       END DO
    END DO

    DO i=1,3
       DO j=1,3
          DO k=1,3
             sf%dd_v(i,j,k) = pf%dd_pi(i,j,k)/pf%rho - pf%d_pi(i,k)*pf%d_rho(j)/pf%rho**2 - pf%d_pi(j,k)*pf%d_rho(i)/pf%rho**2 &
                  - pf%p(k)*pf%dd_rho(i,j)/pf%rho**2 + 2*pf%p(k)*pf%d_rho(i)*pf%d_rho(j)/pf%rho**3
          END DO
       END DO
    END DO

    DO i=1,3
       DO j=1,3
          sf%d2d_v(i,j) = pf%rho**(-4)*(-6*pf%d_rho(i)**2*pf%d_rho(j)*pf%p(j) &
               + 2*pf%rho*( pf%d_rho(i)**2*pf%d_pi(j,j) + pf%dd_rho(i,i)*pf%d_rho(j)*pf%p(j) &
               + pf%d_rho(j)*(pf%d_pi(i,j)*pf%d_rho(j) + pf%dd_rho(i,j)*pf%p(j))  ) &
               - pf%rho**2*(2*pf%dd_rho(i,j)*pf%d_pi(i,j) & 
               +2*pf%dd_pi(i,j,j)*pf%d_rho(i) + pf%dd_rho(i,i)*pf%d_pi(j,j) + pf%dd_pi(i,i,j)*pf%d_rho(j) & 
               + pf%d2d_rho(i,j)*pf%p(j)) & 
               + pf%rho**3*pf%d2d_pi(i,i,j) )
       END DO
    END DO

    DO i=1,3
       DO j=1,3
          sf%dd_kappa(i,j) = pf%dd_sigma(i,j) - pf%d_pi(i,j)*sf%d_v(i,j) - pf%d_pi(i,i)*sf%d_v(j,j) - pf%p(i)*sf%dd_v(i,j,j) &
               - pf%dd_pi(i,j,i)*sf%v(j)
       END DO
    END DO

  END SUBROUTINE secondary_fields

  ! !********************************************************************************************************************

  SUBROUTINE tertiary_fields(a,pf,sf,tf)
    use constants
    implicit none
    TYPE(primaryfields),   INTENT(in) :: pf
    TYPE(secondaryfields), INTENT(in) :: sf

    TYPE(tertiaryfields)  :: tf
    REAL*8, INTENT(in) :: a
    INTEGER :: i,j
    REAL*8 :: sum
    REAL*8 :: rho_bg
    Real*8 Hubble,Hubble0,omega_m0
    ! Real*8  :: m
    ! COMMON /VARS/    m

    omega_m0 = 2.67d-1
    Hubble0  = 100*h !(3.24*h*1.0d-18)
    Hubble = Hubble0*sqrt(omega_m0/(a**3) + 1.0-omega_m0)
    ! rho_bg = (omega_m0/a**3)*((3.d0*Hubble0**2)/(8.d0*pi*G))
    rho_bg = 5.82/a**3
    print *, "rho_bg  =", rho_bg
    tf%delta   = pf%rho/rho_bg - 1
    print *, "delta   =", tf%delta

    sum = 0
    DO i=1,3
       sum = sum + sf%d_v(i,i)
       ! print *, sum
    END DO
    tf%theta   = (-1/(Hubble*a))*sum
    print *, "theta   =", tf%theta

    sum=0
    DO i = 1,3
       DO j=1,3
          sum = sum + sf%dd_kappa(i,j)
          !            print *, "sum j=", sum
       END DO
       sum = sum + sf%dd_phi(i,i)
       !            print *, "sum i=", sf%dd_phi(i,i)      
    END DO
    tf%As      = (1.0/rho_bg)*sum
    print *, "As      =", tf%As
    sum=0
    DO i=1,3
       sum = sum + pf%dd_rho(i,i)
    END DO
    tf%d2delta = (1/rho_bg)*sum
    print *, "d2delta =", tf%d2delta
    sum=0
    DO i=1,3
       DO j=1,3
          sum = sum + sf%d2d_v(i,j)
       END DO
    END DO
    tf%d2theta = (-1/(Hubble*a)) * sum
    print *, "d2theta =", tf%d2theta
    print *, " "

  END SUBROUTINE tertiary_fields

  subroutine calculate_fields(r,snapdata, rank)
    use readsnap
    implicit none
    Integer :: rank
    Real    :: r(3)
    type(primaryfields)   :: pf
    type(secondaryfields) :: sf
    type(tertiaryfields)  :: tf
    type(snapshotdata)    :: snapdata
    ! Real    :: x(:,:),v(:,:)
    ! Real*8  :: a


    call primary_fields(r,snapdata,pf)
    call secondary_fields(pf,snapdata%a,sf)
    call tertiary_fields(snapdata%a,pf,sf,tf)
    call write_to_file(pf,tf, rank)

  end subroutine calculate_fields

  ! !********************************************************************************************************************

    subroutine write_to_file(pf,tf,rank)
      implicit none
      type(tertiaryfields) :: tf
      type(primaryfields)  :: pf
      ! character(len=8) :: fmt, x1
      ! character(len=50) :: filename ! format descriptor
      integer :: rank
      ! integer :: i1

      ! fmt = '(I3.3)' ! an integer of width 5 with zeros at the left
      ! i1  = rank
      ! write (x1,fmt) i1 ! converting integer to string using a 'internal file'
      ! filename = 'data/tertiary_fields'//trim(x1)//'.dat'
      ! open(unit = rank,file=filename)
      write(rank+10,"(/ F15.6, F15.6, E15.6, E15.6, E15.6, E15.6)", advance='no') pf%rho, tf%delta, &
      tf%theta, tf%As, tf%d2delta, tf%d2theta
      ! write(rank+10,*) rank, pf%rho

     ! close(rank) 
    end subroutine write_to_file

  !*********************************************************************************************************************


END MODULE fields
