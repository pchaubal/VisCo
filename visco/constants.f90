    Module constants
    implicit none
	Integer,parameter  :: N_sample = 560000
    Integer,parameter  :: downsample_size = 40000
    Real*8, parameter  :: Box_size = 60.0 ! in Mpc/h


    Real*8, parameter  :: h      = 0.67                                                                
    Real*8, parameter  :: Lambda = (h/3.0) ! in Mpc                                                    
    Real*8, parameter  :: radius = 7.0/(Lambda)                                                        
    Real*8, parameter  :: pi     = 3.14159265358979                                                    
    Real*8, parameter  :: G      = 43.0208   ! (Mpc/M_10sun)*(km/s)^2                                  
    Real*8, parameter  :: c      = 3*10**5 !9.72*1.0d-15   ! Mpc/sec 



	end module constants