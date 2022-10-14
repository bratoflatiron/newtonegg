!*********************************************************************************************
!A module containing all the constant values
!*********************************************************************************************
module constants
    implicit none
    real(kind=8), parameter :: pi   = 4.d0*atan(1.d0); ! Define great pi
    real(kind=8), parameter :: sigma= 7.0d0  ! Active forcing: sigma
    real(kind=8), parameter :: c0   = 12.d0  ! Surface concentration of MT
    real(kind=8), parameter :: beta = 0.8d0  ! Geometric parameter with velocity gradient and constrain force
    real(kind=8), parameter :: drot = 0.0d0  ! Rotational diffusion 
    real(kind=8), parameter :: eps  = 1.d-7  ! Tolerance for BIEFMM


    !! Parameters for the generating curve (CONFORMAL MAP)
    real(kind=8), parameter :: rad  = 1.0d0;  ! Reference circle
    real(kind=8), parameter :: aY   = 0.15d0; ! This is equal to Y for conformal map 
    real(kind=8), parameter :: aZ   = 0.10d0; ! This is equal to Z for conformal map
    real(kind=8), parameter :: Rfa  = 5.d0;   ! Scaling of the system
    real(kind=8), parameter :: Rs   = 1.0d0/Rfa; ! This is a scaling factor for r
    real(kind=8), parameter :: Zs   = 1.0d0/Rfa; ! This is a scaling factor for z
    
    real(kind=8), parameter :: delta= rad*0.85d0/Rfa;  ! Distance from the slip surface


end module constants
    