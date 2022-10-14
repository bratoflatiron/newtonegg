subroutine loadconstant(R, del, kchu, iref, nmid, nch2d, npts2d, nquad, npols, finout)
    use constants 
    real(kind=8), intent(out) :: R, del
    integer, intent(out)      :: kchu, iref, nmid, nch2d, npts2d, nquad, npols, finout

    R = rad; del = delta;

    kchu  = 16;      ! Order of discretization for chunking the generating curve 
    iref  = 0;       ! Number of dyadic refinements at either end points 
    nmid  = 2;       ! Number of panels in the middle section
    nch2d = 2*(iref+1) + nmid; ! number of chunks on the generating curve
    npts2d= nch2d*kchu;
    nquad = 4       ! Order of element for quadrature 
    npols = (nquad+1)*(nquad+2)/2
    finout= 0;     ! exterior vs interior problem 
end subroutine loadconstant
!*****************************************************
! Generating curve for two surfaces
!*****************************************************
subroutine curveB(kchu, iref, nmid, nch2d, npts2d, & 
                & iptype_B2D,ixys_B2D,nptchB,surcoeff_B2D,surf_B2D)
    use constants 
    implicit none 
    integer, intent(in)      :: kchu, iref, nmid, nch2d, npts2d 
    integer, intent(out)     :: iptype_B2D(nch2d), ixys_B2D(nch2d+1), nptchB
    real(kind=8), intent(out):: surcoeff_B2D(6,npts2d), surf_B2D(8,npts2d)
    real(kind=8) :: R_B, rmax

    ! Set rmax and R = rad for the bounding surface
    R_B = rad; rmax = 1.0d0*R_B*Rfa;
    call get_oocyte3d_riemann_chunks(aY,aZ,R_B,kchu,iref,nmid,nch2d,npts2d, & 
        & iptype_B2D,ixys_B2D,surcoeff_B2D,surf_B2D)

    call get_oocyte3d_riemann_fun_mem(aY,aZ,R_B,kchu,iref,nmid,nch2d,rmax,nptchB)

end subroutine curveB
!------------------------------
subroutine curveS(kchu, iref, nmid, nch2d, npts2d, & 
            & iptype_S2D,ixys_S2D,nptchSL,surcoeff_S2D,surf_S2D)
    use constants 
    implicit none 
    integer, intent(in)      :: kchu, iref, nmid, nch2d, npts2d 
    integer, intent(out)     :: iptype_S2D(nch2d), ixys_S2D(nch2d+1),nptchSL
    real(kind=8), intent(out):: surcoeff_S2D(6,npts2d), surf_S2D(8,npts2d)
    real(kind=8) :: R_B, R_SL, rmSL

    ! Estimate rmax and R = rad-delta for the slip surface
    R_B = rad; R_SL = rad-delta; rmSL = rmax*R_SL/R_B;
    
    call get_oocyte3d_riemann_chunks(aY,aZ,R_SL,kchu,iref,nmid,nch2d,npts2d, & 
        & iptype_S2D,ixys_S2D,surcoeff_S2D,surf_S2D)

    call rmxestimate(R_SL, nptchB, kchu, iref, nmid, nch2d, rmSL, nptchSL)

    call get_oocyte3d_riemann_fun_mem(aY,aZ,R_SL,kchu,iref,nmid,nch2d,rmSL,nptchSL)

end subroutine curveS
!*****************************************************
! Setup surface discretization
!*****************************************************
subroutine boundB(kchu,iref,nmid,nch2d,nptB,nptchB, &
                 & nordB, ixyzB, iptB,   surfB,   srBcoeff, wtsB)
    use constants 
    implicit none 
    integer, intent(in)      :: kchu, iref, nmid, nch2d, nptB, nptchB 
    integer, intent(out)     :: nordB(nptchB), ixyzB(nptchB+1),  iptB(nptchB)
    real(kind=8), intent(out):: surfB(12,nptB),srBcoeff(9,nptB), wtsB(nptB)
    real(kind=8) :: R_B, rmax

    !! Set rmax and R = rad for the bounding surface
    R_B = rad; rmax = 1.0d0*R_B*Rfa;

    !! Boundary surface
    call get_oocyte3d_riemann_fun_geom(aY,aZ,R_B,kchu,iref,nmid,nch2d,rmax, &
        & nquad,nptchB,nptB,nordB,ixyzB,iptB,surfB,srBcoeff)

    !! Obtain quadrature weights 
    call get_qwts(nptchB,  nordB,  ixyzB,  iptB,  nptB,    surfB,    wtsB)

end subroutine boundB
!------------------------------
subroutine boundSL(kchu,iref,nmid,nch2d,nptSL,nptchSL, &
    & nordSL, ixyzSL, iptSL, surfslip,srslipcoeff,wtsslip)
    use constants 
    implicit none 
    integer, intent(in)      :: kchu, iref, nmid, nch2d, nptSL, nptchSL 
    integer, intent(out)     :: nordSL(nptchSL),ixyzSL(nptchSL+1),iptSL(nptchSL)
    real(kind=8), intent(out):: surfslip(12,nptSL),srslipcoeff(9,nptSL),wtsslip(nptSL)
    real(kind=8) :: R_B, R_SL, rmSL

    !! Estimate rmax and R = rad-delta for the slip surface
    R_B = rad; R_SL = rad-delta; rmSL = rmax*R_SL/R_B;
    call rmxestimate(R_SL, nptchB, kchu, iref, nmid, nch2d, rmSL, nptchSL)

    !! Slip surface
    call get_oocyte3d_riemann_fun_geom(aY,aZ,R_SL,kchu,iref,nmid,nch2d,rmSL, &
                    & nquad,nptchSL,nptSL,nordSL,ixyzSL,iptSL,surfslip,srslipcoeff)

    !! Obtain quadrature weights 
    call get_qwts(nptchSL,nordSL, ixyzSL, iptSL, nptSL, surfslip, wtsslip)

end subroutine boundSL
!*****************************************************
! Interpolation matrix obtain 
!*****************************************************
subroutine get_interp_mat(kchu,iref,nmid,nch2d,nptB,nptSL,npols,nptchB,nordB,ixyzB,iptB, & 
            & xinterp, matinterp)
    use constants 
    implicit none 
    integer, intent(in)      :: kchu, iref, nmid, nch2d
    integer, intent(in)      :: nptB, nptSL, npols, nptchB, nordB(nptchB), ixyzB(nptchB+1),  iptB(nptchB)
    real(kind=8), intent(out):: xinterp(3,nptSL), matinterp(3*nptSL,3*nptB)
    real(kind=8) :: tchse(nch2d+1)
    real(kind=8) :: R_B, rmax, rmSL
    integer :: iort

    !! Estimate rmax and R = rad-delta for the slip surface
    R_B = rad; R_SL = rad-delta; rmax = 1.0d0*R_B*Rfa;

    !! Compute tchse
    iort = 1; pars(1) = aY; pars(2) = aZ; pars(3) = rad
    call get_oocyte3d_tchse(iref,nmid,nch2d,tchse) 

    !! Compute points for interpolation 
    call normpoints(nptSL, R_SL, R_B, surfslip, xinterp)

    call interpmat(nch2d,tchse,kchu,iort,npols,pars, & 
    & nptB,nptSL,nptchB,nordB,ixyzB,iptB,rmax,xinterp,matinterp)

end subroutine get_interp_mat

!*******************************************************
! March forward in time 
!******************************************************
subroutine timemarch(delt,rmax,nDB,rho,xinterp,LB_grad,Jmat,init,nMT, & 
            & brink,iter,kchu,iref,nmid,nch2d,nptB,nptSL,nptchB,nordB,ixyzB,iptB)
    use constants
    use omp_lib
    implicit none
    integer, intent(in)       :: brink,iter,nptB,nptSL,nptchB,nordB(nptchB),ixyzB(nptchB+1),iptB(nptchB)
    integer, intent(in)       :: kchu,iref,nmid,nch2d
    real(kind=8), intent(in)  :: delt, rmax, nDB(3,nptB), rho(nptB), xinterp(3,nptSL)
    real(kind=8), intent(in)  :: LB_grad(9*nptB,3*nptSL), Jmat(9*nptB,3*nptB)
    real(kind=8), intent(in)  :: init(3,nptB)
    real(kind=8), intent(out) :: nMT(3,nptB)
    !----------Variables for closure---------------------------------------------------- 
    real(kind=8) :: q11(8),  q22(8), q33(8) 
    real(kind=8) :: r113(8),r223(8),r333(8)
    real(kind=8) :: s1111(8), s2222(8), s1122(8),s1133(8),s2233(8),s3333(8)
    real(kind=8) :: ppN(3,nptB),pppNN(3,nptB),ppppNNN(3,nptB),ppNN(nptB),pppvg(3,nptB)
    real(kind=8) :: P3Bing(3*nptB,9*nptB)
    !-----------------Variables for traction and evolve----------------------------
    real(kind=8) :: Fact(3,nptB),Fcons(3,nptB),Ftot(3,nptB)
    real(kind=8) :: Fslip(3,nptSL)
    real(kind=8) :: ndot(3,nptB),v_grad(3,3,nptB),vgrad_vec(9*nptB)
    real(kind=8) :: dot(nptB), pars(3), tchse(nch2d+1)
    external :: funcurve_oocyte_riemann
    integer :: mt, ifc, iort
    integer :: i, k

    CALL OMP_SET_NUM_THREADS(4);


    !! Compute tchse
    iort = 1; pars(1) = aY; pars(2) = aZ; pars(3) = rad
    call get_oocyte3d_tchse(iref,nmid,nch2d,tchse) 

    !! Assign nMT to init 
    nMT = init;

    !! No conformal
    ifc = 0;

    !! Compute terms for calculating closures
    call momcoeff(q11,q22,q33,           &
            &    r113,r223,r333,         &
            &    s1111,s2222,s3333,s1122,s1133,s2233)

    !! Time marching starts here 
    open(3, file='MTvec.txt', status="old", position="append")
    do k = 1,iter
        !! Find the contractions for ndot from closures
        call closure(nptB,nMT,nDB,  &
        &   q11,q22,q33,        &
        &   r113,r223,r333,     &
        &   s1111,s2222,s3333,s1122,s1133,s2233, &
        &   ppNN,ppN,pppNN,ppppNNN,P3Bing)


        !! Compute active stress
        call activestress(nptB,nMT,Fact)

        !! Call constrain force jump 
        call constrain(nptB,nMT,nDB,ppNN,ppN,pppNN,ppppNNN, Fcons); 

        !! Define the net slip stress: note there is a negative sign coming due to convention 
        do mt = 1,nptB 
            Ftot(:,mt) = c0*(Fact(:,mt) - beta*Fcons(:,mt))*rho(mt);
        end do 

        !! Now perform interpolation to compute the traction on the slip surface points 
        call axissym_fun_interp(3,nch2d,tchse,kchu,funcurve_oocyte_riemann,3,pars,  &
                        rmax,iort,nptSL,xinterp,nptchB,nordB,ixyzB,iptB,nptB,Ftot,  &
                        Fslip)


        !! Obtain the RHS to solve for the velocity gradient vector 
        call vgradrhs(nptB,nptSL,LB_grad,Fslip,vgrad_vec,v_grad)


        !! Obtain velocity gradient from gmres solve 
        call vgradgmres(brink,nptB,Jmat,P3Bing,rho,vgrad_vec,v_grad)


        !! Find the contractions for ndot from closures
        call closure(nptB,nMT,nDB,      &
                &   q11,q22,q33,        &
                &   r113,r223,r333,     &
                &   s1111,s2222,s3333,s1122,s1133,s2233, &
                &   ppNN,ppN,pppNN,ppppNNN,P3Bing)

        !! Compute ppp:vgrad
        call gradeval(nptB,P3Bing,vgrad_vec,pppvg)


        !! Compute ndot 
        call alignrotate(nptB,nMT,nDB,v_grad, & 
            & pppvg,ppNN,ppN,pppNN,ppppNNN,ndot)

        !! Time marching by explicit Euler
        nMT = nMT + ndot*delt;

        if (mod(k-1,10) == 0) then 
            do i = 1,nptB 
                write(3,*)  i, real(nMT(:,i))
            end do 
        end if

        ! do i = 1,nptB 
        !     dot(i) = dot_product(nMT(:,i),nDB(:,i));
        !     if (abs(dot(i)) > 1.0d0) then 
        !         nMT(:,i) = nMT(:,i)/norm2(nMT(:,i));
        !     end if
        ! end do 
        ! if (maxval(abs(dot)) > 2.d0) then 
        !     print*, "error please break: code will stop"
        !     stop
        ! end if

    end do 
    close(3)


end subroutine timemarch
!******************************************************************************************
! Estimate rmax for the slip surface such that both the surfaces have same no. of points
!******************************************************************************************
subroutine rmxestimate(R, nptchB, kchu, iref, nmid, nch2d, rmSL, nptchSL)
    use constants 
    implicit none 
    integer, intent(in)       :: nptchB, kchu,iref,nmid,nch2d
    real(kind=8), intent(in)  :: R
    real(kind=8), intent(inout):: rmSL
    integer,      intent(out):: nptchSL
    real(kind=8) :: incr
    integer :: check, count

    !rmSL = (rad-delta);
    incr = 0.001; 
    call get_oocyte3d_riemann_fun_mem(aY,aZ,R,kchu,iref,nmid,nch2d,rmSL,nptchSL)
    check= nptchB-nptchSL; 


    count = 1;
    do while ((check .ne. 0) .and. (count < 50000))
        if (check > 0) then
            rmSL = rmSL-incr; 
        else
            rmSL = rmSL+incr;
        end if
        call get_oocyte3d_riemann_fun_mem(aY,aZ,R,kchu,iref,nmid,nch2d,rmSL,nptchSL)
        check= nptchB-nptchSL; count = count + 1;
    end do 
   
end subroutine rmxestimate
!******************************************************************************************
! Compute points on the boundary that are normal to the slip surface point 
!******************************************************************************************
subroutine normpoints(nptSL, rSL, rB, surfslip, xinterp)
    use constants 
    implicit none 
    integer, intent(in) :: nptSL
    real(kind=8), intent(in) :: rSL, rB, surfslip(12,nptSL)
    real(kind=8), intent(out):: xinterp(3,nptSL) 
    real(kind=8) :: x0(nptSL), y0(nptSL), z0(nptSL) 
    real(kind=8) :: t(nptSL), dis(nptSL), nDx(nptSL), nDy(nptSL), nDz(nptSL)
    real(kind=8) :: vec(3), nv(2)
    integer :: j

    do j = 1,nptSL 
        x0(j) = surfslip(1,j); y0(j) = surfslip(2,j); z0(j) = surfslip(3,j);
        call testimate(rSL,x0(j),y0(j),z0(j),t(j))
    end do   

    !! Compute the point
    dis = 1.d0;
    do j = 1,nptSL 
        x0(j) = surfslip(1,j); y0(j) = surfslip(2,j); z0(j) = surfslip(3,j);
        vec(1)= surfslip(10,j); 
        vec(2)= surfslip(11,j); 
        vec(3)= surfslip(12,j); 
        
        call curvenormal(rSL,t(j),nv)
        call nraph(rB,x0(j),y0(j),z0(j),nv,dis(j),t(j))


        xinterp(1,j) = x0(j)+dis(j)*vec(1); 
        xinterp(2,j) = y0(j)+dis(j)*vec(2); 
        xinterp(3,j) = z0(j)+dis(j)*vec(3); 
    end do

end subroutine normpoints
!******************************************************************************************
! Compute points on the boundary that are normal to the slip surface point 
!******************************************************************************************
subroutine readprecomp(nptB,nptSL,LB_grad,tracD_mat,pdir)
    implicit none 
    integer, intent(in) :: nptB,nptSL
    real(kind=8), intent(out) ::  LB_grad(9*nptB,3*nptSL), tracD_mat(3*nptB,3*nptSL)
    CHARACTER(len=255), intent(in) :: pdir
    real(kind=8) :: a1(9*nptB*3*nptSL), a2(3*nptB*3*nptSL)
    CHARACTER(len=255) :: pwd


    call getcwd(pwd)
    CALL chdir(pdir)

    !! 1. vel_D2T
    open(5, file="LB_grad.bin", form='unformatted', access='stream')
     read(5) a1
    close(5)
    LB_grad = reshape(a1, (/9*nptB, 3*nptSL/)); 

    !! 2. vel_S2D
    open(5, file="tracD_mat.bin", form='unformatted', access='stream')
     read(5) a2
    close(5)
    tracD_mat = reshape(a2, (/3*nptB, 3*nptSL/)); 

    CALL chdir(pwd)


end subroutine readprecomp
!****************************************************************************************
! Compute surface density by prescribing a dependence on theeta-and-phi
!****************************************************************************************
subroutine density(npts,surfB,type,rho)  
    use constants 
    implicit none 
    integer, intent(in)      :: npts
    real(kind=8), intent(in) :: surfB(12,npts)
    character(len=2), intent(in) :: type
    real(kind=8), intent(out):: rho(npts)
    real(kind=8)             :: xyz(3), theta, phi
    integer :: i
    
    open(2,file='density.txt')
        do i = 1,npts
            xyz = surfB(1:3,i); 
            xyz = xyz/norm2(xyz)
            call cartpol(xyz,theta,phi)

            rho(i) = 1.d0;

            !! Dorso-Ventral gradient
            if (type == 'DV') then 
                rho(i) = 1.0d0 + 0.50d0*sin(theta)*sin(phi);
            end if

            !! Anterior-Posterior gradient
            if (type == 'AP') then 
                rho(i) = 1.0d0 + 0.50d0*cos(theta)
            end if

            write(2,*)  real(rho(i)), real(surfB(1:3,i))
        end do 
    close(2)

end subroutine density 
!****************************************************************************************
! n-field defined on the surface of the container 
!****************************************************************************************
subroutine initial(npts,choice,surfB,nDB,nMT)  
    use constants 
    implicit none 
    integer, intent(in)      :: npts, choice
    real(kind=8), intent(in) :: surfB(12,npts)
    real(kind=8), intent(out):: nDB(3,npts), nMT(3,npts)
    real(kind=8) :: coef1,coef2, fac, dot
    real(kind=8) :: t1(3), t2(3), normal(3)
    real(kind=8) :: vec(3), mat(3,3), tang(3)
    integer :: j,k,mt, idum
    real(kind=8), external :: gasdev

    !! Now provide some tangential perturbations 
    idum = -100; fac = 2.0d0;

    !! Define coefficients for projection
    coef1= 1.0d0; 
    coef2= 1.0d0; 

    vec(1) = 1.d0; vec(2) = 0.d0; vec(3) = 0.d0; ! This is a vector dotted with the tangent plane 

    do mt = 1,npts 
        !! Define coefficients randomly
        coef1 = gasdev(idum); coef2 = gasdev(idum);

        !! Define tangent and normal
        t1 = surfB(4:6,mt); t2 = surfB(7:9,mt); normal = surfB(10:12,mt)
        
        !! Use this for tangent plane projection
        do j = 1,3 
            mat(j,j) = 1.d0;
            do k = 1,3
                mat(j,k) = mat(j,k)-normal(j)*normal(k); 
            end do
            tang(j) = mat(j,1)*vec(1) +  mat(j,2)*vec(2) + mat(j,3)*vec(3); 
        end do 

        !! Define MT orientation and normalize
        if (choice == 0) then 
            nMT(:,mt) = (coef1*t1 + coef2*t2)*fac - 2.d0*normal; 
        else
            nMT(:,mt) = tang - 2.d0*normal; 
        end if
        nMT(:,mt) = tang*0.01 - surfB(10:12,mt);


        dot = norm2(nMT(:,mt))
        nMT(:,mt) = nmT(:,mt)/dot
    end do 

    !! Flip sign of the normal to the surface: make them point inwards
    nDB =-surfB(10:12,:);


end subroutine initial 
!*********************************************************************
! INTIAIL CONDITION ON THE ANGLE
!*********************************************************************
subroutine restart(npts,choice,surfB,nDB,nMT)  
    use constants
    implicit none
    integer, intent(in)      :: npts, choice
    real(kind=8), intent(in) :: surfB(12,npts)
    real(kind=8), intent(out):: nDB(3,npts), nMT(3,npts)
    real(kind=8) :: psiar(npts,3)
    integer :: i,j
    
    ! Define spacing and random seed
    open (12, file="restart.txt")
    read(12,*) ((psiar(i,j), j = 1,3), i = 1,npts)
    close(12)
    
    
    do i = 1,npts
    nMT(:,i)   = psiar(i,:)
    end do

    !! Flip sign of the normal to the surface: make them point inwards
    nDB =-surfB(10:12,:);

    
end subroutine restart
!****************************************************************************************
! Compute active traction of the problem 
!****************************************************************************************
subroutine activestress(N,nMT,Fact)
    use constants
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: nMT(3,N)
    real(kind=8), intent(out):: Fact(3,N) 

    Fact =-sigma*nMT;
 
end subroutine activestress
!****************************************************************************************
! Alignment torque that tries to align the microtubule with the internal normal
!****************************************************************************************
subroutine constrain(N,nMT,nD, & 
    & ppNN,ppN,pppNN,ppppNNN, Fcons)
    use constants
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: nMT(3,N),  nD(3,N)
    real(kind=8), intent(in) :: ppNN(N), ppN(3,N), pppNN(3,N), ppppNNN(3,N)
    real(kind=8), intent(out):: Fcons(3,N) 
    real(kind=8) :: pN(N)
    real(kind=8) :: conx, cony, conz, term
    integer :: i,mt
     
    !! Now compute the contractions of high order moments with surface normal vector 
    pN = 0.d0; 
    do mt = 1,N 
        do i = 1,3 
            pN(mt) = pN(mt) + nD(i,mt)*nMT(i,mt); ! ni nDi 
        end do 
    end do 
    
    !! Compute the constrain stress now
    Fcons = 0.d0;
    do mt = 1,N 
        ! Initialize 
        conx = 0.d0; cony = 0.d0; conz = 0.d0;
        ! Scalar for first term 
        term = pi/2.d0 - pN(mt) + pi/4.d0*ppNN(mt); 
        ! Components
        conx = nD(1,mt)*term - (pi/2.d0*ppN(1,mt) - pppNN(1,mt) + pi/4.d0*ppppNNN(1,mt));
        cony = nD(2,mt)*term - (pi/2.d0*ppN(2,mt) - pppNN(2,mt) + pi/4.d0*ppppNNN(2,mt));
        conz = nD(3,mt)*term - (pi/2.d0*ppN(3,mt) - pppNN(3,mt) + pi/4.d0*ppppNNN(3,mt));
        ! Assign to vector
        Fcons(1,mt) = conx; Fcons(2,mt) = cony; Fcons(3,mt) = conz;
    end do 

end subroutine constrain 
!****************************************************************************************
! Alignment torque that tries to align the microtubule with the internal normal
!****************************************************************************************
subroutine alignrotate(N,nMT,nD,v_grad,& 
    & pppvg,ppNN,ppN,pppNN,ppppNNN, ndot)
    use constants
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: nMT(3,N),  nD(3,N)
    real(kind=8), intent(in) :: v_grad(3,3,N)
    real(kind=8), intent(in) :: pppvg(3,N), ppN(3,N), pppNN(3,N), ppppNNN(3,N)
    real(kind=8), intent(out):: ndot(3,N) 
    real(kind=8) :: pN(N), ppNN(N)
    real(kind=8) :: pvg(3,N)
    real(kind=8) :: nxdot, nydot, nzdot, term
    integer :: i,j,mt
        
        !! Now compute the contractions of high order moments with surface normal vector 
        pN = 0.d0; pvg= 0.d0;
        do mt = 1,N 
            do i = 1,3 
                pN(mt) = pN(mt) + nD(i,mt)*nMT(i,mt); ! ni nDi 
                do j = 1,3
                    pvg(i,mt) = pvg(i,mt) + nMT(j,mt)*v_grad(i,j,mt) ! nk ui/xk
                end do 
            end do 
        end do 
    
        !! Now create the entries for ndot arising due to the torque 
        ! Initialize ndot 
        ndot = 0.d0;
        do mt = 1,N 
        ! Initialize 
        nxdot = 0.d0; nydot = 0.d0; nzdot = 0.d0;
        ! Scalar for first term 
        term = pi/2.d0 - pN(mt) + pi/4.d0*ppNN(mt); 
        ! Components
        nxdot = nD(1,mt)*term - (pi/2.d0*ppN(1,mt) - pppNN(1,mt) + pi/4.d0*ppppNNN(1,mt));
        nydot = nD(2,mt)*term - (pi/2.d0*ppN(2,mt) - pppNN(2,mt) + pi/4.d0*ppppNNN(2,mt));
        nzdot = nD(3,mt)*term - (pi/2.d0*ppN(3,mt) - pppNN(3,mt) + pi/4.d0*ppppNNN(3,mt));
        ! Assign to vector
        ndot(1,mt) = nxdot; ndot(2,mt) = nydot; ndot(3,mt) = nzdot;
        end do 
    
        !! Account for velocity gradient 
        do mt = 1,N 
        ! Initialize 
        nxdot = 0.d0; nydot = 0.d0; nzdot = 0.d0;
        ! Compute the contribution due to velocity gradient
        nxdot = pvg(1,mt)-pppvg(1,mt);
        nydot = pvg(2,mt)-pppvg(2,mt);
        nzdot = pvg(3,mt)-pppvg(3,mt);
        ! Add to vectors 
        ndot(1,mt) = ndot(1,mt) + nxdot;
        ndot(2,mt) = ndot(2,mt) + nydot;
        ndot(3,mt) = ndot(3,mt) + nzdot;
        end do 
    
    !! Add rotational diffusion
    ndot = ndot - 2.d0*drot*nMT;
    
    end subroutine alignrotate 
!****************************************************************************************
! Compute the contractions of the moments
!****************************************************************************************
subroutine closure(N,nMT,nD, &
        &   q11,q22,q33,     &
        &   r113,r223,r333,  &
        &   s1111,s2222,s3333,s1122,s1133,s2233, &
        &   ppNN,ppN,pppNN,ppppNNN,P3Bing)
    use constants
    use omp_lib
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: nMT(3,N),  nD(3,N)
    real(kind=8), intent(in) :: q11(8),  q22(8), q33(8) 
    real(kind=8), intent(in) :: r113(8),r223(8),r333(8)
    real(kind=8), intent(in) :: s1111(8), s2222(8), s3333(8), s1122(8),s1133(8),s2233(8)
    real(kind=8), intent(out):: ppNN(N), ppN(3,N), pppNN(3,N), ppppNNN(3,N), P3Bing(3*N,9*N)
    real(kind=8) :: bingmat(3,9)
    integer      :: mt

    P3Bing = 0.d0;
    !$OMP DO PRIVATE(bingmat)
    do mt = 1,N 
    call pppmat(nMT(:,mt),nD(:,mt),r113,r223,r333,bingmat)
    P3Bing(3*(mt-1)+1:3*(mt-1)+3,9*(mt-1)+1:9*(mt-1)+9) = bingmat;
    end do 
    !$OMP END DO


    !! Compute the contractions using the Bingham closure 
    !$OMP DO 
    do mt = 1,N 
    call contract(nMT(:,mt),nD(:,mt),                & 
            &   q11,q22,q33,                         &
            &   r113,r223,r333,                      &
            &   s1111,s2222,s3333,s1122,s1133,s2233, &
            &   ppNN(mt),ppN(:,mt),pppNN(:,mt),ppppNNN(:,mt))
    end do 
    !$OMP END DO


end subroutine closure
!*******************************************************
! Subroutine for cartesian to polar
!******************************************************
subroutine cartpol(xyz,theta,phi) 
    implicit none 
    real(kind=8), intent(in) :: xyz(3)
    real(kind=8), intent(out):: theta, phi 
    real(kind=8) :: pi
    real(kind=8) :: x, y
  
    !! Define the constant pi
    pi = 4.d0*atan(1.d0);
    
    !! Angle with the z axis
    theta = acos(xyz(3)); 
    
    !! Compute phi = [0 2*pi]
    x = xyz(1); y = xyz(2);
    if (x > 0 .and. y .ge. 0) then 
      phi   = atan(y/x);
    elseif (x > 0 .and. y < 0) then 
      phi   = 2.d0*pi + atan(y/x);
    elseif (x < 0) then 
      phi   = atan(y/x) + pi; 
    end if
  
    if (x == 0 .and. y > 0) then 
      phi = pi/2.d0;
    elseif (x == 0 .and. y < 0) then 
      phi = 3.d0*pi/2.d0; 
    end if
  
end subroutine cartpol
    








