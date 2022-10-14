!******************************************************************************************
! Compute the parameter value
!******************************************************************************************
subroutine testimate(rSL,x0,y0,z0,tpar)
    use constants 
    implicit none 
    real(kind=8), intent(in) :: rSL,x0,y0,z0
    real(kind=8), intent(out):: tpar 
    real(kind=8) :: area, w, rtest
    real(kind=8) :: tt(2500), r(2500), z(2500), arr(2500)
    real(kind=8) :: r2, dr, err, v, d(2)
    integer :: i,N, loc
    
    !! Define the tval 
    N = 2500;
    do i = 1,N 
        tt(i) = (i-1.d0)*pi/(N-1.d0); 
    end do 
    
    !! Compute the generating curve with high accuracy
    area = pi*rSL**2; 
    w  =   sqrt(area/pi + aY**2 + aZ**2)
    r  =   w*sin(tt) - aY*sin(tt) - aZ/sqrt(2.d0)*sin(2*tt);     
    z  =   w*cos(tt) + aY*cos(tt) + aZ/sqrt(2.d0)*cos(2*tt);   
    r  = r/Rs; z = z/Zs; 

    !! Compute tpar as a minimum location as it will get rectified 
    rtest = sqrt(x0**2+y0**2); 
    do i = 1,2500 
        d(1) = r(i)-rtest; 
        d(2) = z(i)-z0;
        arr(i) = norm2(d);
    end do 
    loc   = minloc(abs(arr),1);
    tpar  = tt(loc); 
    
    !! Do a newton raphson
    err = 1; 
    do while (err > 1e-8)
        r2 = w*sin(tpar) - aY*sin(tpar) - aZ/sqrt(2.d0)*sin(2*tpar);   r2 = r2/rS;
        dr = w*cos(tpar) - aY*cos(tpar) - 2*aZ*cos(2*tpar)/sqrt(2.d0); dr = dr/rS;
        v  = r2-rtest;
        tpar = tpar - v/dr;
        err = abs(v);
    end do 

end subroutine testimate
!******************************************************************************************
! Compute normal vector 
!******************************************************************************************
subroutine curvenormal(rSL,theta,nn)
    use constants 
    implicit none 
    real(kind=8), intent(in) :: rSL,theta
    real(kind=8), intent(out):: nn(2)
    real(kind=8) :: area, w, f, fp,dst,d2st,xt,yt,xt2,yt2
    real(kind=8) :: a0,a1,a2 

    !! Compute the generating curve with high accuracy
    area = pi*rSL**2; 
    w  =   sqrt(area/pi + aY**2 + aZ**2)

    a0 = w; a1 = aY; a2 = aZ;  

    !! Now calculate all the required derivatives to calculate the normal
    f   = a0**2 + (a1)**2 + 2*a2**2 - 2*a0*(a1)*cos(2*theta) + 2*sqrt(2.d0)*a2*(a1*cos(theta) - a0*cos(3*theta)); f = sqrt(f);
    fp  = 4*a0*(a1)*sin(2*theta) + 2*sqrt(2.d0)*a2*(-a1*sin(theta)  + 3*a0*sin(3*theta)); fp = fp/(2*f);
    dst = 1/f;
    d2st=-fp/f**3;

    xt = -a0*sin(theta) - a1*sin(theta) -   sqrt(2.d0)*a2*sin(2*theta);
    yt =  a0*cos(theta) - a1*cos(theta) -   sqrt(2.d0)*a2*cos(2*theta);
    xt2= -a0*cos(theta) - a1*cos(theta) - 2*sqrt(2.d0)*a2*cos(2*theta);
    yt2= -a0*sin(theta) + a1*sin(theta) + 2*sqrt(2.d0)*a2*sin(2*theta);

    !! Now finally evaluate the normal (arc length derivative)
    nn(1) = xt2*dst**2 + xt*d2st;
    nn(2) = yt2*dst**2 + yt*d2st;
    nn    = -nn/norm2(nn);

end subroutine curvenormal
!******************************************************************************************
! Compute the actual point
!******************************************************************************************
subroutine nraph(rbound,x0,y0,z0,nv,dis,tt)
    use constants 
    implicit none 
    real(kind=8), intent(in)   :: rbound,x0,y0,z0,nv(2)
    real(kind=8), intent(inout):: dis,tt 
    real(kind=8) :: area, w
    real(kind=8) :: r, z, dr, dz
    real(kind=8) :: xB, yB, v(2)
    real(kind=8) :: Amat(2,2), Ainv(2,2), det
    real(kind=8) :: zB, rB, Jd
    real(kind=8) :: guess(2), err

    !! Compute the generating curve with high accuracy
    area = pi*rbound**2; 
    w  =   sqrt(area/pi + aY**2 + aZ**2)

    err= 1.d0; guess(1) = dis; guess(2) = tt;


    do while (err > 1e-5) 
        r  =   w*sin(tt) - aY*sin(tt) - aZ/sqrt(2.d0)*sin(2*tt); r = r/Rs;    
        z  =   w*cos(tt) + aY*cos(tt) + aZ/sqrt(2.d0)*cos(2*tt); z = z/Zs;   
        rB =   sqrt(x0**2+y0**2); zB = z0;
         

        !! Compute derivatives 
        dr = w*cos(tt) - aY*cos(tt) - 2*aZ*cos(2*tt)/sqrt(2.d0); dr = dr/rS;
        dz =-w*sin(tt) - aY*sin(tt) - 2*aZ*sin(2*tt)/sqrt(2.d0); dz = dz/zS;

        !! Define minimizing vector 
        v(1) = dis*nv(1)+zB-z; v(2) = dis*nv(2)+rB-r;
         
        !! Define the Jacobian matrix 
        Amat(1,1) = nv(1); Amat(1,2) = -dz; 
        Amat(2,1) = nv(2); Amat(2,2) = -dr;

        det = Amat(1,1)*Amat(2,2)-Amat(1,2)*Amat(2,1);

        !! Compute inverse of the Jacobian 
        Ainv(1,1) = Amat(2,2); Ainv(1,2) =-Amat(1,2); 
        Ainv(2,1) =-Amat(2,1); Ainv(2,2) = Amat(1,1);
        Ainv = Ainv/det;
        
        !! Calculate new guess 
        guess = guess - matmul(Ainv,v); 

        !! Define error 
        err = norm2(v); 

        !! Update things 
        dis = guess(1); tt = guess(2);

    end do 
    !print*, dis, real(tt)
    ! if (dis > 0.d0) then 
    !     print*, real(dis), real(tt)
    ! end if



end subroutine nraph