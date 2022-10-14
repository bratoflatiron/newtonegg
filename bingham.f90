!********************************************************
! Chebyshev polynomials 
!********************************************************
subroutine chebyshev(x,chebT)
    implicit none 
    real(kind=8), intent(in) :: x
    real(kind=8), intent(out):: chebT(8)

    chebT(1) = 1.d0; 
    chebT(2) = x; 
    chebT(3) = 2.d0*x**2   - 1.d0;
    chebT(4) = 4.d0*x**3  -  3.d0*x; 
    chebT(5) = 8.d0*x**4  -  8.d0*x**2 + 1.d0; 
    chebT(6) = 16.d0*x**5 - 20.d0*x**3 + 5.d0*x; 
    chebT(7) = 32.d0*x**6 - 48.d0*x**4 + 18.d0*x**2 - 1.d0; 
    chebT(8) = 64.d0*x**7 -112.d0*x**5 + 56.d0*x**3 - 7.d0*x; 

end subroutine chebyshev 
!********************************************************
! Rodrigue's formula 
!********************************************************
subroutine rodrigue(n,W,Wt)
    implicit none 
    real(kind=8), intent(in) :: n(3)
    real(kind=8), intent(out):: W(3,3), Wt(3,3)
    real(kind=8) :: vu(3), ez(3), theta
    real(kind=8) :: kv(3), K1(3,3), K2(3,3), eye(3,3)

    !! Define the unit vector and z axis 
    vu = n/norm2(n)
    ez = 0.d0; ez(3) = 1.d0;

    !! Define the angle of rotation 
    theta = acos(vu(3)); 

    !! Cross product: (vu x ez)
    kv = 0.d0; 
    kv(1) = vu(2); kv(2) =-vu(1); 
    kv = kv/norm2(kv)

    !! Define the cross product matrix 
    K1 = 0.d0; 
    K1(1,2) = -kv(3); K1(1,3) = kv(2);
    K1(2,1) =-K1(1,2);
    K1(2,3) =-kv(1);
    K1(3,1) =-K1(1,3); K1(3,2) = -K1(2,3);

    K2 = matmul(K1,K1)

    !! Define the identity matrix 
    eye = 0.d0; eye(1,1) = 1.d0; eye(2,2) = 1.d0; eye(3,3) = 1.d0; 

    !! Define the final rotation matrix and its transpose 
    W = eye + sin(theta)*K1 + (1-cos(theta))*K2;
    Wt= transpose(W)

end subroutine rodrigue 
!********************************************************
! Compute the Q-tensor from the Bingham closure 
!********************************************************
subroutine Qtensor(nMT,q11,q22,q33,Q)
    implicit none 
    real(kind=8), intent(in) :: nMT(3)
    real(kind=8), intent(in) :: q11(8),  q22(8), q33(8) 
    real(kind=8), intent(out):: Q(3,3)
    real(kind=8) :: W(3,3), Wt(3,3), QR(3,3)
    real(kind=8) :: nmag, cMT, chebT(8)
    real(kind=8) :: qn11,  qn22, qn33 

    !! Obtain rodrigue's rotation matrices 
    call rodrigue(nMT,W,Wt)

    !! Now compute the chebyshev input 
    nmag = norm2(nMT); 
    cMT  = 2.d0*nmag-1.d0; 

    !! Compute Chebyshev polynomial interpolant value 
    call chebyshev(cMT,chebT)

    !! Compute the moments 
    call sumcoeff(chebT,q11,qn11)
    call sumcoeff(chebT,q22,qn22)
    call sumcoeff(chebT,q33,qn33)

    !! Compute the Q-tensor in the rotated frame 
    QR = 0.d0; 
    QR(1,1) = qn11; QR(2,2) = qn22; QR(3,3) = qn33; 

    !! Now Rotate it to the physical frame
    Q = matmul(QR,W); 
    Q = matmul(Wt,Q);

end subroutine Qtensor
!********************************************************
! Contractions with unit normal vectors on the surface  
!********************************************************
subroutine contract(nMT,nD,         & 
                &   q11,q22,q33,    &
                &   r113,r223,r333, &
                &   s1111,s2222,s3333,s1122,s1133,s2233, &
                &   T2,T3,T4,T5)
    implicit none 
    real(kind=8), intent(in) :: nMT(3), nD(3)
    real(kind=8), intent(in) :: q11(8),  q22(8), q33(8) 
    real(kind=8), intent(in) :: r113(8),r223(8),r333(8)
    real(kind=8), intent(in) :: s1111(8), s2222(8), s3333(8), s1122(8),s1133(8),s2233(8)
    real(kind=8), intent(out):: T2, T3(3), T4(3), T5(3) 
    real(kind=8) :: nr(3), W(3,3), Wt(3,3)
    real(kind=8) :: nmag, cMT, chebT(8)
    real(kind=8) :: qn11,  qn22, qn33 
    real(kind=8) :: rn113,rn223,rn333
    real(kind=8) :: sn1111, sn2222, sn3333, sn1122, sn1133, sn2233

    !! Obtain rodrigue's rotation matrices 
    call rodrigue(nMT,W,Wt)

    !! Rotate the surface normal 
    nr = matmul(W,nD)

    !! Now compute the chebyshev input 
    nmag = norm2(nMT); 
    cMT  = 2.d0*nmag-1.d0; 

    !! Compute Chebyshev polynomial interpolant value 
    call chebyshev(cMT,chebT)

    !! Compute the moments 
    call sumcoeff(chebT,q11,qn11)
    call sumcoeff(chebT,q22,qn22)
    call sumcoeff(chebT,q33,qn33)

    call sumcoeff(chebT,r113,rn113)
    call sumcoeff(chebT,r223,rn223)
    call sumcoeff(chebT,r333,rn333)

    call sumcoeff(chebT,s1111,sn1111)
    call sumcoeff(chebT,s2222,sn2222)
    call sumcoeff(chebT,s3333,sn3333)
    call sumcoeff(chebT,s1122,sn1122)
    call sumcoeff(chebT,s1133,sn1133)
    call sumcoeff(chebT,s2233,sn2233)



    !! Compute the relevant contractions 
    ! pp:nd nd
    T2 = qn11*nr(1)**2 +  qn22*nr(2)**2 + qn33*nr(3)**2;
    
    ! pp:nd
    T3 = 0.d0; 
    T3(1) = qn11*nr(1); 
    T3(2) = qn22*nr(2); 
    T3(3) = qn33*nr(3); 
    T3 = matmul(Wt,T3)
    
    ! ppp:nd nd
    T4 = 0.d0; 
    T4(1) = 2*rn113*nr(1)*nr(3);
    T4(2) = 2*rn223*nr(2)*nr(3);
    T4(3) = rn113*nr(1)**2 + rn223*nr(2)**2 + rn333*nr(3)**2;
    T4 = matmul(Wt,T4)
   
    ! pppp:nd nd nd 
    T5 = 0.d0; 
    T5(1) = sn1111*nr(1)**3 + 3*sn1122*nr(1)*nr(2)**2 + 3*sn2233*nr(1)*nr(3)**2;
    T5(2) = sn2222*nr(2)**3 + 3*sn1122*nr(2)*nr(1)**2 + 3*sn2233*nr(2)*nr(3)**2;
    T5(3) = sn3333*nr(3)**3 + 3*sn2233*nr(3)*nr(1)**2 + 3*sn2233*nr(3)*nr(2)**2;
    T5 = matmul(Wt,T5)

end subroutine contract
!**********************************************************************
! Obtain matrix representation for the Bingham matrix for calculation
!**********************************************************************
subroutine pppmat(nMT,nD,r113,r223,r333,bingmat)
    implicit none 
    real(kind=8), intent(in) :: nMT(3), nD(3)
    real(kind=8), intent(in) :: r113(8),r223(8),r333(8)
    real(kind=8), intent(out):: bingmat(3,9)
    real(kind=8) :: nr(3), W(3,3), Wt(3,3)
    real(kind=8) :: nmag, cMT, chebT(8)
    real(kind=8) :: rn113,rn223,rn333
    real(kind=8) :: R(3,3,3), R_rot(3,3,3), t(9,1)
    integer :: i,j,k

    !! Obtain rodrigue's rotation matrices 
    call rodrigue(nMT,W,Wt)

    !! Rotate the surface normal 
    nr = matmul(W,nD)

    !! Now compute the chebyshev input 
    nmag = norm2(nMT); 
    cMT  = 2.d0*nmag-1.d0; 

    !! Compute Chebyshev polynomial interpolant value 
    call chebyshev(cMT,chebT)

    !! Compute the moments 
    call sumcoeff(chebT,r113,rn113)
    call sumcoeff(chebT,r223,rn223)
    call sumcoeff(chebT,r333,rn333)

    !! Compute the Bingham matrix for mapping in the rotated coordinate
    R = 0.d0; 
    do i = 1,3
        if (i == 1) then 
            do j = 1,3 
                do k = 1,3
                    R(i,k,j) = rn113*(W(1,j)*Wt(k,3)+W(3,j)*Wt(k,1));
                end do 
            end do 
        elseif (i == 2) then 
            do j = 1,3 
                do k = 1,3
                    R(i,k,j) = rn223*(W(2,j)*Wt(k,3)+W(3,j)*Wt(k,2));
                end do 
            end do 
        elseif (i == 3) then 
            do j = 1,3 
                do k = 1,3
                    R(i,k,j) = rn333*(W(3,j)*Wt(k,3))+rn113*(W(1,j)*Wt(k,1))+rn223*(W(2,j)*Wt(k,2));
                end do 
            end do 
        end if 
    end do 

    !! Rotate back to physical space
    R_rot = R; R = R*0;
    do i = 1,3
        do j = 1,3
            R(i,:,:) = R(i,:,:) +  wt(i,j)*R_rot(j,:,:)
        end do 
    end do 
    
    !! Reshape data to obtain output
    do i = 1,3
        t  =  reshape(R(i,:,:),(/9,1/));
        bingmat(i,:) = t(:,1); 
    end do

end subroutine pppmat






