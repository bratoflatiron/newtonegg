!********************************************************************************************
! Subroutine for LU decomposition from NR
! The subroutine takes the following inputs
! a = matrix NXN
! n = size of the matrix a
! np = physical dimension that are used for memory allocation, 500 for 's code style
! indx = this is an integer array that records the permutation order of the LU subst
! d = it takes the value of either +- 1 depending on the odd-even permutation
!*********************************************************************************************
subroutine ludcmp(a,n,np,indx,d)
    implicit none
    integer, intent(in) :: n,np
    integer, intent(out) :: indx(n)
    integer :: NMAX
    real(kind=8), intent(inout) :: a(np,np)
    real(kind=8), intent(out) :: d
    real(kind=8) :: TINY
    parameter (NMAX=10000,TINY=1.d-20)
    integer   :: i,imax,j,k
    real(kind=8) :: aamax,dum,sum,vv(NMAX)

    ! Define d
    d=1.d0
    
    ! Check for singularity
    do i=1,n
        aamax=0.d0
        do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        !if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
    enddo

    do j=1,n
        do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
        end do
        aamax=0.d0
        do i=j,n
            sum=a(i,j)
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
                imax=i
                aamax=dum
            end if
        enddo
        if (j.ne.imax)then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
        end if
        indx(j)=imax
        if(a(j,j).eq.0.d0) a(j,j)=TINY
        if(j.ne.n) then
            dum=1.d0/a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            enddo
        end if
    enddo

end subroutine ludcmp
!*************************************************************************************
! SUBROUTINE for LU backsubstitution
! This routine has to be use in conjunction with ludcmp
! a = the input is the matrix from output of lubksb
! n = size of the matrix a
! np = physical dimension that are used for memory allocation, 500 for 's code style
! indx = this is an integer array that records the permutation order of the LU subst
! b = This is input/output of the subroutine. b comes as the input of the RHS of the eqn Ax = b.
!     Then it gets modified and the solution x of the equation is stored in b
!*************************************************************************************
subroutine lubksb(a,n,np,indx,b)
    integer, intent(in) :: n,np,indx(n)
    real(kind=8), intent(in) :: a(np,np)
    real(kind=8), intent(inout) ::  b(n)
    real(kind=8) :: sum
    integer :: i,ii,j,ll

    ii=0
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            enddo
        elseif (sum.ne.(0.d0)) then
            ii=i
        endif
        b(i)=sum
    enddo

    do i=n,1,-1
        sum=b(i)
        do j=i+1,n
            sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
    enddo

end subroutine lubksb
! ***********************************************************************
! * GASDEV returns double precision figures from a Gaussian Distribution
! * with mean 0 and variance 1. 'IDUM' is set to a negative number to
! * initialize and then not changed. Reset IDUM every time, with a new
! * negative number (iteration number X -1, for instance) to avoid
! * repetition of saem sequence
! ***********************************************************************
 real(kind=8) function gasdev(idum)
 integer, intent(inout):: idum
 integer :: iset
 real(kind=8) :: fac,gset,rsq,v1,v2
real(kind=8), external :: ran1

 !common/rndm/idum
 SAVE iset,gset
 DATA iset/0/
 if (idum.lt.0) iset=0
 if (iset.eq.0) then
1	 	v1=2.d0*ran1(idum)-1.d0
     v2=2.d0*ran1(idum)-1.d0
     rsq=v1**2.d0+v2**2.d0
     if(rsq>=1.d0 .or. rsq==0.d0) go to 1
     fac=sqrt(-2.d0*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
 else
     gasdev=gset
     iset=0
 endif
 return
 END
!--------------------------------------------------------------------------------------
 real(kind=8) function ran1(idum)
 !common/rndm/idum
 INTEGER, intent(inout):: idum
 integer, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB
real(kind=8), parameter ::  AM=1.d0/IM,EPS=1.2d-7,RNMX=1.d0-EPS
 INTEGER :: j,k,iv(NTAB),iy
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
 if (idum<=0 .or. iy==0) then
     idum=max(-idum,1)
     do j=NTAB+8,1,-1
         k=idum/IQ
         idum=IA*(idum-k*IQ)-IR*k
         if (idum<0) idum=idum+IM
         if (j<=NTAB) iv(j)=idum
     enddo
     iy=iv(1)
 endif
 k=idum/IQ
 idum=IA*(idum-k*IQ)-IR*k
 if(idum<0) idum=idum+IM
 j=1+iy/NDIV
 iy=iv(j)
 iv(j)=idum
 ran1=min(AM*iy,RNMX)
 return
 END