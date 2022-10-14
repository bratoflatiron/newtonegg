!*********************************************************************************************
! Compute the interpolation matrix 
!*********************************************************************************************
subroutine interpmat(nch2d,tchse,kchu,iort,npols,pars, & 
                   & nptB,ntarg,nptchB,nordB,ixyzB,iptB,rmax,xinterp,matint)
    use constants 
    use omp_lib
    implicit none 
    integer, intent(in)        :: nch2d, tchse, kchu, npols, iort
    integer, intent(inout)     :: nptB, ntarg, nptchB, nordB(nptchB), ixyzB(nptchB+1), iptB(nptchB)
    real(kind=8), intent(in)   :: pars(3), rmax, xinterp(3,ntarg)
    real(kind=8), intent(out)  :: matint(3*ntarg,3*nptB)
    real(kind=8) :: uvs_targ(2,ntarg)
    real(kind=8), allocatable :: xmattarg(:)
    external :: funcurve_oocyte_riemann
    integer  :: ipatchtarg(ntarg), ixmattarg(ntarg+1)
    integer  :: ipatch, istart, istart2, npols_use, idim
    integer  :: i, l, lmem

    !! Compute local ellipsoid coordinatre
    call axissym_fun_local_coord_targ(nch2d,tchse,kchu, &
                         &  funcurve_oocyte_riemann,    & 
                         &  3,pars,rmax,iort,ntarg,xinterp,nptchB,nordB,ixyzB,iptB, &
                         &  nptB,ipatchtarg,uvs_targ)


    !! Compute lmem length
    lmem = 0;
    call get_surf_interp_mat_targ_mem(nptchB,ixyzB,ntarg,ipatchtarg,lmem)


    !! Compute xmattarg
    allocate(xmattarg(lmem))
    call get_surf_interp_mat_targ(nptchB,nordB,ixyzB,iptB,nptB,ntarg,ipatchtarg,uvs_targ,lmem,xmattarg,ixmattarg)

    !! Form the interpolation matrix (size 3*nptSL x 3*nptB)
    do i=1,ntarg
        ipatch = ipatchtarg(i)
        istart = ixyzB(ipatch)
        istart2 = ixmattarg(i)
        npols_use = ixyzB(ipatch+1) - ixyzB(ipatch) 
        do idim = 1,3
            do l=1,npols
                matint(3*(i-1)+idim,3*(istart+l-2)+idim) =  xmattarg(istart2+l-1)
            end do 
        end do
    end do

end subroutine interpmat 
! !*******************************************************************************************
! ! EVALUATE GRADIENT USING GMRES 
! !*******************************************************************************************
subroutine vgradgmres(brink,nptB,Jmat,P3Bing,rho,vgrad_vec,vgrad)
    use constants 
    use omp_lib
    implicit none 
    integer, intent(inout)     :: brink
    integer, intent(in)        :: nptB
    real(kind=8), intent(in)   :: Jmat(9*nptB,3*nptB), P3Bing(3*nptB,9*nptB), rho(nptB)
    real(kind=8), intent(inout):: vgrad_vec(9*nptB), vgrad(3,3,nptB) 
    real(kind=8) :: fac, rhs(9*nptB), errs(50), res, soln(9*nptB)
    real(kind=8) :: P3Bing_rho(3*nptB,9*nptB)
    integer :: mt
    integer :: i, maxit, niter
    
    do mt = 1,nptB
        P3Bing_rho(3*(mt-1)+1:3*(mt-1)+3,9*(mt-1)+1:9*(mt-1)+9) = rho(mt)*P3Bing(3*(mt-1)+1:3*(mt-1)+3,9*(mt-1)+1:9*(mt-1)+9);
    end do 

    if (brink == 1) then 
        rhs = vgrad_vec; fac = c0*beta/2.d0;
        
        !! Use GMRES to compute solution
        maxit = 50;
        call dgmres_stok_blas(nptB,fac*Jmat,P3Bing_rho,1.d0,maxit,rhs,1.d-6,niter,errs,res,soln)
        !print*, "iter:", niter, real(res)
        vgrad_vec = soln;

        !! Transpose with negative sign to compute vgrad
        vgrad = 0.d0;
        !$OMP DO
        do i = 1,nptB 
            vgrad(:,:,i) = -transpose(reshape(vgrad_vec(9*(i-1)+1:9*(i-1)+9),(/3,3/)))
        end do 
        !$OMP END DO 

        !! If the GMRES did not converge to right tolerence thn don't perform/include the Brinkman term: brink = 0
        if (res > 1.d-6) then 
            brink = 0; 
        end if 


        open(2,file='param.txt',status="old", position="append")
            write(2,*) "iter:", niter, real(res), "brink:", brink
        close(2)
    end if

end subroutine vgradgmres
! !*******************************************************************************************
! ! EVALUATE GRADIENT CORRECTLY 
! !*******************************************************************************************
subroutine lhsform(nptB,Jmat,P3Bing)!,LHS)
    use omp_lib
    use constants 
    implicit none 
    integer, intent(in)      :: nptB
    real(kind=8), intent(in) :: Jmat(9*nptB,3*nptB), P3Bing(3*nptB,9*nptB)
    !real(kind=8), intent(out):: LHS(9*nptB,9*nptB)
    real(kind=8) :: LHS(9*nptB,9*nptB), test(3,3)
    integer :: i,j,reclen
    real(kind=8) :: ostart, oend

    ! Multiply A*B = C; (A is MxK), (B is KxN), (C is MxN)
    ! dgemm	('No transpose','No transpose',
    !     integer 	M,
    !     integer 	N,
    !     integer 	K,
    !     double precision 	1.0,
    !     double precision, dimension(lda,*) 	A,
    !     integer 	M,
    !     double precision, dimension(ldb,*) 	B,
    !     integer 	K,
    !     double precision 	0.d0,
    !     double precision, dimension(ldc,*) 	C,
    !     integer 	M 
    !     )	

    !! Compute Matrix-Matrix multiplication to get the B matrix: Jmat(9*nptB,3*nptB) x P3Bing(3*nptB,9*nptB)
    call dgemm('N','N',9*nptB,9*nptB,3*nptB,1.d0,Jmat,9*nptB,P3Bing,3*nptB,0.d0,LHS,9*nptB);

    !! Add the identity term 
    LHS = c0*beta*LHS/2.d0;
    do i = 1,9*nptB 
        LHS(i,i) = 1.d0 + LHS(i,i); 
    end do    

    inquire(iolength=reclen) LHS
    open(5,file='tracD_mat.bin',form='unformatted',status='replace',access='stream')
        write(5) real(LHS)
    close(5)
    
    ! test = 0.d0; test(1,3) = 4.d0; test(3,1) = -1;
    ! inquire(iolength=reclen) test
    ! open(5,file='tracD_mat.bin',form='unformatted',access='stream')
    !     write(5) real(test)
    ! close(5)


end subroutine lhsform
! !*******************************************************************************************
! ! EVALUATE GRADIENT CORRECTLY 
! !*******************************************************************************************
subroutine lhsformprec(nptB,Jmat,P3Bing)!,LHS)
    use constants 
    implicit none 
    integer, intent(in)      :: nptB
    real(kind=8), intent(in) :: Jmat(9*nptB,3*nptB), P3Bing(3*nptB,9*nptB)
    !real(kind=8), intent(out):: LHS(9*nptB,9*nptB)
    real(kind=8) :: LHS(9*nptB,9*nptB), test(3,3)
    integer :: i,j,reclen

    ! Multiply A*B = C; (A is MxK), (B is KxN), (C is MxN)
    ! dgemm	('No transpose','No transpose',
    !     integer 	M,
    !     integer 	N,
    !     integer 	K,
    !     double precision 	1.0,
    !     double precision, dimension(lda,*) 	A,
    !     integer 	M,
    !     double precision, dimension(ldb,*) 	B,
    !     integer 	K,
    !     double precision 	0.d0,
    !     double precision, dimension(ldc,*) 	C,
    !     integer 	M 
    !     )	

    !! Compute Matrix-Matrix multiplication to get the B matrix: Jmat(9*nptB,3*nptB) x P3Bing(3*nptB,9*nptB)
    call dgemm('N','N',9*nptB,9*nptB,3*nptB,1.d0,Jmat,9*nptB,P3Bing,3*nptB,0.d0,LHS,9*nptB);

    !! Add the identity term 
    LHS = c0*beta*LHS/2.d0;
    do i = 1,9*nptB 
        LHS(i,i) = 1.d0 + LHS(i,i); 
    end do    

    inquire(iolength=reclen) LHS
    open(5,file='tracD_prec.bin',form='unformatted',status='replace',access='stream')
        write(5) real(LHS)
    close(5)
    
    ! test = 0.d0; test(1,3) = 4.d0; test(3,1) = -1;
    ! inquire(iolength=reclen) test
    ! open(5,file='tracD_mat.bin',form='unformatted',access='stream')
    !     write(5) real(test)
    ! close(5)


end subroutine lhsformprec