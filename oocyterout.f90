      subroutine get_oocyte3d_riemann_chunks(a,b,rad,k,iref,nmid, &
        nch2d,npts2d,iptype2d,ixys2d,srccoefs2d,srcvals2d)
      implicit real *8 (a-h,o-z)
      real *8 t,p1,p2
      integer k,nch2d,npts2d
      integer iptype2d(nch2d),ixys2d(nch2d+1)
      real *8 srccoefs2d(6,npts2d),srcvals2d(8,npts2d)
      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8 pars(3)
      complex *16 zpars
      integer ipars
      real *8, allocatable :: tchse(:)
      
      allocate(tchse(nch2d+1))

      call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)
      
      !call prin2_long('tchuse=*',tchse,nch2d+1)

      itype = 2
      allocate(ts(k),ws(k),umat(k,k),vmat(k,k))
      call legeexps(itype,k,ts,umat,vmat,ws)

      do i=1,nch2d
        ixys2d(i) = (i-1)*k + 1
        iptype2d(i) = 1
      enddo

      ixys2d(nch2d+1) = npts2d+1
      pars(1) = a
      pars(2) = b
      pars(3) = rad
      np = 3
      
      do ich=1,nch2d
        tchstart = tchse(ich) 
        tchend = tchse(ich+1) 

        hh = (tchend-tchstart)/2.0d0

        do j=1,k
          ipt = (ich-1)*k + j
          tt = tchstart + (ts(j)+1.0d0)/2.0d0*(tchend-tchstart)

          call funcurve_oocyte_riemann(tt,np,pars,x,y,dx,dy,dx2,dy2)
          srcvals2d(1,ipt) = x
          srcvals2d(2,ipt) = y
          srcvals2d(3,ipt) = dx*hh
          srcvals2d(4,ipt) = dy*hh
          srcvals2d(5,ipt) = dx2*hh**2
          srcvals2d(6,ipt) = dy2*hh**2
          ds = sqrt(dx**2 + dy**2)
          srcvals2d(7,ipt) = dy/ds
          srcvals2d(8,ipt) = -dx/ds
        enddo
        do j=1,k
          jpt = (ich-1)*k + j
          srccoefs2d(1:6,jpt) = 0
          do l=1,k
            lpt = (ich-1)*k + l
            srccoefs2d(1:6,jpt) = srccoefs2d(1:6,jpt) + &
               umat(j,l)*srcvals2d(1:6,lpt)
          enddo
        enddo
      enddo


      return
      end subroutine get_oocyte3d_riemann_chunks
!
!


      subroutine get_oocyte3d_tchse(iref,nmid,nch2d,tchse)
      implicit real *8 (a-h,o-z)
      integer iref,nmid,nch2d
      real *8 tchse(nch2d+1)


      tchse(1) = 0.0d0
      hpan = 1.0d0/(nmid+2)
      rpan = hpan*(0.5d0)**(iref)
      tchse(2) = tchse(1) + rpan
      do i=2,iref+1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan*2
      enddo

      do i=iref+2,iref+1+nmid
        tchse(i+1) = tchse(i) + hpan
      enddo

      rpan = hpan/2
      do i=iref+2+nmid,nch2d-1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan/2
      enddo
      tchse(nch2d+1) = 1.0d0


      return
      end subroutine get_oocyte3d_tchse



      subroutine funcurve_oocyte_riemann(t,np,pars,x,y,dxdt, &
        dydt,d2xdt2,d2ydt2)
      use constants
      implicit real *8 (a-h,o-z)
      real *8 t,pars(np)

      a    = pars(1); b = pars(2); r = pars(3); 
      area = pi*r**2; 
      w    = sqrt(area/pi + a**2 + b**2)

      tt     = pi*t
      x      =        w*sin(tt) - a*sin(tt) - b/sqrt(2.d0)*sin(2*tt);            x = x/Rs;
      dxdt   =    pi*(w*cos(tt) - a*cos(tt) - b*sqrt(2.d0)*cos(2*tt));        dxdt = dxdt/Rs;
      d2xdt2 =-pi*pi*(w*sin(tt) - a*sin(tt) - b*sqrt(2.d0)*2.d0*sin(2*tt)); d2xdt2 = d2xdt2/Rs;
      
      y      =        w*cos(tt) + a*cos(tt) + b/sqrt(2.d0)*cos(2*tt);           y = y/Zs;
      dydt   =   -pi*(w*sin(tt) + a*sin(tt) + b*sqrt(2.d0)*sin(2*tt));       dydt = dydt/Zs;
      d2ydt2 =-pi*pi*(w*cos(tt) + a*cos(tt) + b*sqrt(2.d0)*2.d0*cos(2*tt)); d2ydt2 = d2ydt2/Zs; 

      return
      end subroutine funcurve_oocyte_riemann
      
!
!
!
!
!


      subroutine get_oocyte3d_riemann_fun_mem(a,b,rad,k,iref,nmid,nch2d, &
         rmax,npatches)
!
!  This subroutine estimates the number of patches required
!  for parametrizing an oocyte like geometry parameterized via
!  the parameters a,b
!
!  r(t) = w/2*sin(pi*t) - a/2*sin(pi*t) - b/2/sqrt(2)*sin(2*pi*t)
!  z(t) = w*cos(pi*t) + a*cos(pi*t)  + b/sqrt(2)*cos(2*pi*t),
!
!  t\in [0,1] and w = sqrt(1+a^2 + b^2)
!
!
!
!  The generating curve is discretized via the parameters,
!  iref, nmid, leading to a total number of chunks
!  2*(iref+1) + nmid. The interval [0,1] is first uniformly
!  discretized using nmid+2 patches, and the first and last
!  patch are then dyadically refined for iref levels
!
!  Input arguments:
!    - a: real *8
!           parameter a above
!    - b: real *8
!           parameter b above
!    - k: integer
!        order of discretization for chunking the generating curve
!    - iref: integer
!        number of dyadic refinements at either end points
!    - nmid: integer
!        number of panels in the middle section
!    - nch2d: integer
!        number of chunks on the generating curve = 2*(iref+1) + nmid
!    - rmax: real *8
!        max size of patches in the 3d discretization
!  
!  Output arguments:
!    - npatches: integer
!        number of patches in the 3d discretization
!
!
!
      implicit real *8 (a-h,o-z)
      real *8 t,p1,p2,rmax
      integer nch2d,npatches,k
      real *8, allocatable :: tchse(:)
      real *8 pars(3)
      external funcurve_oocyte_riemann

      np = 3
      pars(1) = a
      pars(2) = b
      pars(3) = rad

      allocate(tchse(nch2d+1))
      call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)

      call get_axissym_fun_mem(nch2d,tchse,k,funcurve_oocyte_riemann, &
        np,pars,rmax,npatches)
      
      return
      end subroutine get_oocyte3d_riemann_fun_mem
!
!
!
!
      subroutine get_oocyte3d_riemann_fun_geom(a,b,rad,k,iref,nmid,nch2d, &
         rmax,norder,npatches,npts,norders,ixyzs,iptype,srcvals, &
         srccoefs)
!
!  This subroutine discretizes an oocyte like geometry parameterized via
!  the parameters a,b
!
!  r(t) = w/2*sin(pi*t) - a/2*sin(pi*t) - b/2/sqrt(2)*sin(2*pi*t)
!  z(t) = w*cos(pi*t) + a*cos(pi*t)  + b/sqrt(2)*cos(2*pi*t),
!
!  t\in [0,1] and w = sqrt(1+a^2 + b^2)
!
!  The generating curve is discretized via the parameters,
!  iref, nmid, leading to a total number of chunks
!  2*(iref+1) + nmid. The interval [0,1] is first uniformly
!  discretized using nmid+2 patches, and the first and last
!  patch are then dyadically refined for iref levels
!
!  Input arguments:
!    - a: real *8
!           parameter a above
!    - b: real *8
!           parameter b above
!    - k: integer
!        order of discretization for chunking the generating curve
!    - iref: integer
!        number of dyadic refinements at either end points
!    - nmid: integer
!        number of panels in the middle section
!    - nch2d: integer
!        number of chunks on the generating curve = 2*(iref+1) + nmid
!    - rmax: real *8
!        max size of patches in the 3d discretization
!  
!  Output arguments:
!    - npatches: integer
!        number of patches in the 3d discretization
!
!
!
      implicit real *8 (a-h,o-z)
      real *8 t,p1,p2,rmax
      integer nch2d,npatches,k
      integer npts,norder
      integer norders(npatches),ixyzs(npatches+1),iptype(npatches)
      real *8 srcvals(12,npts),srccoefs(9,npts)

      real *8 pars(3)
      real *8, allocatable :: tchse(:)
      external funcurve_oocyte_riemann
      
      np = 3
      pars(1) = a; pars(2) = b; pars(3) = rad;
      allocate(tchse(nch2d+1))
      call get_oocyte3d_tchse(iref,nmid,nch2d,tchse)
      
      iort = 1
      call get_axissym_fun_geom(nch2d,tchse,k,funcurve_oocyte_riemann, &
        np,pars,rmax,iort,norder,npatches,npts,norders,ixyzs,iptype, &
        srcvals,srccoefs)
      
      return
      end subroutine get_oocyte3d_riemann_fun_geom
