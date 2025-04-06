      module global
        use inrtype, only: sp
        implicit none

        integer, parameter :: kr=4,kz=4,kp=4
        real(sp), allocatable :: knotr(:),knotz(:),knotx(:),knoty(:)

        real(sp), allocatable :: bscoefscal1(:,:),bscoefscal2(:,:),bscoefscal3(:,:),bscoefscal4(:,:), &
                                 bscoefpsi(:,:),bscoefbphi(:,:),bscoefscalb(:,:),&
                                 bscoefxy2R(:,:),bscoefxy2Z(:,:),bscoeftheta(:,:)
        real(sp), allocatable :: bscoefvec1r(:,:),bscoefvec1z(:,:),bscoefvec1p(:,:), &
                                bscoefvec2r(:,:),bscoefvec2z(:,:),bscoefvec2p(:,:), &
                                bscoefvec3r(:,:),bscoefvec3z(:,:),bscoefvec3p(:,:), &
                                bscoefvec4r(:,:),bscoefvec4z(:,:),bscoefvec4p(:,:), &
                                bscoefvec5r(:,:),bscoefvec5z(:,:),bscoefvec5p(:,:)
        ! physics constant
        !real(sp), parameter :: mp = 1.6726485d-27, me=9.1094d-31,ee =
        !1.6021766208d-19, mu0 = 4._dp*pi*1.d-7
      end module global
!
!
      module vars1
         use inrtype, only: sp
         implicit none
         integer :: nxefit, nyefit,nxefit1
         real(sp), allocatable :: fold(:,:),gridpsi(:),gridr(:),gridz(:),gridpsi1(:),&
                                  gridx(:),gridy(:),gridx1(:)
         
         real(sp)  :: xdim,zdim,rcentr,rgrid1,zmid,rmagx,zmagx,simagx,sibdry,bcentr,cpasma,xdum
         real(sp), allocatable :: fpol(:), pres(:),workk1(:),workk2(:),qpsi(:),qpsi1(:)
         real(sp) :: signj
         real(sp) :: tmpfpol,tmpq
         real(sp) :: Rmax, Rmin, Zmax, Zmin, Rmax1, Rmin1, Zmax1, Zmin1
         integer :: nrmax, nrmin, nzmax, nzmin
         real(sp), allocatable :: vpsi(:)
         real(sp) :: dr, dz
         !lines(:,i,j)记录第j条等高线上第i个点的r,z,deltapsi,b^2值
         !theta(:,i,j)表示第j条等高线上，第i个点的theta,r,z的值。
         real(sp), allocatable ::  lines(:,:,:), theta(:,:,:)
         real(sp), allocatable :: knotpsi(:)
         real(sp), allocatable ::cscoeffpol(:,:),cscoefq(:,:),&
                                 cscoefrpsi(:,:)
         real(sp), allocatable :: cscoefZ(:,:),cscoefR(:,:),&
                                  cscoefphi(:,:),cscoeftheta(:,:),&
                                  cscoefpsir(:,:),cscoeftime(:,:),&
                                  cscoefvparal(:,:)
         real(sp), allocatable :: cscoefZ1(:,:),cscoefR1(:,:),&
                                  cscoefphi1(:,:),cscoeftheta1(:,:),&
                                  cscoefpsir1(:,:),cscoeftime1(:,:)

         real(sp), allocatable :: scal_b(:,:), scal2(:,:)
         real(sp), allocatable :: vec_b_r(:,:),vec_b_z(:,:),vec_b_p(:,:),&
                                  vec_gradb_r(:,:),vec_gradb_z(:,:),vec_gradb_p(:,:), &
                                  vec_gradbr_r(:,:),vec_gradbr_z(:,:),vec_gradbr_p(:,:),&
                                  vec_gradbz_r(:,:),vec_gradbz_z(:,:),vec_gradbz_p(:,:), &
                                  vec_gradrbp_r(:,:),vec_gradrbp_z(:,:),vec_gradrbp_p(:,:) 
        real(sp) rmag,bmag,vmag,tmag,rminor
        ! rminor: minor radius of wall
        real(sp), allocatable :: rzbdry(:,:),rzlim(:,:)
!
        real(sp) :: rhom_r
        real(sp), allocatable ::vpsir(:),vrpsi(:),smallr(:,:),Rgeom(:),r_rho(:)
        real(sp), allocatable ::knotrpsi(:),knottheta(:),bscoefrr(:,:),bscoefzz(:,:),PR(:,:,:),PZ(:,:,:),JJ(:,:,:),thet(:),tpR(:,:),tpZ(:,:),cscoefqpsi(:,:),qq(:),nablapsi2(:,:,:),nablapt(:,:,:)
        real(sp), allocatable :: VRr(:), VZz(:), psiRZ(:,:), BRZ(:,:),gRZ(:,:)
        real(sp) :: Rrmax,Zzmax,Rrmin,Zzmin
        real(sp), allocatable::vectorpsir(:),vectortheta(:),knotRZpsir(:),knotRZtheta(:),cscoef_rR(:,:),cscoef_rZ(:,:)

        real(sp), allocatable :: alpha(:,:),alphadum(:,:)
        !real(sp) :: La,Lb,Ea,Eb,E0,Ed,Ec,L0,Ld,rr0,rrd
        !integer numL,nE,nt,pa,pb, nrho
!
        integer nbdry,nlim
        integer,parameter :: npsi=122,ntheta=361,ii=10
        character(len=50) gfile
        character(len=50) profiles
        character(len=50) profile
        character(len=100) str
        logical alive
        integer nl,nr,nl2,nr2,nz
        !real(sp) :: rhoh, ne0, v_A0
        !real(sp) :: Ec32,E32,tE,tL,Cn

      end module vars1
!
      module orbit_eq_mod
         use splines
         use global
         use inrtype
         use vars
         use vars1
         use shared_mod
!
         implicit none
         !real(kind = 8) :: alpha, beta, gamma, n1, n2, sigma
         real(kind = 8) :: ctime0, ctime1
         integer, parameter, private :: r8 = selected_real_kind(12,100)
         save
         integer, parameter :: nvar = 6, nmax = 10,  num = 2e4, na = 1
         integer :: nstep
         real(kind=8), parameter :: pi = 3.14159265358979
         ! alpha : kx1/kz1, beta : kz2/kz1, gamma : kx2/kz1
         ! n1: kz1/kzbar, n2: kz2/kzbar, here kzbar = kz1

         ! num : maximum number of particles computed


         contains
           ! Runge Kutta method
           subroutine rk4(y,dydx,n,x,h,yout)
              implicit none
              real(8) :: y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
              real(8) :: hh,h6,xh,x,h
              integer :: i,n
              hh=h*.5
              h6=h/6.0
              xh=x+hh
              do i=1,n
                 yt(i)=y(i)+hh*dydx(i)
              enddo

                 call derivs1(xh,yt,dyt)

              do i=1,n
                 yt(i)=y(i)+hh*dyt(i)
              enddo

                 call derivs1(xh,yt,dym)

              do i=1,n
                  yt(i)=y(i)+h*dym(i)
                  dym(i)=dyt(i)+dym(i)
              enddo

                  call derivs1(x+h,yt,dyt)

              do i=1,n
                 yout(i)=y(i)+h6*(dydx(i)+2.*dym(i)+dyt(i))
              enddo
          end subroutine rk4

          ! Equations of motion, here,
          ! y(1) : R, y(2) : Z, y(3) : phi, y(4): v_paral, 
          ! y(5) : E, y(6) : mube, y(7) ; x, y(8) : y .
          ! dydx(1,2,3,4,5,6,7,8) : derivatives of the above variables with
          ! respect to time
          ! b_nor: B field normalised to B(0) at axis
          subroutine derivs(x,y,dydx)
              implicit none
              !common/aspect/e
              !real(kind=8) :: x,e
              real(kind=8) :: x
              real(kind=8) :: br,bz,bt,b2,b_nor,fac2,fac3,bpz,bpr
              real(kind=8) :: brpz,brpr,bzpr,bzpz,rbtpz,rbtpr,fac4,fac5
              real(kind=8) :: rpx,zpx,rpy,zpy
              real(kind=8) :: y(*),dydx(*)
!              if(y(1)==0.0D0) y(1)=1.0D-a6
              br = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1r,y(1),y(2),0,0)
              bz = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1z,y(1),y(2),0,0)
              bt = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1p,y(1),y(2),0,0)
              b2 = br**2+bz**2+bt**2
              b_nor = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,y(1),y(2),0,0)
              brpz = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec3z,y(1),y(2),0,0)
              brpr = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec3r,y(1),y(2),0,0)
              bzpr = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec4r,y(1),y(2),0,0)
              bzpz = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec4z,y(1),y(2),0,0)
              bpz = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec2z,y(1),y(2),0,0)
              bpr = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec2r,y(1),y(2),0,0)
              rbtpr = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec5r,y(1),y(2),0,0)
              rbtpz = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec5z,y(1),y(2),0,0)
              fac2 = b_nor+y(4)/b2*bt*((brpz-bzpr)-&
                     1.0/b_nor*(br*bpz-bz*bpr+&
                     bz*bpr-br*bpz))  
              fac3 = 1/fac2 
              fac4 = 1/b_nor
              fac5 = 1/y(1)
!
              dydx(1) = fac3*y(4)*&
                        (br+y(4)/b2*bt*bpz)+&
                        fac3*y(6)*y(5)*fac4*bt*bpz
              dydx(2) = fac3*y(4)*&
                        (bz-y(4)/b2*bt*bpr)-&
                        fac3*y(6)*y(5)*fac4*bt*bpr
              dydx(3) = fac3*y(4)*(bt*fac5+&
                        y(4)*fac4*fac5*(brpz-bzpr-&
                        fac4*(br*bpz-bz*bpr)))+&
                        fac3*y(6)*y(5)*fac4*fac5*&
                        (bz*bpr-br*bpz)
              dydx(4) = -fac3*y(6)*y(5)*(br*bpr+bz*bpz)
              dydx(5) = 0.0
              dydx(6) = 0.0
                 return
          end subroutine derivs
!
         ! Equations of motion, here,
          ! y(1) : x, y(2) : theta, y(3) : phi, y(4): v_paral, 
          ! y(5) : E, y(6) : mube .
          ! dydx(1,2,3,4,5,6) : derivatives of the above variables with
          ! respect to time
          ! b_nor: B field normalised to B(0) at axis
!
          subroutine derivs1(x,y,dydx)
              implicit none
              !common/aspect/e
              !real(kind=8) :: x,e
              real(kind=8) :: x
              real(kind=8) :: bb,bb1,bb2,jj1,jj2,qq,bpt,bpx
              real(kind=8) :: rr,rr1,rpx,rpt,jj3
              real(kind=8) :: gg,gg0,gg1,gg2,gg3
              real(kind=8) :: fac1,fac2,fac3,dum,dum1
              real(kind=8) :: y(*),dydx(*)
!              if(y(1)==0.0D0) y(1)=1.0D-a6
              rr=rf(y(2),y(1))
              rr1=rr**2
              !rpx=rfpr(y(2),y(1))
              rpx=rfpx(y(2),y(1))
              rpt=rfpt(y(2),y(1))
              jj1=jf(y(2),y(1))
              jj2=jr2(y(1))
              !jj2=jr20(y(2),y(1))
              jj3=jr2p(y(1))
              !jj3=jr20p(y(2),y(1))
              bb=b_nor(y(2),y(1))
              bb1=b_star(y(2),y(1),y(4))
              bb2=1.0/bb
              bpt=b_norpt(y(2),y(1))
              bpx=b_norpx(y(2),y(1))
              qq = qfun(y(1))
              gg = gzz(y(2),y(1))
              gg0 = grt(y(2),y(1))
              gg1 = grr(y(2),y(1))
              gg2 = gtrpt(y(2),y(1))*jj1+jfpt(y(2),y(1))*grt(y(2),y(1))
              gg3 = grrpx(y(2),y(1))*jj1+jfpx(y(2),y(1))*grr(y(2),y(1))
              fac1 = 1.0/bb**2
              fac2 = 1.0/bb1
              dum = y(4)**2 + y(6)*y(5)*bb
              dum1 = jj3/qq-jj2/qq/qq*dqfun(y(1))
!
              dydx(1) = -fac1*fac2/jj1/e/e*dum*bpt
              dydx(2) = fac2*y(4)/qq/rr1+&
                        fac1*fac2/jj1/e/e*dum*bpx
              dydx(3) = fac2*y(4)*gg+&
                        fac2*y(4)**2/bb*gg*&
                        (1.0/qq/rr1*(gg3+gg2)-&
                        jj2/qq/rr*(gg1*rpx+gg0*rpt)+&
                        gg1*dum1)-&
                        fac1*fac2*dum*jj2/qq*gg*(gg1*bpx+gg0*bpt)
              dydx(4) = -fac2*y(6)*y(5)/rr1/qq*bpt
              dydx(5) = 0.0
              dydx(6) = 0.0
                 return
          end subroutine derivs1

!
!
!
      subroutine contourpoints!(psi)
        implicit none
        integer :: i, j, kk, nz0, nr0, m(1), n 
        real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
        real(sp), allocatable :: tmpr(:), tmpz(:), tmppsi(:),tmpcscoef(:,:), tmpl(:)
        real(sp), allocatable :: vr(:), vz(:), cscoef1dr(:,:,:),cscoef1dz(:,:,:)
        real(sp) :: dpsidr, dpsidz, xr, xz

        integer :: nrmagx, nzmagx
        real(sp) :: dtheta, tt, tt1, tt2, tt3, tt4
        real(sp) :: dl, rr, zz, ll, npsi1
!

        tt1 = atan((zmax-zmagx)/(rmax-rmagx))
        tt2 = atan((zmax-zmagx)/(rmin-rmagx)) + pi
        tt3 = atan((zmin-zmagx)/(rmin-rmagx)) + pi
        tt4 = atan((zmin-zmagx)/(rmax-rmagx)) + two_pi

        !print *, tt1, tt2, tt3, tt4
        !pause 1234

        nrmagx = inrlocate(gridr, rmagx)
        nzmagx = inrlocate(gridz, zmagx)
        !print *, nrmagx,nzmagx
        nr0 = nrmax - nrmin + 1
        nz0 = nzmax - nzmin + 1

        allocate (vr(nr0), vz(nz0))
        allocate (cscoef1dr(4,nr0-1,nz0), cscoef1dz(4,nz0-1,nr0))

        vpsi = vpsi / bmag / rmag / rmag  ! normlise

        do i = 1, nr0
                vr(i) = gridr(nrmin+i-1)
        end do  
        do i = 1, nz0
                vz(i) = gridz(nzmin+i-1)
        end do

        dtheta = two_pi/(ntheta-1)
        
       ! print *, rmagx, zmagx
       ! normlise grid (r,z)
        gridr = gridr / rmag
        gridz = gridz / rmag
        n = npsi*2
        allocate(tmppsi(n),tmpcscoef(4,n-1),tmpl(n))

        do j = 1, ntheta
                lines(1,j,1) = rmagx
                lines(2,j,1) = zmagx
        end do
                
        do j = 1, ntheta
                tt = (j-1)*dtheta

                if ( (tt >= 0 .and. tt< tt1) .or. (tt >= tt4 ) .or. (tt>= tt2 .and. tt< tt3) ) then 
                        
                        if (tt >= tt2 .and. tt< tt3) then
                                xr = rmin
                        else                            
                                xr = rmax
                        end if

                        xz = zmagx + tan(tt)*(xr-rmagx)
                        dr = (xr - rmagx)/real(n-1,sp)
                        dz = (xz - zmagx)/real(n-1,sp)
                        dl =sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
        
                        tmpl(1) = 0.
                        tmppsi(1) = 0. 
                

                        do i = 2, n
                                tmpl(i) = (i-1)*dl
                                rr = rmagx + (i-1)*dr
                                zz = zmagx + (i-1)*dz
                               ! use units same to sub check_eq 
                                rr = rr/rmag
                                zz = zz/rmag
                                tmppsi(i) =cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
                                if (abs(tmppsi(i)) < abs(tmppsi(i-1))) then
                                         !print *, 'warning for tt 1',
                                         !i,j
                                        do kk= i, n
                                                tmppsi(kk) = tmppsi(i-1)+ (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
                                        end do
                                        exit
                                end if
                        end do
                        !pause
                        call inrcsnak(tmppsi,tmpl,tmpcscoef)
!                        
                        do i = 2, npsi
                                ll =inrcsval(tmppsi,tmpcscoef,vpsi(i),0)
                                rr = rmagx + ll*cos(tt)
                                zz = zmagx + ll*sin(tt)
                                lines(1,j,i) = rr
                                lines(2,j,i) = zz
                        end do          
                end if
!
                if ( (tt >= tt1 .and. tt< tt2) .or. (tt >= tt3 .and. tt<tt4)  ) then

                        if  (tt >= tt1 .and. tt< tt2) then
                                xz = zmax
                        else
                                xz = zmin
                        end if

                        xr = cotan(tt)*(xz-zmagx) + rmagx

                        dr = (xr - rmagx)/real(n-1,sp)
                        dz = (xz - zmagx)/real(n-1,sp)
                        dl =sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
        
                        tmpl(1) = 0.
                        tmppsi(1) = 0.          

                        do i = 2, n
                                tmpl(i) = (i-1)*dl
                                rr = rmagx + (i-1)*dr
                                zz = zmagx + (i-1)*dz
                                !use units in sub check_eq
                                rr = rr/rmag
                                zz = zz/rmag
                                tmppsi(i) =cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
                                if (abs(tmppsi(i)) < abs(tmppsi(i-1))) then
                                         !print *, 'warning for tt 1',
                                         !i,j
                                        do kk= i, n
                                                tmppsi(kk) = tmppsi(i-1)+ (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
                                        end do
                                        exit
                                end if
                        end do
                        !pause
                        call inrcsnak(tmppsi,tmpl,tmpcscoef)
                       
                        do i = 2, npsi
                                ll =inrcsval(tmppsi,tmpcscoef,vpsi(i),0)
                                rr = rmagx + ll*cos(tt)
                                zz = zmagx + ll*sin(tt)
                                lines(1,j,i) = rr
                                lines(2,j,i) = zz
                                
                        end do          
                end if

        end do

                

        open(111,file='contourpsi.dat',status='replace', action='write')
                ! write fluxa
                write(111,101) npsi-ii, ntheta
                do i = 1, npsi-ii     ! revised by lmyu
                        write(111,103) vpsi(i)
                end do

                pw = vpsi(npsi-ii)

                do i = 1, npsi-ii
                        do j = 1, ntheta
                                write(111,104) lines(1:2,j,i) 
                        end do
                end do
         npsi1 = npsi - ii
         Rmax1 = maxval(lines(1,:,npsi1))
         Rmin1 = minval(lines(1,:,npsi1))
         Zmax1 = maxval(lines(2,:,npsi1)) 
         Zmin1 = minval(lines(2,:,npsi1)) 
         !print *, 'Rmax1,Rmin1=', Rmax1, Rmin1, Zmax1, Zmin1
        close(111)
        100 FORMAT (I5)
        101 FORMAT (2I5)
        103 FORMAT (ES25.16)
        104 FORMAT (2ES25.16)
        105 FORMAT (3ES25.16)
        106 FORMAT (4ES25.16)
        !pause 11
         return 
      end subroutine contourpoints
!
      subroutine contourpointsr!(psi)
        implicit none
        integer :: i, j, kk, nz0, nr0, m(1), n 
        real(sp), allocatable :: tmpr(:), tmpz(:), tmppsi(:), tmpcscoefr(:,:),tmpcscoef(:,:), tmpl(:)
        real(sp), allocatable :: vr(:), vz(:), cscoef1dr(:,:,:),cscoef1dz(:,:,:),tmppsir(:)
        real(sp) :: dpsidr, dpsidz, xr, xz, rpsi

        integer :: nrmagx, nzmagx
        real(sp) :: dtheta, tt, tt1, tt2, tt3, tt4
        real(sp) :: dl, rr, zz, ll
        real(sp) :: dum, dum1
!

        tt1 = atan((zmax-zmagx)/(rmax-rmagx)) !按矩形算，gfile是矩形区域
        tt2 = atan((zmax-zmagx)/(rmin-rmagx)) + PI_D
        tt3 = atan((zmin-zmagx)/(rmin-rmagx)) + PI_D
        tt4 = atan((zmin-zmagx)/(rmax-rmagx)) + TWOPI_D 

        !print *, tt1, tt2, tt3, tt4
        !pause 1234

        nrmagx = inrlocate(gridr, rmagx) !定位磁轴
        nzmagx = inrlocate(gridz, zmagx)

        nr0 = nrmax - nrmin + 1 !r,z格点数目
        nz0 = nzmax - nzmin + 1
        n = npsi*2 !用于延拓

        allocate (vr(nr0), vz(nz0))
        allocate (cscoef1dr(4,nr0-1,nz0), cscoef1dz(4,nz0-1,nr0))
        allocate(vrpsi(npsi),vpsir(npsi),Rgeom(npsi),r_rho(npsi))
        allocate(tmppsi(n),tmpcscoef(4,n-1),tmpcscoefr(4,npsi-1),tmpl(n),tmppsir(npsi))
        allocate(cscoefrpsi(4,npsi-1))

        do i = 1, nr0
                vr(i) = gridr(nrmin+i-1)
        end do
        do i = 1, nz0
                vz(i) = gridz(nzmin+i-1)
        end do

        dtheta = TWOPI_D/(ntheta-1)

        !print *, rmagx, zmagx
        vpsir(1) = 0
        tmppsir(1) = 0
        do j = 1, ntheta
                lines(1,j,1) = rmagx !psi0
                lines(2,j,1) = zmagx
        end do

        !计算r
        tt = PI_D
        xr = rmin
        xz = zmagx
        dr = (xr - rmagx)/real(n-1,sp)
        dz = 0.
        dl = sqrt((xr-rmagx)**2)/real(n-1,sp) !

        tmpl(1) = 0.
        tmppsi(1) = 0.

        do i = 2, n
                tmpl(i) = (i-1)*dl
                rr = rmagx + (i-1)*dr
                zz = zmagx + (i-1)*dz
                rr = rr/rmag
                zz = zz/rmag
                tmppsi(i) = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
                !print *, 'tmppsi=', rr,zz,tmpl(i),tmppsi(i)
                if (abs(tmppsi(i)) < abs(tmppsi(i-1))) then
                             !print *, 'warning for tt 1', i,j
                        do kk= i, n
                                tmppsi(kk) = tmppsi(i-1) +(kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
                        end do
                        exit
                end if
        end do
       !pause
        call inrcsnak(tmppsi,tmpl,tmpcscoef)

        do i = 2, npsi
                ll = inrcsval(tmppsi,tmpcscoef,vpsi(i),0)!利用tmppsi插值(vpsi(i)，theta)对应的ll

                tmppsir(i) = ll
        end do
        tt = 0
        xr = rmax
        xz = zmagx
        dr = (xr - rmagx)/real(n-1,sp)
        dz = 0.
        dl = sqrt((xr-rmagx)**2)/real(n-1,sp) !

        tmpl(1) = 0.
        tmppsi(1) = 0.

        do i = 2, n
                tmpl(i) = (i-1)*dl
                rr = rmagx + (i-1)*dr
                zz = zmagx + (i-1)*dz
                rr = rr/rmag
                zz = zz/rmag
                tmppsi(i) = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
                if (abs(tmppsi(i)) < abs(tmppsi(i-1))) then!边界以外延拓？
                               !print *, 'warning for tt 1', i,j
                        do kk= i, n
                                tmppsi(kk) = tmppsi(i-1) +(kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
                        end do
                        exit
                end if
        end do
       !pause
        call inrcsnak(tmppsi,tmpl,tmpcscoef)
        do i = 2, npsi
                ll =inrcsval(tmppsi,tmpcscoef,vpsi(i),0)!利用tmppsi插值(vpsi(i))对应的ll
                tmppsir(i) = (ll+tmppsir(i))*0.5
        end do
        ! coefs of r for given psi by lmyu
        call inrcsnak(vpsi,tmppsir,cscoefrpsi)
        !求均匀r对应的R,Z
        call inrcsnak(tmppsir,vpsi,tmpcscoefr)
        vrpsi(1) = 0.
        do i = 2 , npsi
                vpsir(i) = (i-1)*tmppsir(npsi)/real(npsi-1,sp)
                vrpsi(i) = inrcsval(tmppsir,tmpcscoefr,vpsir(i),0)
        end do
        do j = 1, ntheta
                tt = (j-1)*dtheta

                if ( (tt >= 0 .and. tt< tt1) .or. (tt >= tt4 ) .or. (tt>= tt2 .and. tt< tt3) ) then

                        if (tt >= tt2 .and. tt< tt3) then
                                xr = rmin
                        else
                                xr = rmax
                        end if

                        xz = zmagx + tan(tt)*(xr-rmagx)
                        dr = (xr - rmagx)/real(n-1,sp)
                        dz = (xz - zmagx)/real(n-1,sp)
                        dl =sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)

                        tmpl(1) = 0.
                        tmppsi(1) = 0.


                        do i = 2, n
                                tmpl(i) = (i-1)*dl
                                rr = rmagx + (i-1)*dr
                                zz = zmagx + (i-1)*dz
                                rr = rr/rmag
                                zz = zz/rmag
                                tmppsi(i) =cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
                                if (abs(tmppsi(i)) < abs(tmppsi(i-1))) then
                                         !print *, 'warning for tt 1',
                                         !i,j
                                        do kk= i, n
                                                tmppsi(kk) = tmppsi(i-1)+ (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
                                        end do
                                        exit
                                end if
                        end do
                        !pause
                        call inrcsnak(tmppsi,tmpl,tmpcscoef)
                        if(j == 1)Rgeom(1) = 0.
                        do i = 2, npsi
                                ll =inrcsval(tmppsi,tmpcscoef,vrpsi(i),0) !利用tmppsi插值(vpsi(i)，theta)对应的ll
                                rr = rmagx + ll*cos(tt)!求(vpsi(i)，theta)对应的(R,Z)
                                zz = zmagx + ll*sin(tt)
                                lines(1,j,i) = rr
                                lines(2,j,i) = zz
                                smallr(j,i) = ll
       if(j == 1) Rgeom(i) = rr
       if(j == (ntheta-1)/2) Rgeom(i) = (Rgeom(i) + rr)/2
                        end do
                       !pause
                end if
                if ( (tt >= tt1 .and. tt< tt2) .or. (tt >= tt3 .and. tt<tt4)  ) then

                        if  (tt >= tt1 .and. tt< tt2) then
                                xz = zmax
                        else
                                xz = zmin
                        end if

                        xr = cotan(tt)*(xz-zmagx) + rmagx

                        dr = (xr - rmagx)/real(n-1,sp)
                        dz = (xz - zmagx)/real(n-1,sp)
                        dl =sqrt((xr-rmagx)**2+(xz-zmagx)**2)/real(n-1,sp)
                        tmpl(1) = 0.
                        tmppsi(1) = 0.

                        do i = 2, n
                                tmpl(i) = (i-1)*dl
                                rr = rmagx + (i-1)*dr
                                zz = zmagx + (i-1)*dz
                                rr = rr/rmag
                                zz = zz/rmag
                                tmppsi(i) =cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,rr,zz,0,0)
                                if (abs(tmppsi(i)) < abs(tmppsi(i-1)))then
                                        !print *, 'warning for tt 1',
                                        !i,j
                                        do kk= i, n
                                                tmppsi(kk) = tmppsi(i-1)+ (kk-i+1)*(tmppsi(i-1)-tmppsi(i-2))
                                        end do
                                        exit
                                end if
                        end do
                       !pause
                        call inrcsnak(tmppsi,tmpl,tmpcscoef)

                        do i = 2, npsi
                                ll =inrcsval(tmppsi,tmpcscoef,vrpsi(i),0) !利用tmppsi插值(vpsi(i)，theta)对应的ll
                                rr = rmagx + ll*cos(tt)!求(vpsi(i)，theta)对应的(R,Z)
                                zz = zmagx + ll*sin(tt)
                                lines(1,j,i) = rr
                                lines(2,j,i) = zz
                                smallr(j,i) = ll
                        end do
                end if
        end do

        do i = 1 , npsi
                r_rho(i) = vpsir(i)*sqrt(Rgeom(i)/Rgeom(npsi))
        end do
     deallocate(tmpcscoef)
        allocate (tmpcscoef(4,npsi-2))
        call inrcsnak(r_rho(2:npsi),Vpsir(2:npsi),tmpcscoef)
       rhom_r = inrcsval(r_rho(2:npsi),tmpcscoef,0.5*r_rho(npsi),0)
     deallocate(tmpcscoef)

        open(111,file='contourrpsi.dat',status='replace',action='write')
        open(112,file='qprofile.dat',status='replace',action='write')
                ! write flux
                write(111,101) npsi, ntheta
                do i = 1, npsi
                        write(111,103) vpsir(i)
                end do
                ! minor radius: a by lmyu
                rminor = inrcsval(vpsi,cscoefrpsi,vpsi(npsi-ii),0)
                print *, 'rminor=',rminor
                do i = 1, npsi-ii
                        do j = 1, ntheta
                                write(111,104) lines(1:2,j,i)
                        end do
                        dum = inrcsval(vpsi,cscoefrpsi,vpsi(i),0)/rminor
                        dum1 = inrcsval(gridpsi,cscoefqpsi,vpsi(i),0)
                        write(112,104) dum, dum1
                end do
        close(111)
        close(112)

        100 FORMAT (I5)
        101 FORMAT (2I5)
        103 FORMAT (ES25.16)
        104 FORMAT (2ES25.16)
        105 FORMAT (3ES25.16)
        106 FORMAT (4ES25.16)
        return
       end subroutine contourpointsr

!
       subroutine check_eq ! computing BR,BZ et.al from gfile 
        implicit none
        !common/aspect/e!aspect rate:e=a/R_0
        real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
        integer,parameter :: n=201,na=361
        real(r8) :: R(na),Z(na)
        integer :: i,j
        !real(r8) :: e,dt,dx,Phi,ddum
        real(r8) :: dt,dx,Phi,ddum
        real(r8) :: x(n),t(na),Phi1(na),Phi2(n)



!
!         gfile='g055851.003703'
         !gfile='gfile' ! filename
         gfile='g055851.003703'
!          gfile = 'g063887.002100.q.circ'
!          gfile = 'g048605.04686 '
         ! gfile = 'g098840.012264'
        inquire(file=trim(gfile),exist=alive)
        if(.not.alive) then
                print *, "gfile:" // trim(gfile) // " does not exist."
                stop
        end if

        open(10, file=trim(gfile),status="old",action="read")
                read(10,"(a100)") str !The last two elements should be read into nxefit, nyefit !读取第一行
                nr=len_trim(str)
                nl=index(str(1:nr),' ',.true.)
                nr2=len_trim(str(1:nl))
                nl2=index(str(1:nr2),' ',.true.)
                read(str(nl2:nr),*) nxefit,nyefit !65,65
                allocate(fpol(nxefit),pres(nxefit),workk1(nxefit),workk2(nxefit),fold(nxefit,nyefit),qpsi(nxefit))

                read(10,"(5e16.9)") xdim,zdim,rcentr,rgrid1,zmid !line 2
                read(10,"(5e16.9)") rmagx,zmagx,simagx,sibdry,bcentr!line 3
                read(10,"(5e16.9)") cpasma,Simagx,Xdum,Rmagx,Xdum!line 4 
                read(10,"(5e16.9)") zmagx,xdum,sibdry,xdum,xdum!line 5 
                !xdim,zdim      ; Size of the domain in meters;Horizontal dimension in meter of computational box, Vertical dimension in meter of computational box
                !rcentr,bcentr  ; Reference vacuum toroidal field (m,T); R in meter of vacuum toroidal magnetic fieldBCENTR, Vacuum toroidal magnetic field in Tesla at rcentr
                !rgrid1         ; R of left side of domain; Minimum R in meter of rectangular computational box
                !zmid           ; Z at the middle of the domain; Z of center of computational box in meter
                !rmagx,zmagx    ; Location of magnetic axis
                !simagx         ; Poloidal flux at the axis (Weber /rad);
                !sibdry         ; Poloidal flux at plasma boundary (Weber / rad)
                !cpasma  ; Plasma current in Ampere
                if (sibdry - simagx > 0) then
                        signj = 1.0
                else
                        signj = -1.0
                end if
                read(10,"(5e16.9)") fpol(1:nxefit) !line 6 to 31 Poloidal current function on uniform flux grid
                                                   !fpol=-g(psi) Poloidal current function in m-T, F = RB_phi on flux grid
                !if (fpol(1) < 0.) fpol = - fpol !使q>0 !!!必要吗？ pfz
                read(10,"(5e16.9)") pres(1:nxefit) !line 32 to 57 Plasma pressure in nt/m^2 on uniform flux grid P
                read(10,"(5e16.9)") workk1(1:nxefit)!line 58 to 83 FF'(psi) in (mT)^2 / (Weber /rad) on uniform flux grid
                read(10,"(5e16.9)") workk2(1:nxefit)!line 84 to 109 P'(psi) in (nt /m^2) / (Weber /rad) on uniform flux grid
                read(10,"(5e16.9)") fold(1:nxefit,1:nyefit)!line 110 to 3438 Poloidal flux in Weber/rad on grid points  Poloidal
                                                           !flux in Weber / rad on the rectangular grid points; psi
                read(10,"(5e16.9)") qpsi(1:nxefit) !line 3439 to 3464 q values on uniform flux grid ; q values on uniform
                !pause

                read(10,*) nbdry,nlim !number of plasma boundary points and wall boundary/limiter points
                if(nbdry>0) then
                        allocate(rzbdry(2,nbdry))
                        read(10,"(5e16.9)") rzbdry(1:2,1:nbdry) !line 3466 to 3502 (R,Z) of Plasma boundary points
                else
                        allocate(rzbdry(2,1))
                        rzbdry = 0
                end if
                if(nlim>0) then
                        allocate(rzlim(2,nlim))
                        read(10,"(5e16.9)") rzlim(1:2,1:nlim) !line 3503 to 3523 (R,Z) of Wall boundary points
                else
                        allocate(rzlim(2,1))
                        rzlim = 0
                end if
        close(10)
!      
      nxefit1=int(nxefit/3)
      !write(*,*) 'nxefit1=',nxefit1
      allocate(qpsi1(nxefit1))
      allocate(vpsi(npsi))
      allocate (lines(4,ntheta,npsi))
      allocate(smallr(ntheta,npsi))
      allocate(gridr(nxefit),gridz(nyefit),gridpsi(nxefit),gridpsi1(nxefit1),knotr(nxefit+kr),knotz(nyefit+kz))
      allocate(cscoeffpol(4,nxefit-1),cscoefqpsi(4,nxefit-1),cscoefq(4,nxefit1-1))
      allocate(bscoefpsi(nxefit,nyefit),bscoefbphi(nxefit,nyefit),bscoefscalb(nxefit,nyefit),bscoefscal2(nxefit,nyefit))
      allocate(vec_b_r(nxefit,nyefit),vec_b_z(nxefit,nyefit),vec_b_p(nxefit,nyefit))
      allocate(vec_gradb_r(nxefit,nyefit),vec_gradb_z(nxefit,nyefit),vec_gradb_p(nxefit,nyefit))
      allocate(vec_gradbr_r(nxefit,nyefit),vec_gradbr_z(nxefit,nyefit),vec_gradbr_p(nxefit,nyefit))
      allocate(vec_gradbz_r(nxefit,nyefit),vec_gradbz_z(nxefit,nyefit),vec_gradbz_p(nxefit,nyefit))
      allocate(vec_gradrbp_r(nxefit,nyefit),vec_gradrbp_z(nxefit,nyefit),vec_gradrbp_p(nxefit,nyefit))
      allocate(scal_b(nxefit,nyefit),scal2(nxefit,nyefit))
      !generate the grid points on R and Z direction, the rectangular
      !grid points
       do i=1,nxefit
                gridr(i) = rgrid1 + xdim*(i-1)/(nxefit-1)
       end do
        do j=1,nyefit
                !gridz(j) = (zmid-0.5*zdim) + zdim*(j-1)/(nyefit-1)!(zmid-0.5*zdim)=zbelow
                gridz(j) = (zmid-0.5*zdim) + zdim*(j-1)/(nyefit-1)-zmagx!shift flux for zmagx =0 by lmyu
        end do
        ! shift boundary for z_magx=0 by lmyu
        do i=1,nbdry
           rzbdry(2,i) = rzbdry(2,i) - zmagx
        end do
!!!!!!
        Rmax = maxval(rzbdry(1,:)) !最外侧点
        Rmin = minval(rzbdry(1,:)) !最内侧点
        Zmax = maxval(rzbdry(2,:)) !最高点
        Zmin = minval(rzbdry(2,:)) !最低点
        !print *, Rmax, Rmin, Zmax, Zmin
        !pause
        zmagx = 0.0 ! set z_magx=0 by lmyu
        nrmax = inrlocate(gridr,rmax) + 1 !定位最大最小值点
        nrmin = inrlocate(gridr,rmin)
        nzmax = inrlocate(gridz,zmax) + 1
        nzmin = inrlocate(gridz,zmin)


        do i=1,npsi
          vpsi(i) = (sibdry-simagx)*(i-1)/real(npsi-1,sp)  ! revised by lye
        end do

!!!!!!
        ! set psi = 0 at magnetic axis, by lye
        do i = 1, nxefit
                do j = 1, nyefit
                        fold(i,j) = fold(i,j) - simagx
                        !if (fold(i,j) < 0.) print *, i, j , fold(i,j)
                end do
        end do
!               print *, simagx!, fold(:,65)
        
        !pause
        !generate the uniform flux grid points
        do i=1,nxefit
                !gridpsi(i) = simagx+(sibdry-simagx)*(i-1)/(nxefit-1)
                gridpsi(i) = (sibdry-simagx)*(i-1)/real(nxefit-1,sp)  !revised by lye
        end do
        do i=1,nxefit1
                !gridpsi(i) = simagx+(sibdry-simagx)*(i-1)/(nxefit-1)
                gridpsi1(i) = (sibdry-simagx)*(i-1)/real(nxefit1-1,sp)
        end do
         j=1
        do i=2,nxefit,3
           qpsi1(j)=qpsi(i)
           write(33,*) qpsi(i)
           j=j+1
        enddo
        
        ! write for matlab plot
        ! write RZ grids
        open(111,file='RZ_plt.dat',status='replace', action='write')
                write(111,101) nxefit
                write(111,101) nyefit
                write(111,103) rmagx
                write(111,103) zmagx
                write(111,103) gridr(1)
                write(111,103) gridr(nxefit) 
                write(111,103) gridz(1)
                write(111,103) gridz(nyefit)
                do j = 1, nyefit
                        do i = 1, nxefit
                                write(111,103) gridr(i)
                        end do
                end do
                do j = 1, nyefit
                        do i = 1, nxefit
                                write(111,103) gridz(j)
                        end do
                end do
        close(111)

        100 FORMAT (I5)
        101 FORMAT (2I5)
        103 FORMAT (ES25.16)
        104 FORMAT (2ES25.16)
        105 FORMAT (3ES25.16)
        106 FORMAT (4ES25.16)

        open(111,file='flux2d.dat',status='replace', action='write')
                ! write flux
                write(111,101) nxefit, nyefit
                do j = 1, nyefit
                        do i = 1, nxefit
                                write(111,103) fold(i,j) 
                        end do
                end do
                ! write boundary
                write(111, "(4i5)") nbdry,nlim
                do i=1,nbdry
                        write(111,104) rzbdry(1:2,i)
                end do
                do i=1,nlim
                        write(111,104) rzlim(1:2,i)
                end do
                ! Write magnetic axis
        close(111)
! contourpoints

                     
!normalise
        rmag = rmagx       ! major radius of magnetic axis (m)
        !rmaj = rmag*1e2    ! major axis in cm- given in equilibrium file
        rmaj = 1000.0 ! n=6 TAE,ITPA, A. Konies et al., Nucl. Fusion 58, 126027 (2018).
        !rmaj = 165.0 ! hl-2a
        !rmaj = 300.0 ! n=3 TAE, Z.Y. Liu GMEC
                     ! https://w3.pppl.gov/~ngorelen/TAE_lin.html
        !rmaj = 173.8 ! d3d
        !bmag = abs(fpol(1))/rmag !magnetic field at axis (T)
        bmag = 3.0 ! n=6 TAE, A. Könies et al., Nucl. Fusion 58, 126027 (2018)
        !bmag = 1.33 ! hl-2a
        !bmag = 1.0 ! n=3 TAE Z.Y. Liu GMEC
                   ! https://w3.pppl.gov/~ngorelen/TAE_lin.html
        !bmag = 1.958 ! d3d
        bkg = bmag*1e1 !  b at the magnetic axis in kgauss
        gridr = gridr /rmag
        gridz = gridz / rmag
        gridpsi = gridpsi / bmag / rmag / rmag
        gridpsi1 = gridpsi1 / bmag / rmag / rmag
        fpol = abs(fpol) / bmag /rmag ! set fpol positive. negative of fpol just shows in the sign of the moving direction of particles.
        fold = fold / bmag /rmag /rmag
         write(8,*) 'bmag=',bmag
      !generate the knot sequence on R and Z
      !direction;对于某个网格点Vpsi(i)，寻找出其等高线上的点(R1,
      !Z1)……(Rn, Zn)
        call cdbbsnak(gridr,kr,knotr) !生成4阶r节点
        call cdbbsnak(gridz,kz,knotz) !生成4阶z节点

        !calculate coef of fpol on uniform flux grid points
        call inrcsnak(gridpsi,fpol,cscoeffpol)!求g的psi各点一维三次样条插值系数 !!!
        ! calculate coef of qprofile on uniform flux grid points
        call inrcsnak(gridpsi1,qpsi1,cscoefq)
        call inrcsnak(gridpsi,qpsi,cscoefqpsi)
        !calculate coef of psi on (R,Z)
        call cdbbscoef2d(gridr,gridz,fold,knotr,knotz,kr,kz,bscoefpsi) !求出psi B样条基的系数
      do i=1,nxefit
                do j=1,nyefit
                        vec_b_r(i,j) = -cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),0,1)/gridr(i) ! BR=-PpsiPZ/R
                        vec_b_z(i,j) =  cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),1,0)/gridr(i) !BZ=PpsiPR/R
                        tmpfpol =inrcsval(gridpsi,cscoeffpol,fold(i,j),0) !在psi网格上插值g gives value of g in (R,Z) grid
                        vec_b_p(i,j) = tmpfpol/gridr(i)   !Bphi=g/R
                        scal_b(i,j) = sqrt(vec_b_r(i,j)**2+ vec_b_z(i,j)**2+vec_b_p(i,j)**2) !scal_b stands for scale magnetic field on (R,Z)
                        scal2(i,j) = tmpfpol
                end do
        end do
!
        allocate(bscoefvec1r(nxefit,nyefit),bscoefvec1z(nxefit,nyefit),bscoefvec1p(nxefit,nyefit))
!
        call cdbbscoef2d(gridr,gridz,vec_b_p,knotr,knotz,kr,kz,bscoefvec1p)
        call cdbbscoef2d(gridr,gridz,vec_b_r,knotr,knotz,kr,kz,bscoefvec1r)
        call cdbbscoef2d(gridr,gridz,vec_b_z,knotr,knotz,kr,kz,bscoefvec1z)
        call cdbbscoef2d(gridr,gridz,scal_b,knotr,knotz,kr,kz,bscoefscalb)
        call cdbbscoef2d(gridr,gridz,scal2,knotr,knotz,kr,kz,bscoefscal2)
!
        do i=1,nxefit
                do j=1,nyefit
                        vec_gradb_r(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,gridr(i),gridz(j),1,0)
                        vec_gradb_z(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,gridr(i),gridz(j),0,1)
                        vec_gradb_p(i,j) = 0.0
                        vec_gradbr_r(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1r,gridr(i),gridz(j),1,0)
                        vec_gradbr_z(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1r,gridr(i),gridz(j),0,1)
                        vec_gradbr_p(i,j) = 0.0
                        vec_gradbz_r(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1z,gridr(i),gridz(j),1,0)
                        vec_gradbz_z(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefvec1z,gridr(i),gridz(j),0,1)
                        vec_gradbz_p(i,j) = 0.0
                        vec_gradrbp_r(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefscal2,gridr(i),gridz(j),1,0)
                        vec_gradrbp_z(i,j) = cdbbsval2d(knotr,knotz,kr,kz,bscoefscal2,gridr(i),gridz(j),0,1)
                        vec_gradrbp_p(i,j) = 0.0
                end do
        end do
!
        allocate(bscoefvec2r(nxefit,nyefit),bscoefvec2z(nxefit,nyefit),bscoefvec2p(nxefit,nyefit))
        allocate(bscoefvec3r(nxefit,nyefit),bscoefvec3z(nxefit,nyefit),bscoefvec3p(nxefit,nyefit))
        allocate(bscoefvec4r(nxefit,nyefit),bscoefvec4z(nxefit,nyefit),bscoefvec4p(nxefit,nyefit))
        allocate(bscoefvec5r(nxefit,nyefit),bscoefvec5z(nxefit,nyefit),bscoefvec5p(nxefit,nyefit))


        !coef of vec_2
        call cdbbscoef2d(gridr,gridz,vec_gradb_r,knotr,knotz,kr,kz,bscoefvec2r)
        call cdbbscoef2d(gridr,gridz,vec_gradb_z,knotr,knotz,kr,kz,bscoefvec2z)
        !coef of vec_3
        call cdbbscoef2d(gridr,gridz,vec_gradbr_r,knotr,knotz,kr,kz,bscoefvec3r)
        call cdbbscoef2d(gridr,gridz,vec_gradbr_z,knotr,knotz,kr,kz,bscoefvec3z)
        !coef of vec_4
        call cdbbscoef2d(gridr,gridz,vec_gradbz_r,knotr,knotz,kr,kz,bscoefvec4r)
        call cdbbscoef2d(gridr,gridz,vec_gradbz_z,knotr,knotz,kr,kz,bscoefvec4z)
        !coef of vec_5
        call cdbbscoef2d(gridr,gridz,vec_gradrbp_r,knotr,knotz,kr,kz,bscoefvec5r)
        call cdbbscoef2d(gridr,gridz,vec_gradrbp_z,knotr,knotz,kr,kz,bscoefvec5z)

          
!
      deallocate(fpol)
      deallocate(pres)
      deallocate(workk1)
      deallocate(workk2)
      deallocate(fold)
      deallocate(qpsi)
      if (allocated(rzbdry)) then
        deallocate(rzbdry)
      endif
      if (allocated(rzlim)) then
        deallocate(rzlim)
      endif
      deallocate(qpsi1)
      deallocate(vpsi)
      deallocate(lines)
      deallocate(smallr)
      deallocate(gridr)
      deallocate(gridz)
      deallocate(gridpsi)
      deallocate(gridpsi1)
      deallocate(knotr)
      deallocate(knotz)
      deallocate(cscoeffpol)
      deallocate(cscoefqpsi)
      deallocate(cscoefq)
      deallocate(bscoefpsi)
      deallocate(bscoefbphi)
      deallocate(bscoefscalb)
      deallocate(bscoefscal2)
      deallocate(vec_b_r)
      deallocate(vec_b_z)
      deallocate(vec_b_p)
      deallocate(vec_gradb_r)
      deallocate(vec_gradb_z)
      deallocate(vec_gradb_p)
      deallocate(vec_gradbr_r)
      deallocate(vec_gradbr_z)
      deallocate(vec_gradbr_p)
      deallocate(vec_gradbz_r)
      deallocate(vec_gradbz_z)
      deallocate(vec_gradbz_p)
      deallocate(vec_gradrbp_r)
      deallocate(vec_gradrbp_z)
      deallocate(vec_gradrbp_p)
      deallocate(scal_b)
      deallocate(scal2)
      deallocate(bscoefvec1r)
      deallocate(bscoefvec1z)
      deallocate(bscoefvec1p)
      deallocate(bscoefvec2r)
      deallocate(bscoefvec2z)
      deallocate(bscoefvec2p)
      deallocate(bscoefvec3r)
      deallocate(bscoefvec3z)
      deallocate(bscoefvec3p)
      deallocate(bscoefvec4r)
      deallocate(bscoefvec4z)
      deallocate(bscoefvec4p)
      deallocate(bscoefvec5r)
      deallocate(bscoefvec5z)
      deallocate(bscoefvec5p)
      end subroutine check_eq
!
        function expE(E)
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          real(r8):: E,expE
          real(r8) :: Ec32,E32,tE
          !Ec32=(Ec)**1.5
          !E32=E**1.5+Ec32
          !tE=(E-E0)/(Ed)
          !tE=-1.0*tE*tE
          !expE=2.0*exp(tE)/sqrt(pi)/Ed/E32 ! slow down
          expE=exp(-E) ! maxwell,isotropic
        endfunction expE
!
        function erfcE(E)
          implicit none
          real(r8):: E,erfcE
          real(r8) :: Ec32,E32,tE
          Ec32=(Ec)**1.5
          E32=E**1.5+Ec32
          tE=(E-E0)/(Ed)
          erfcE=erfc(tE)/E32
        endfunction erfcE
!
        function experfcE(E)
          implicit none
          real,parameter ::pi=3.14159265358979,two_pi=6.28318530717959d0
          real(r8):: E,experfcE
          real(r8) :: tE
          tE=(E-E0)/Ed
          experfcE=2.0*exp(-1.0*tE*tE)/&
                   erfc(tE)/sqrt(pi)/Ed
        endfunction experfcE

!
        function sqrtE(E)
          implicit none
          real(r8):: E,sqrtE
          real(r8) :: Ec32,E32
          Ec32=Ec**1.5
          E32=E**1.5+Ec32
          sqrtE=3.0*sqrt(E)/2.0/E32
        endfunction sqrtE

!
        function expL(Lamb)
          implicit none
          real(r8):: Lamb,expL
          real(r8):: tL
          tL = (Lamb - L0)/Ld
          tL = -1.0*tL*tL
          expL= exp(tL)
        endfunction expL

!
        function expL1(E,Lamb)
          implicit none
          real(r8):: Lamb,E,expL1
          real(r8):: tL
          tL = Lamb - L0
          if (E==0.0) E=1e-6
          expL1=-2.0*Lamb*tL/Ld/Ld/E
        endfunction expL1
!
        function expLE(E,Lamb)
          implicit none
          real(r8):: Lamb,E,expLE
          expLE=expL1(E,Lamb)*erfcE(E)
        endfunction expLE

!
       function exppsi(psi)
          implicit none
          real(r8):: psi,exppsi
          real(r8) :: tpsi
!         normalised psi: pol/pw
!         rr0 = pol0/pw, rrd= Delta_psi/pw
          tpsi=(psi-rr0)/rrd
          tpsi=-1.0*tpsi
          exppsi=exp(tpsi)
        endfunction exppsi
!
       function erfcE32(E)
          implicit none
          real(r8) E
          real(r8):: erfcE32
          erfcE32=sqrtE(E)*erfcE(E)
        endfunction erfcE32
!
       function omega_star(Ek,Lamb)
        implicit none
        !common/aspect/e!aspect rate:e=a/R_0
        !real(r8) Ek, psi,e,Lamb
        real(r8) Ek, psi,Lamb
        real(r8) :: omega_star
!       psi = pol/pw 
        omega_star = 1.0/2.0*rhoh_bar*&
                     sqrt(2.0)/e/e/rrd/&
                     (experfcE(Ek)+sqrtE(Ek)+&
                      expL1(Ek,Lamb))
      endfunction omega_star
!
      function pfpe(Ek,Lamb,psi)
       implicit none
       real,parameter ::pi=3.14159265358979,two_pi=6.28318530717959d0
       !common/coef_dist/Cn
       !real(r8) Ek, Lamb,psi,Cn
       real(r8) Ek, Lamb,psi
       real(r8) :: pfpe
!      psi = pol/pw
!      Ek = E/ekev
       !pfpe = - 8.0/sqrt(pi)/Cn*&
       !        expL(Lamb)*exppsi(psi)*&
       !        (expE(Ek)+erfcE32(Ek)+&
       !        expLE(Ek,Lamb))
       ! pfpe = -1.0/pi**1.5*exppsi(psi)*&
       !        expE(Ek) ! maxwell, isotropic
       pfpe = -1.0/pi**1.5*den_h(psi)*&
                expE(Ek) ! maxwell, isotropic, itpa
       !pfpe = - 1.0/Cn*(expE(Ek)+erfcE32(Ek))*&
       !               exppsi(psi) ! slow down, isotropic
      endfunction pfpe
!
!
      function pfpz(Ek,Lamb,psi) ! partial F / partial Pz
       implicit none
       real,parameter ::pi=3.14159265358979,two_pi=6.28318530717959d0
       !common/coef_dist/Cn
       !real(r8) Ek, Lamb,psi,Cn
       real(r8) :: pfpz
       real(r8) Ek, Lamb,psi
!      psi = pol/pw
!      Ek = E/ekev
       !pfpz = - 8.0/sqrt(pi)/Cn/rrd/pw*&
       !        expL(Lamb)*exppsi(psi)*erfcE(Ek)
       !pfpz = -1.0/rrd/pw/pi**1.5*&
       !       exppsi(psi)*expE(Ek) ! maxwell,isotropic
       pfpz = 1.0/pw/pi**1.5*&
              denp_h(psi)*expE(Ek) ! maxwell,isotropic itpa
       !pfpz = - 1.0/Cn/rrd/pw*&
       !           exppsi(psi)*erfcE(Ek) ! slow down, isotropic
      endfunction pfpz
!
       function pfpz1(Ek,Lamb,psi)
        implicit none
        !common/aspect/e!aspect rate:e=a/R_0
        !real(r8) Ek, psi,e,Lamb
        real(r8) Ek, psi,Lamb
        real(r8) :: pfpz1
!       psi = pol/pw
!      Ek = E/ekev

        pfpz1 = 1.0/2.0*rhoh_bar*&
                sqrt(2.0)/e/e*&
                pfpz(Ek,Lamb,psi)
      endfunction pfpz1
!

!      subroutine setgrid
!      implicit none
!      integer n,i,nog
!      parameter(n=800)
!      real(r8) :: xg(-3:n+3),delta,h,xend,xini,q0,q1,cq1
!      common/xgrida/xg
!      common/xgrid1/delta
!      common/nogr/nog
!      common/qfundat/q0,q1,cq1,xini,xend
!
      !write(0,*)"input nog"
      !read(5,*) nog
      
!  
!      h= (xend-xini)/float(nog)
!       h=(x2-x1)/float(nog)
!      delta=h

!      do i=-3,nog+3
!      xg(i)=h*float(i)+xini
!      enddo
!      return
!      end subroutine setgrid
!
!       function den(x)
!       implicit none
!       integer i,nfit
!       parameter(nfit=5)
!       common/dendat/cden(nfit)
!       common/denf/dalpha,dedge,dedge1
!       real(r8) :: x,xx,den,cden,dalpha,x1,dedge,dedge1
!       den=cden(1)
!       x1=x
!       if(x.gt.1.0D0)x1=1.0D0
!       xx=x1**2
!       do i=1,nfit-1
!       den=den+cden(i+1)*xx**i
!       enddo
!       den=den-dedge*exp((x1-1.0)/dalpha)
!       den=1.0D0

!       end function den
!
!       function denp(x)
!       implicit none
!       integer i,nfit
!       parameter(nfit=5)
!       common/dendat/cden(nfit)
!       common/denf/dalpha,dedge,dedge1
!       real(r8) :: x,xx,denp,cden,dalpha,dedge,dedge1
!       denp=0.0
!       xx=x**2
!       do i=1,nfit-1
!       denp=denp+(2*i)*cden(i+1)*x**(2*i-1)
!       enddo
!       denp=denp-dedge*exp((x-1.0)/dalpha)/dalpha
!       if(x.gt.1.0D0)denp=0.0D0
!       denp=0.0D0

!       end function denp

!
!      function b_nor(t,x) ! B/B(0)
!       implicit none
!       common/aspect/e
!       real(r8) :: e, x, t
!       real(r8) ::b_nor
       
!       b_nor = sqrt(1.0+(e*jr2(x)/qfun(x))**2*grr(t,x))/rf(t,x)
      ! b_nor = sqrt(1.0+(e*jr2(x)/qfun(x))**2/rf(t,x)/rf(t,x))/rf(t,x)

!      end function b_nor
         
!
!
!      function b_norpt(t,x)
!       implicit none
!       common/aspect/e
!       real(r8) :: e, x, t
!       real(r8) ::b_norpt
!       real(r8) :: rr,rr1,gr,gr1
!       real(r8) :: dum,dum1,qq
!       rr=rf(t,x)
!       rr1=rfpt(t,x)
!       gr=grr(t,x)
!       gr1=grrpt(t,x)
!       dum=e*jr2(x)
!       qq=qfun(x)
!       dum1=sqrt(1.0+dum**2/qq**2*gr)
!       b_norpt = -1.0/rr/rr*rr1*dum1+&
!                 0.5/rr/dum1*dum**2/qq**2*gr1
       !b_norpt = kappa_t(t,x)*b_nor(t,x)
!      end function b_norpt

!
!
!      function b_norpx(t,x)
!       implicit none
!       common/aspect/e
!       real(r8) :: e, x, t
!       real(r8) ::b_norpx
!       real(r8) :: rr,rr1,rr2,dum,dum1,qq2
!       real(r8) :: dum2, dum3
!       real(r8) :: gr,gr1
!       rr=rf(t,x)
!       rr1=rfpr(t,x)
!       rr2=rf(t,x)**2
!       dum=(e*jr2(x))**2
!       qq2=qfun(x)**2
!       gr = grr(t,x)
!       gr1 = grrpx(t,x)
!       dum1=sqrt(1.0+dum/qq2*gr)
!       dum2 = 2.0*jr2(x)/qfun(x)*&
!              (jr2p(x)/qfun(x)-jr2(x)/qq2*dqfun(x))
!       dum3 = jr2(x)**2/qq2
!       b_norpx = -1.0/rr2*rr1*dum1+&
!                  0.5/rr/dum1*e**2*(gr*dum2+gr1*dum3)
       !b_norpx = kappa_r(t,x)*b_nor(t,x)
!      end function b_norpx

!

!
!      function b_star(t,x,v4) ! b**(theta,r,v||)
!       implicit none
!       common/aspect/e
!       real(r8) :: t,x,v4,e,b2,g
!       real(r8) :: g1,g2,g3,g4
!       real(r8) :: j1,j2,j3,j4,j5,j6
!       real(r8) :: rr,rr1,rr2,qq,qq1
!       real(r8) :: b_star,bb
!       bb=b_nor(t,x)
!       b2=b_nor(t,x)**2
!       g=gzz(t,x)
!       g1=grrpx(t,x)
!       g2=grr(t,x)
!       g3=gtrpt(t,x)
!       g4=grt(t,x)
!       j1=jf(t,x)
!       j2=jfpx(t,x)
!       j3=jfpt(t,x)
!       j4=jr2(x)
!       j5=jr2p(x)
!       rr=rf(t,x)
!       rr1=rfpr(t,x)
!       rr2=rfpt(t,x)
!       qq=qfun(x)
!       qq1=dqfun(x)
!       b_star = bb + v4*g*j4/b2/qq*&
!                (g1+j2/j1*g2+g3+j3/j1*g4)-&
!                v4/b2/qq/rr*j4**g*(g2*rr1+g4*rr2)+&
!                v4/b2*(j5/qq-j4*qq1/qq**2)*g*g2
       !b_star = b_nor(t,x)
!      end function b_star
!

!
!      function fem(i,m,x)
!      implicit none
!      integer n,i,m,nog
!      parameter(n=800)
!      common/nogr/nog
!      real(r8) :: fem,x
!        if(i.eq.1) then
!          fem=felst(m,x)
!        else if(i.eq.(nog-1))then
!          fem=felend(m,x)
!        else
!          fem=fel(i,x)
!        endif
!       end function fem
!
!      function femp(i,m,x)
!      implicit none
!      integer n,i,m,nog
!      parameter(n=800)
!      common/nogr/nog
!      real(r8) :: x,femp
!      if(i.eq.1)then
!      femp=felstp(m,x)
!      else if(i.eq.(nog-1))then
!      femp=felendp(m,x)
!      else
!      femp=felp(i,x)
!      endif
!      return
!      end function femp

!
!      function fel(i,x)
!      implicit none
!      integer n
!      parameter(n=800)
!      real(r8) :: xg(-3:n+3)
!      common/xgrida/xg
!      integer i,i0,im1,im2,ip1,ip2
!      real(r8) :: fel,x,xbar
!      i0=i
!      im1=i-1
!      im2=i-2
!      ip1=i+1
!      ip2=i+2
      
!      if(x.gt.xg(im2).and.x.le.xg(im1))then
!      xbar=(x-xg(im2))/(xg(im1)-xg(im2))
!      fel=xbar**3
!      else if(x.gt.xg(im1).and.x.le.xg(i0))then
!  
!      xbar=(x-xg(i0))/(xg(i0)-xg(im1))
!      fel=4.0D0-xbar**2*(6.0D0+3.0D0*xbar)
!  
!      else if(x.gt.xg(i0).and.x.le.xg(ip1))then
!      xbar=(x-xg(i0))/(xg(ip1)-xg(i0))
!      fel=4.0D0-xbar**2*(6.0D0-3.0D0*xbar)

!      else if(x.gt.xg(ip1).and.x.le.xg(ip2))then
!      xbar=(x-xg(ip2))/(xg(ip2)-xg(ip1))
!      fel=-xbar**3
!      else
!      fel=0.0D0
!      endif
!
!      return
!      end function fel
!
!      function felend(m,x)
!  
!      implicit none
!      integer n,m,nog
!      parameter(n=800)
!      real(r8) :: xg(-3:n+3)
!      common/xgrida/xg
!      real(r8) :: felend,x,delta,cm1,cn,deltax,augden
!      common/xgrid1/delta
!      common/nogr/nog

!      augden=1.0D0
!      deltax=delta*(denp(augden)/den(augden)+1.0D0)
!      cm1=(3.0D0+deltax)/(3.0D0-deltax)
!      cn = -0.5D0*deltax/(3.0D0-deltax)
!      felend=cm1*fel(nog-1,x)+cn*fel(nog,x)-fel(nog+1,x)
!!      return
!      end function felend
!
!      function felst(m,x)

!      implicit none
!      integer m
!      real(r8) :: felst,x
!
!      if(m.ne.1)then
!      felst=fel(1,x)+fel(-1,x)-0.5D0*fel(0,x)
!      else
!      felst=fel(1,x)-fel(-1,x)
!      endif
!      return
!      end function felst
!
!      function felstp(m,x)
!      implicit none
!      integer m
!      real(r8) :: felstp,x
!      if(m.ne.1)then
!      felstp=felp(1,x)+felp(-1,x)-0.5D0*felp(0,x)
!      else
!      felstp=felp(1,x)-felp(-1,x)
!      endif
!      return
!      end function felstp
!
!     function felendp(m,x)
!      implicit none
!      integer n,m,nog
!      parameter(n=800)
!      real(r8) :: xg(-3:n+3)
!      common/xgrida/xg
!      real(r8) :: felendp,x,delta,cm1,cn,deltax,augden
!      common/xgrid1/delta
!      common/nogr/nog
!      augden=1.0D0
!      deltax=delta*(denp(augden)/den(augden)+1.0D0)
!      cm1=(3.0D0+deltax)/(3.0D0-deltax)
!      cn = -0.5D0*deltax/(3.0D0-deltax)
!      felendp=cm1*felp(nog-1,x)+cn*felp(nog,x)-felp(nog+1,x)
!      return
!      end function felendp
!
!      function felp(i,x)
!      implicit none
!      integer n
!      parameter(n=800)
!      real(r8) :: xg(-3:n+3)
!      common/xgrida/xg
!      integer i,i0,im1,im2,ip1,ip2
!      real(r8) :: felp,x,xbar
!      i0=i
!      im1=i-1
!      im2=i-2
!      ip1=i+1
!      ip2=i+2
!
!      if(x.gt.xg(im2).and.x.le.xg(im1))then
!      xbar=(x-xg(im2))/(xg(im1)-xg(im2))
!      felp=3.0D0*xbar**2/(xg(im1)-xg(im2))
!  
!      else if(x.gt.xg(im1).and.x.le.xg(i0))then
!  
!      xbar=(x-xg(i0))/(xg(i0)-xg(im1))
!      felp=-(12.0D0*xbar+9.0D0*xbar**2)/(xg(i0)-xg(im1))
!  
!      else if(x.gt.xg(i0).and.x.le.xg(ip1))then
!      xbar=(x-xg(i0))/(xg(ip1)-xg(i0))
!      felp=-(12.0D0*xbar-9.0D0*xbar**2)/(xg(ip1)-xg(i0))

!      else if(x.gt.xg(ip1).and.x.le.xg(ip2))then
!      xbar=(x-xg(ip2))/(xg(ip2)-xg(ip1))
!      felp=-3.0D0*xbar**2/(xg(ip2)-xg(ip1))
!      else
!      felp=0.0D0
!      endif
!  
!      return 
!      end function felp
!
!     subroutine xrange(j,i,xa,xb)
!      implicit none
!      integer n,nq,i,j,nx,ny,nog
!      parameter(n=800,nq=4)
!      common/nogr/nog
!      real(r8) :: xa,xb,delta,xg(-3:n+3)
!      real(r8) :: eps
!      real(r8) :: q0,q1,cq1,xini,xend
!      common/xgrida/xg
!      common/xgrid1/delta
!      common/qfundat/q0,q1,cq1,xini,xend


!      eps=1.0D-6
!      eps=xini+eps
     ! if(j.eq.1)then ! fixed 20230304 yulm
!      if(j.eq.0)then  
!      xa=eps
!      xb=xg(j+2)
!      else if(j.eq.(nog-1))then
!      xa=xg(j-2)
!      xb=xg(nog)
!      else
!      xa=xg(j-2)
!      xb=xg(j+2)
!      if(i.lt.j)xb=xg(i+2) ! fixed  20230304 yulm
!      if(i.gt.j)xa=xg(i-2)
!      if(xa.lt.eps)xa=eps
!      endif

!      if(i.lt.j)xb=xg(i+2)
!      if(i.gt.j)xa=xg(i-2)
!      if(xa.lt.eps)xa=eps

      !ny=(xb-xa)/delta+0.01D0
      !nx=ny*nq
!      return
!      end subroutine xrange



!     mesh grid
      function fun_alpha(x,y,lam,sigma)
        implicit none
        !common/aspect/e!
        !real(r8) :: e, fun_alpha, x, y, lam
        real(r8) ::  fun_alpha, x, y, lam
        real(r8) :: rr,zz, sigma
        !rr = 1.0 + e*x
        !zz = 0.0
        if (y >=0) then
           zz = 0.0
        else
           zz = PI_D
        endif

        !b_nor = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,rr,zz,0,0)
        if (sigma==1.0) then
           fun_alpha = asin(sqrt(abs(1.0-lam*b_nor(zz,x)))) ! fixed negative value with alpha=0
        else
           fun_alpha = -asin(sqrt(abs(1.0-lam*b_nor(zz,x))))
        end if
      end function
!
       function fun_alphap(x,y,lam,sigma)
        implicit none
        !common/aspect/e
        !real(r8) :: fun_alphap,x,y, lam, e
        real(r8) :: fun_alphap,x,y,lam
        real(r8) :: rr,zz, sigma,bnor1,bnorpr1,bnorpt1
        real(r8) :: dum
        !rr = 1.0 + e*x
        !zz = 0.0
        if (y >=0) then
           zz = 0.0
        else
           zz = PI_D
        endif
        !b_nor = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,rr,zz,0,0)
        !b_nor_pr=cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,rr,zz,1,0)
        bnor1 = b_nor(zz,x)
        !bnorpr1 = b_norpx(zz,x)/rfpr(zz,x)
        bnorpr1 = b_norpx(zz,x)/rfpx(zz,x)
        bnorpt1 = b_norpt_bar_rfpt(zz,x) ! fixed rfpt(0,x)=0
        dum = bnorpr1 + bnorpt1
        if (sigma==1.0) then
           fun_alphap = 1.0/sqrt(lam*bnor1)*0.5/sqrt(abs(1.0-lam*bnor1))*&
                    (-1.0)*lam*e*dum ! fixed negative value with alpha=0
        else
           fun_alphap = 1.0/sqrt(lam*bnor1)*0.5/sqrt(abs(1.0-lam*bnor1))*&
                    lam*e*dum

        end if
      end function
!
!     arc length derive: ds/dx
      function fun_s(x,y,lam,sigma)
        implicit none
        real(r8) :: fun_s, x,y,lam,sigma
!
        fun_s = sqrt(1.0 + fun_alphap(x,y,lam,sigma)**2)
      end function
!
      subroutine check_mesh
        implicit none
        ! grids set up in (X, alpha) 
        ! nx: X at alpha = 0, na: alpha at X=Xmin
        ! num: X for alpha as function of X
        ! num1: X for arc length integral
        integer,parameter :: nx=14,num=101,na=7,num1=21
        !common/aspect/e
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,temp,vbg
        ! s: arc length for X at alpha = 0
        real(r8) ::  s(nx-1),sdum(num)
        ! lamb/lamdum: mube
        real(r8) :: xgrid0(nx+1),xgrid(nx),lamb(nx),&
                    xx(num),xxx(num1),fi(num1)
        ! aa: alpha for X=Xmin, Z>0
        real(r8) :: lamdum(na),aa0(na+1),aa(na)
        real(r8), allocatable :: ss(:)
        real(r8), allocatable :: cscoefsx(:,:)
        real(r8), allocatable :: en0(:)
        integer :: n,n1,n2,n3,i,j,k,m,l,nn,kk
        integer :: nlam
!
730      format('# XR, alpha(pitch angle), mube, kth')
         open(unit = 6,file = 'mesh_xal_arc.dat')
         write(6,730)
731      format('# XR, alpha(pitch angle), mube, sigma')
         open(unit = 7,file = 'mesh_xal.dat')
         write(7,731)

441      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,I6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
!
        allocate(alpha(nx,num),alphadum(na,num))
!        allocate(en0(nE))
!
        adum = (Rmax1 -Rmin1)/2.0
        xmin = (Rmin1 -rmagx)/adum
        xmax = (Rmax1 - rmagx)/adum
!       xtp > xmin
        n=size(xgrid0)
        n1=size(xx)
        n2=size(xxx)
        call linspace(xmin,xmax,n,xgrid0)
        call linspace(Ea,Eb,nE,en0)
        dx = xgrid0(2)-xgrid0(1)
!        dE = en0(2)-en0(1)
        do i=1,n-1
           xgrid(i)=xgrid0(i) + dx/2.0
        enddo
!        do i=1,nE-1
!           en0(i) = en0(i) + dE/2.0
!        enddo
        zz = 0.0
        sdum0 = 0.0
        do i=1,nx ! mube for alpha = 0
           rr = 1.0 + e*xgrid(i)
           bnor = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,rr,zz,0,0)
           lamb(i) = 1.0/bnor
        enddo
        allocate(cscoefsx(4,num-1))
        kk = 0
        do i=1,2 ! sigma
           if (i==1) then
                 sigma = 1.0
           else
                 sigma = -1.0
           endif
         
           do j= 1,nx-1 ! alpha curve for same arc length
              call linspace(xgrid(j),xmax,n1,xx)
              do k=1,num
                 ! alpha curve
                 alpha(j,k)=fun_alpha(xx(k),xx(k),lamb(j),sigma)
                 write(7,442) xx(k), alpha(j,k), &
                              lamb(j), sigma
                 if (k>1) then ! arc length
                    call linspace(xx(1),xx(k),n2,xxx)
                    dx = abs(xxx(2)-xxx(1))
                    do l=2,num1
                       fi(l) = fun_s(xxx(l),xxx(l),lamb(j),sigma)
                    enddo
                    fi(1) = -fi(3)+2.0*fi(2) ! fix Inf for dalpha/dx
                    call simp(n2,dx,fi,sdum(k))
                 else
                     sdum(k) = 0.0
                 endif
              enddo
              !print *,'xx',xx, 'sdum',sdum
!
              call inrcsnak(sdum,xx,cscoefsx)
!
              nn = max(2,floor(sdum(num)/(xmax-xmin)*real(nx))+1)
!             
              allocate(ss(nn))
              call linspace(sdum0,sdum(num),nn,ss)
              do l=1,nn ! X for equal arc length
                 kk = kk + 1
                 xdum = inrcsval(sdum,cscoefsx,ss(l),0)
                 write(6,441) xdum, fun_alpha(xdum,xdum,lamb(j),sigma), &
                              lamb(j), kk
              enddo
              deallocate(ss)
            enddo
        enddo
!
       !else
!      xtp < xmin
       n3=size(aa0)
       vbg = 0.0
       call linspace(vbg,PIO2,n3,aa0)
       call linspace(xmin,xmax,n1,xx)
       da = aa0(2)-aa0(1)
       do i=1,n3-1
       !do i=n3-1,n3-1 ! alpha = 1
          aa(i)=aa0(i)+da/2.0
       enddo
       do i=1,na ! alpha at X=Xmin, alpha \= 0
       !do i=n3-1,n3-1 ! lambda = 0
          rr = 1.0 + e*xmin
          bnor = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,rr,zz,0,0)
          lamdum(i) = (1.0- sin(aa(i))**2)/bnor
       enddo
!        
       do i=1,2 
       !do i=1,1 ! sigma = 1
          if (i==1) then
             sigma = 1.0
          else
             sigma = -1.0
          endif
!
          do j= 1,na
          !do j=n3-1,n3-1 ! lambda = 0
              do k=1,num
                 alphadum(j,k)=fun_alpha(xx(k),xx(k),lamdum(j),sigma)
                 write(7,442) xx(k), alphadum(j,k), &
                              lamdum(j), sigma
                 if (k>1) then
                    call linspace(xx(1),xx(k),n2,xxx)
                    dx = abs(xxx(2)-xxx(1))
                    do l=1,num1
                       fi(l) = fun_s(xxx(l),xxx(l),lamdum(j),sigma)
                    enddo
                    call simp(n2,dx,fi,sdum(k))
                  else
                     sdum(k) = 0.0
                 endif
              enddo
!
              !print *, sdum
              call inrcsnak(sdum,xx,cscoefsx)
!
              nn = max(2,floor(sdum(num)/(xmax-xmin)*real(nx))+1)
              !print *,'nn=',nn
              allocate(ss(nn))
              call linspace(vbg,sdum(num),nn,ss)
              do l=1,nn ! X for equal arc length
                 kk = kk + 1
                 xdum = inrcsval(sdum,cscoefsx,ss(l),0)
                 write(6,441) xdum, fun_alpha(xdum,xdum,lamdum(j),sigma), &
                              lamdum(j), kk
              enddo
              deallocate(ss)
          enddo
       enddo
 
      ! end if

       end subroutine check_mesh
!
      subroutine check_mesh1 ! Lambda = 0 numerical equalibrium
        implicit none
        ! grids set up in (X, alpha) 
        ! nx: X at alpha = 0, na: alpha at X=Xmin
        ! num: X for alpha as function of X at Lambda = 0
        ! num1: X for arc length integral
        integer,parameter :: nx=9,num=101,na=51,num1=21
        !common/aspect/e
        !rr=R/R0, zz=Z/R0=0
        !da: dalpha=pi/2/na
        !dE: grid of energy
        !xmin=Rmin-R0,xmax=Rmax-R0
        !Rmin:high field side major radius
        !Rmax:low filed side major radius 
        !R0: magnetic axis major radius
        !real(r8) :: e,lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,temp,vbg
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        ! s: arc length for X at alpha = 0
        real(r8) ::  s(nx-1),sdum(num)
        ! lamb/lamdum: mube
        ! xx: grids of [Xmin,Xmax], X=(R-R0)/rminor
        ! xxx: grids of the ith X grid 
        real(r8) :: xgrid0(nx+1),xgrid(nx),lamb(nx),&
                    xx(num),xxx(num1),fi(num1)
        ! aa: grids of alpha for X=Xmin, Z=0
        real(r8) :: lamdum(na),aa0(na+1),aa(na)
        ! ss: arc length
        ! cscoefx: spline interp coefs between X and s
        real(r8), allocatable :: ss(:),xdum0(:)
        real(r8), allocatable :: cscoefsx(:,:)
        real(r8), allocatable :: en0(:)
        integer :: n,n1,n2,n3,i,j,k,m,l,nn,kk,iE
        integer :: nlam
!
730      format('# XR, alpha(pitch angle), mube, E(keV),dx(a),dE, kth')
         open(unit = 6,file = 'mesh_xal_arc.dat')
         write(6,730)
731      format('# XR, alpha(pitch angle), mube, E(keV) sigma')
         open(unit = 7,file = 'mesh_xal.dat')
         write(7,731)
732      format('# Ei(keV)')
         open(unit = 8, file = 'grid_E.dat')
         write(8,732)

441      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,&
               1F12.6,1x,1F12.6,1x,I6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
!
        allocate(alpha(nx,num),alphadum(na,num))
        allocate(en0(nE))

!
        adum = (Rmax1 -Rmin1)/2.0
        xmin = (Rmin1 -rmagx)/rminor
        xmax = (Rmax1 - rmagx)/rminor
        
        n=size(xgrid0)
        n1=size(xx)
        n2=size(xxx)
        zz = 0.0
        sdum0 = 0.0

        allocate(cscoefsx(4,num-1))
        kk = 0
!      xtp < xmin
        n3=size(aa0)
        vbg = 0.0
       call linspace(vbg,PIO2,n3,aa0)
       call linspace(xmin,xmax,n1,xx)
       da = aa0(2)-aa0(1)
       aa(1)=aa0(n3-1)+da/2.0 ! v||/v ~1.0
       rr = 1.0 + e*xmin
       bnor = cdbbsval2d(knotr,knotz,kr,kz,bscoefscalb,rr,zz,0,0)
       lamdum(1) = (1.0- sin(aa(1))**2)/bnor
!      
       call linspace(Ea,Eb,nE,en0)
       dE = en0(2)-en0(1)
       do i=1,nE-1
           en0(i) = en0(i) + dE/2.0
           write(8,443) en0(i), dE
       enddo
!        
      do i=1,1
          if (i==1) then
             sigma = 1.0
          else
             sigma = -1.0
          endif
!
        do iE = 1, nE-1 ! energy grids should be nE-1?
          do j= 1,1     ! pitch/mube grids
              do k=1,num
                 ! Given lambda, obtaining alphas at X grids
                 alphadum(j,k)=fun_alpha(xx(k),xx(k),lamdum(j),sigma)
                 write(7,442) xx(k), alphadum(j,k), &
                              lamdum(j), en0(iE), sigma
                 if (k>1) then
                    ! xxx of grids [Xmin, Xk]
                    call linspace(xx(1),xx(k),n2,xxx)
                    dx = abs(xxx(2)-xxx(1))
                    do l=1,num1
                       fi(l) = fun_s(xxx(l),xxx(l),lamdum(j),sigma)
                    enddo
                    call simp(n2,dx,fi,sdum(k))
                  else
                     sdum(k) = 0.0
                 endif
              enddo
!
             ! arc length s corresponding to X
              call inrcsnak(sdum,xx,cscoefsx)
             !real(nx) should be real(num)?
              nn = max(2,floor(sdum(num)/(xmax-xmin)*real(nx))+1)
              !
              allocate(ss(nn),xdum0(nn))
              call linspace(vbg,sdum(num),nn,ss)
              ! 
            ds = ss(2) - ss(1)
            do l = 1, nn
               xdum0(l) = inrcsval(sdum,cscoefsx,ss(l),0) ! X for equal arc len, used to calculate deltax x(xdum_i)-x(xdum_i-1)
            enddo
              do l=1,nn-1 ! X for equal arc length should be nn-1
                 kk = kk + 1 ! ion index
                 tmps = ss(l) + ds/2.0
                 ! xdum: Xi at si+ds/2
                 xdum = inrcsval(sdum,cscoefsx,tmps,0) 
                 ! deltax=x(Xi+1)-x(Xi)
                 rr1 = 1.0 + e*xdum0(l)
                 rr2 = 1.0 + e*xdum0(l+1)
                 pol1 = cdbbsval2d(knotr,knotz,kr,kz,&
                       bscoefpsi,rr1,zz,0,0)   ! flux psi of psi(R,Z)
                 pol2 = cdbbsval2d(knotr,knotz,kr,kz,&
                       bscoefpsi,rr2,zz,0,0)   ! flux psi of psi(R,Z)

                 x1 = inrcsval(vpsi,cscoefrpsi,pol1,0)
                 x2 = inrcsval(vpsi,cscoefrpsi,pol2,0)
                 dxdum = abs(x2-x1)/rminor
                 write(6,441) xdum, fun_alpha(xdum,xdum,lamdum(j),sigma), &
                              lamdum(j), en0(iE), dxdum, dE, kk
              enddo
              deallocate(ss)
              deallocate(xdum0)
          enddo
        enddo
      enddo


       end subroutine check_mesh1
!
      subroutine check_mesh2 ! Lambda = 0 analytic equalibrium
        implicit none
        ! grids set up in (X, alpha) 
        ! nx: X at alpha = 0, na: alpha at X=Xmin
        ! num: X for alpha as function of X at Lambda = 0
        ! num1: X for arc length integral
        integer,parameter :: num=101,na=51,num1=21
        !common/aspect/e
        !rr=R/R0, zz=Z/R0=0
        !da: dalpha=pi/2/na
        !dE: grid of energy
        !xmin=Rmin-R0,xmax=Rmax-R0
        !Rmin:high field side major radius
        !Rmax:low filed side major radius 
        !R0: magnetic axis major radius
        !real(r8) :: e,lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,vbg,rdum
        real(r8) :: tmpfpol,temp,temp1,temp2
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        real(r8) :: dl,  x0, ddx, thetadum, phidum
        real(r8) :: thetadum1,thetadum2, dpzdum1,dpzdum2
        real(r8) :: dum, dum1, dum2, dum3, dum4
        ! s: arc length for X at alpha = 0
        real(r8) ::  s(nrho-1),sdum(num)
        ! lamb/lamdum: mube
        ! xx: grids of [Xmin,Xmax], X=(R-R0)/rminor
        ! xxx: grids of the ith X grid 
        real(r8) :: xgrid0(nrho+1),xgrid(nrho),lamb(nrho),&
                    xx(num),xx1(num),xx2(num),xxx(num1),fi(num1)
        ! aa: grids of alpha for X=Xmin, Z=0
        real(r8) :: lamdum(na),aa0(na+1),aa(na)
        real(r8) :: tini,tfinal,a0,a1
        ! ss: arc length
        ! cscoefx: spline interp coefs between X and s
        real(r8), allocatable :: ss(:),xdum0(:),rx(:)
        real(r8), allocatable :: cscoefsx(:,:),cscoefrx(:,:),&
                                 cscoefrx1(:,:)
        real(r8), allocatable :: en0(:)
        integer :: n,n1,n2,n3,i,j,k,m,l,nn,kk,iE
        integer :: nlam,istep,ik
!
730      format('# x, theta, phi, XR, alpha(pitch angle), &
               mube, E(keV),dx(a),dE, kth')
         open(unit = 6,file = 'mesh_xal_arc.dat')
         write(6,730)
731      format('# XR, alpha(pitch angle), mube, E(keV) sigma')
         open(unit = 7,file = 'mesh_xal.dat')
         write(7,731)
732      format('# Ei(keV)')
         open(unit = 8, file = 'grid_E.dat')
         write(8,732)
733      format('# nprt a0 a1 tini tfinal')

441      format(1F12.8,1x,1F12.8,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,&
                1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,I6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
444      format(I6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
!
        allocate(alpha(nrho,num),alphadum(na,num),rx(num))
        allocate(en0(nE))
        rr = 0.0
        rr1 = PI_D 
        zz = 1.0 ! LCFS flux
!      
        !adum = (Rmax1 -Rmin1)/2.0
        adum = (rf(rr,zz) - rf(rr1,zz))/2.0
        !xmin = (Rmin1 -rmagx)/rminor
        xmin = (rf(rr1,zz)-1.0)/e + 0.3
        !print *, 'xmin=',xmin
        !xmax = (Rmax1 - rmagx)/rminor
        xmax = (rf(rr,zz) - 1.0)/e - 0.3
        !print *, 'xmax=',xmax
        
        n=size(xgrid0)
        n1=size(xx) ! n1=num
        n2=size(xxx)
        !zz = 0.0
        sdum0 = 0.0

        allocate(cscoefsx(4,num-1),cscoefrx(4,num-1),cscoefrx1(4,num-1))
        kk = 0
!      xtp < xmin
        n3=size(aa0)
        vbg = 0.0
       call linspace(vbg,PIO2,n3,aa0)
       call linspace(xmin,vbg,n1,xx)
       call linspace(vbg,xmax,n1,xx1)
       dl = 1e-6
       do j=1,n1
          dum1 = 1.0 + e*xx(j)
          if (xx(j) < 0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr1,dum1)
            rx(j) = x0
            write(1003,*) xx(j),x0
          else
            x0 = 0.0
            rx(j) = x0
            write(1003,*) xx(j),x0
          endif
     
       enddo

       call inrcsnak(xx,rx,cscoefrx) ! X<=0, given X, obtain x=r/a
!     
       do j=1,n1
          dum1 = 1.0 + e*xx1(j)
          if (xx1(j) > 0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr,dum1)
            rx(j) = x0
            write(1003,*) xx1(j),x0
          else
            x0 = 0.0
            rx(j) = x0
            write(1003,*) xx1(j),x0
          endif
       enddo
       
       call inrcsnak(xx1,rx,cscoefrx1) ! X>=0, given X, obtain x=r/a
!
       da = aa0(2)-aa0(1)
       aa(1)=aa0(n3-1)+da/2.0 ! v||/v ~1.0
       bnor = b_nor(rr1,zz)
       lamdum(1) = (1.0- sin(aa(1))**2)/bnor ! mube
!      
       call linspace(Ea,Eb,nE,en0)
       dE = en0(2)-en0(1)
       do i=1,nE-1
           en0(i) = en0(i) + dE/2.0
           write(8,443) en0(i), dE
       enddo
!     
!     (E,mube,X) grid set up
      do i=1,1
          if (i==1) then
             sigma = 1.0
          else
             sigma = -1.0
          endif
!
          call linspace(xmin,xmax,n1,xx2)
        do iE = 1, nE-1 ! energy grids should be nE-1?
          do j= 1,1     ! pitch/mube grids
              do k=1,num ! X grids with equal dX
                 ! Given lambda, obtaining alphas at X grids
                 if (xx2(k) >=0) then
                      rdum = inrcsval(xx1,cscoefrx1,xx2(k),0)
                 else
                      rdum = inrcsval(xx,cscoefrx,xx2(k),0)
                 endif
                 alphadum(j,k)=fun_alpha(rdum,xx2(k),lamdum(j),sigma)
                 write(7,442) xx2(k), alphadum(j,k), &
                              lamdum(j), en0(iE), sigma
                 if (k>1) then
                    ! xxx of grids [Xmin, Xk]
                    call linspace(xx2(1),xx2(k),n2,xxx)
                    dx = abs(xxx(2)-xxx(1))
                    do l=1,num1
                       if (xxx(l) >=0) then
                          rdum = inrcsval(xx1,cscoefrx1,xxx(l),0)
                       else
                          rdum = inrcsval(xx,cscoefrx,xxx(l),0)
                       endif

                       fi(l) = fun_s(rdum,xxx(l),lamdum(j),sigma)
                    enddo
                    call simp(n2,dx,fi,sdum(k))  ! arc length of
                                                 ! range [Xmin,xx2(k)]
                  else
                     sdum(k) = 0.0
                 endif
              enddo
!
              ! arc length s corresponding to X
              call inrcsnak(sdum,xx2,cscoefsx) ! given s, obtain X
!             real(nx) should be real(num)?
              nn = max(2,floor(sdum(num)/(xmax-xmin)*real(nrho))+1)
              !
              allocate(ss(nn),xdum0(nn))
              call linspace(vbg,sdum(num),nn,ss)
              ! 
            ds = ss(2) - ss(1)
            do l = 1, nn
               xdum0(l) = inrcsval(sdum,cscoefsx,ss(l),0) ! X for equal arc len, used to calculate deltax x(xdum_i)-x(xdum_i-1)
            enddo
              do l=1,nn-1 ! X grids with equal arc length, should be nn-1
                 kk = kk + 1 ! ion index
                 tmps = ss(l) + ds/2.0
                 ! xdum: Xi at si+ds/2
                 xdum = inrcsval(sdum,cscoefsx,tmps,0) 
                 ! rdum: r/a
                 if (xdum >=0) then
                   rdum = inrcsval(xx1,cscoefrx1,xdum,0)
                   thetadum = 0.0
                 else
                   rdum = inrcsval(xx,cscoefrx,xdum,0)
                   thetadum = PI_D
                 endif
                 ! dx=x(Xi+1)-x(Xi)
                 if (xdum0(l+1)>=0) then
                    x1 = inrcsval(xx1,cscoefrx1,xdum0(l+1),0) 
                    thetadum1 = 0.0
                 else
                    x1 = inrcsval(xx,cscoefrx,xdum0(l+1),0)
                    thetadum1 = PI_D
                 endif
                ! write(1003,*) xdum0(l+1),x1
                 if (xdum0(l)>=0) then
                    x2 = inrcsval(xx1,cscoefrx1,xdum0(l),0)
                    thetadum2 = 0.0
                 else
                    x2 = inrcsval(xx,cscoefrx,xdum0(l),0)
                    thetadum2 = PI_D
                 endif
                 write(1003,*) xdum0(l),x2
                 dxdum = abs(x1-x2)
!
                 phidum = 0.0
!
                 dum = fun_alpha(rdum,xdum,lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 !dum4 = e*(b_norpx(thetadum,rdum)/&
                 !          rfpr(thetadum,rdum)+&
                 !          b_norpt(thetadum,rdum)/&
                 !         rfpt(thetadum,rdum))
                  dum4 = e*(b_norpx(thetadum,rdum)/&
                           rfpx(thetadum,rdum)+&
                           b_norpt(thetadum,rdum)/&
                           rfpt(thetadum,rdum))
                 !dpzdum1 = (e**3*jr2(rdum)/rfpr(thetadum,rdum)-&
                 !         e/rf(thetadum,rdum)**2*&
                 !         (dum2/dum3/&
                 !         b_nor(thetadum,rdum)*lamdum(j)+&
                 !         dum3/b_nor(thetadum,rdum)**2))*&
                 !         (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
                 dpzdum1 = (e**3*jr2(rdum)/rfpx(thetadum,rdum)-&
                          e/rf(thetadum,rdum)**2*&
                          (dum2/dum3/&
                          b_nor(thetadum,rdum)*lamdum(j)+&
                          dum3/b_nor(thetadum,rdum)**2))*&
                          (xdum0(l+1)-xdum0(l))    ! analytic delta Pz 
!
                 !dpzdum1 = (e**3*jr2(rdum)/rfpr(thetadum,rdum)/e+&
                 !         dum4*(dum2/dum3/&
                 !         b_nor(thetadum,rdum)*lamdum(j)+&
                 !         dum3/b_nor(thetadum,rdum)**2))*&
                 !         (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
!
                 dum = fun_alpha(x1,xdum0(l+1),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum1,x1)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x1,0)
!                particle toroidal canonical momentum
                 temp2 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl+1

                 dum = fun_alpha(x2,xdum0(l),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum2,x2)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x2,0)
!                particle toroidal canonical momentum
                 temp1 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl
                 dpzdum2 = temp2 - temp1       ! delta Pz
!
                 print *, 'dpzdum=',dpzdum1, dpzdum2
                 print *, 'E=', en0(iE), 'lambda=',lamdum(j)
!
                 write(6,441) rdum, thetadum, phidum, xdum, &
                              fun_alpha(rdum,xdum,lamdum(j),sigma), &
                              lamdum(j), en0(iE), dxdum, dE, kk
              enddo
              deallocate(ss)
              deallocate(xdum0)
          enddo
        enddo
      enddo
      ik = 1 ! ik: 0 turn off fort.1900 writing 
             ! ik: 1 for turn on fort.1900 writing
      a0 = 0.0
      a1 = 1.0 
      tini = 0.0
      tfinal = 50000
      write(1900,733)
      if (ik.eq.1) then
         write(1900,444) kk,a0,a1,tini,tfinal
       else
          kk=20000
          write(1900,444) kk,a0,a1,tini,tfinal
       endif
       close(1900)
      end subroutine check_mesh2
!
      subroutine check_mesh3 !  analytic equalibrium
        implicit none
        ! grids set up in (X, alpha, E)
        ! nx: X at alpha = 0, na: alpha at X=Xmin
        ! num: X for alpha as function of X at Lambda = 0
        ! num1: X for arc length integral
        integer,parameter :: num=101,na=51,num1=21
        !common/aspect/e
        !rr=R/R0, zz=Z/R0=0
        !da: dalpha=pi/2/na
        !dE: grid of energy
        !xmin=Rmin-R0,xmax=Rmax-R0
        !Rmin:high field side major radius
        !Rmax:low filed side major radius
        !R0: magnetic axis major radius
        !real(r8) :: e,lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,vbg,rdum
        real(r8) :: almin, almax
        real(r8) :: tmpfpol,temp0,temp,temp1,temp2
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        real(r8) :: dl,  x0, ddx, thetadum, phidum
        real(r8) :: thetadum1,thetadum2, dpzdum1,dpzdum2
        real(r8) :: dum, dum1, dum2, dum3, dum4
        ! s: arc length for X at alpha = 0
        real(r8) ::  s(nrho-1),sdum(num)
        ! lamb/lamdum: mube
        ! xx: grids of [Xmin,Xmax], X=(R-R0)/rminor
        ! xxx: grids of the ith X grid
        real(r8) :: xgrid0(nrho+1),xgrid(nrho),lamb(nrho),&
                    xx(num),xx1(num),xx2(num),xxx(num1),fi(num1)
        ! aa: grids of alpha for X=Xmin, Z=0
        !real(r8) :: lamdum(na),aa0(na+1),aa(na),dlam(na)
        real(r8) :: tini,tfinal,a0,a1
        ! ss: arc length
        ! cscoefx: spline interp coefs between X and s
        real(r8), allocatable :: ss(:),xdum0(:),rx(:)
        real(r8), allocatable :: cscoefsx(:,:),cscoefrx(:,:),&
                                 cscoefrx1(:,:)
        real(r8), allocatable :: en0(:),lamdum(:),aa0(:),aa(:),dlam(:)
        integer :: n,n1,n2,n3,i,j,k,m,l,nn,kk,iE
        integer :: nlam,istep,ik
!
730      format('# x, theta, phi, XR, alpha(pitch angle), &
               mube, Pz, E(keV), dL, dP, dE, kth')
         open(unit = 6,file = 'mesh_xal_arc.dat')
         write(6,730)
731      format('# XR, alpha(pitch angle), mube, E(keV) sigma')
         open(unit = 7,file = 'mesh_xal.dat')
         write(7,731)
732      format('# Ei(keV)')
         open(unit = 8, file = 'grid_E.dat')
         write(8,732)
733      format('# nprt a0 a1 tini tfinal')

441      format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                1E16.9,1x,I8)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
444      format(I8,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
448        format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                  1E16.9,1x,1E16.9,1x,I6)

!
        allocate(alpha(nrho,num),alphadum(numL,num),rx(num))
        allocate(en0(nE))
        allocate(lamdum(numL),aa0(numL+1),aa(numL),dlam(numL))
        rr = 0.0
        rr1 = PI_D
        zz = 1.0 ! LCFS flux
!
        !adum = (Rmax1 -Rmin1)/2.0
        adum = (rf(rr,zz) - rf(rr1,zz))/2.0
        !xmin = (Rmin1 -rmagx)/rminor
        xmin = (rf(rr1,zz)-1.0)/e + 1e-9
        !print *, 'xmin=',xmin
        !xmax = (Rmax1 - rmagx)/rminor
        xmax = (rf(rr,zz) - 1.0)/e - 1e-9
        !print *, 'xmax=',xmax

        n=size(xgrid0)
        n1=size(xx) ! n1=num
        n2=size(xxx)! n2=num1
        !zz = 0.0
        sdum0 = 0.0

        allocate(cscoefsx(4,num-1),cscoefrx(4,num-1),cscoefrx1(4,num-1))
        kk = 0
!      xtp < xmin
        bnor = b_nor(rr1,zz)
        almax = asin(sqrt(1.0 - bnor*La)) ! max alpha at Xmin
        almin = asin(sqrt(1.0 - bnor*Lb)) ! min alpha at Xmin
        n3=size(aa0)
        vbg = 0.0
       !call linspace(vbg,PIO2,n3,aa0)
!      aa0: pitch angle,i.e. alpha=sigma_paral*sin^(-1)(abs(v_paral/v))
       call linspace(almin, almax, n3, aa0)
       call linspace(xmin,vbg,n1,xx)
       call linspace(vbg,xmax,n1,xx1)
       dl = 1e-6
       do j=1,n1
          dum1 = 1.0 + e*xx(j)
          if (xx(j) < 0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr1,dum1)
            rx(j) = x0
            write(1003,*) xx(j),x0
          else
            x0 = 0.0
            rx(j) = x0
            write(1003,*) xx(j),x0
          endif

       enddo

       call inrcsnak(xx,rx,cscoefrx) ! X<=0, given X, obtain x=r/a
                                     ! X = (R-R0)/a
!
       do j=1,n1
          dum1 = 1.0 + e*xx1(j)
          if (xx1(j) > 0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr,dum1)
            rx(j) = x0
            write(1003,*) xx1(j),x0
          else
            x0 = 0.0
            rx(j) = x0
            write(1003,*) xx1(j),x0
          endif
       enddo

       call inrcsnak(xx1,rx,cscoefrx1) ! X>=0, given X, obtain x=r/a
!
       da = aa0(2)-aa0(1)
       do i=1,n3-1
          aa(i)=aa0(i)+da/2.0
       enddo
       do i=1,n3-1 ! alpha at X=Xmin, alpha \= 0

          lamdum(i) = (1.0- sin(aa(i))**2)/bnor ! mube corresponds to alpha
          dum = (1.0 - sin(aa0(i))**2)/bnor
          dum1 = (1.0 - sin(aa0(i+1))**2)/bnor         
          dlam(i) = dum - dum1

       enddo
!
       call linspace(Ea,Eb,nE,en0)
       dE = en0(2)-en0(1)
       do i=1,nE-1
           en0(i) = en0(i) + dE/2.0
           write(8,443) en0(i), dE
       enddo
!
!     (E,mube,X) grid set up
      do i=1,1
          if (i==1) then
             sigma = 1.0
          else
             sigma = -1.0
          endif
!
          call linspace(xmin,xmax,n1,xx2) ! X grids
        do iE = 1, nE-1 ! energy grids should be nE-1?
          do j= 1,numL     ! pitch/mube grids
              do k=1,num ! X grids with equal dX
                 ! Given lambda, obtaining alphas at X grids
                 if (xx2(k) >=0) then
                      rdum = inrcsval(xx1,cscoefrx1,xx2(k),0) ! x=r/a
                 else
                      rdum = inrcsval(xx,cscoefrx,xx2(k),0)
                 endif
                 alphadum(j,k)=fun_alpha(rdum,xx2(k),lamdum(j),sigma) ! pitch angle
                 write(7,442) xx2(k), alphadum(j,k), &
                              lamdum(j), en0(iE), sigma
                 if (k>1) then
                    ! xxx of grids [Xmin, Xk]
                    call linspace(xx2(1),xx2(k),n2,xxx)
                    dx = abs(xxx(2)-xxx(1))
                    do l=1,num1
                       if (xxx(l) >=0) then
                          rdum = inrcsval(xx1,cscoefrx1,xxx(l),0)
                       else
                          rdum = inrcsval(xx,cscoefrx,xxx(l),0)
                       endif

                       fi(l) = fun_s(rdum,xxx(l),lamdum(j),sigma)
                    enddo
                    call simp(n2,dx,fi,sdum(k))  ! arc length of
                                                 ! range [Xmin,xx2(k)]
                  else
                     sdum(k) = 0.0
                 endif
              enddo
!
              ! arc length s corresponding to X
              call inrcsnak(sdum,xx2,cscoefsx) ! given s, obtain X
!             real(nx) should be real(num)?
              nn = max(2,floor(sdum(num)/(xmax-xmin)*real(nrho))+1)
              !
              allocate(ss(nn),xdum0(nn))
              call linspace(vbg,sdum(num),nn,ss)
              !
            ds = ss(2) - ss(1)
            do l = 1, nn
               xdum0(l) = inrcsval(sdum,cscoefsx,ss(l),0) ! X for equal arc len, used to calculate deltax x(xdum_i)-x(xdum_i-1)
            enddo
              do l=1,nn-1 ! X grids with equal arc length, should be nn-1
                 kk = kk + 1 ! ion index
                 tmps = ss(l) + ds/2.0
                 ! xdum: Xi at si+ds/2
                 xdum = inrcsval(sdum,cscoefsx,tmps,0)
                 ! rdum: r/a
                 if (xdum >=0) then
                   rdum = inrcsval(xx1,cscoefrx1,xdum,0)
                   thetadum = 0.0
                 else
                   rdum = inrcsval(xx,cscoefrx,xdum,0)
                   thetadum = PI_D
                 endif
                 ! dx=x(Xi+1)-x(Xi)
                 if (xdum0(l+1)>=0) then
                    x1 = inrcsval(xx1,cscoefrx1,xdum0(l+1),0)
                    thetadum1 = 0.0
                 else
                    x1 = inrcsval(xx,cscoefrx,xdum0(l+1),0)
                    thetadum1 = PI_D
                 endif
                ! write(1003,*) xdum0(l+1),x1
                 if (xdum0(l)>=0) then
                    x2 = inrcsval(xx1,cscoefrx1,xdum0(l),0)
                    thetadum2 = 0.0
                 else
                    x2 = inrcsval(xx,cscoefrx,xdum0(l),0)
                    thetadum2 = PI_D
                 endif
                 write(1003,*) xdum0(l),x2
                 dxdum = abs(x1-x2)
!
                 phidum = 0.0
!
                 dum = fun_alpha(rdum,xdum,lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 !dum4 = e*(b_norpx(thetadum,rdum)/&
                 !          rfpr(thetadum,rdum)+&
                 !          b_norpt(thetadum,rdum)/&
                 !          rfpt(thetadum,rdum))
                 !dpzdum1 = (e**3*jr2(rdum)/rfpr(thetadum,rdum)-&
                 !         e/rf(thetadum,rdum)**2*&
                 !         (dum2/dum3/&
                 !         b_nor(thetadum,rdum)*lamdum(j)+&
                 !         dum3/b_nor(thetadum,rdum)**2))*&
                 !         (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
                 dum4 = e*(b_norpx(thetadum,rdum)/&
                           rfpx(thetadum,rdum)+&
                           b_norpt(thetadum,rdum)/&
                           rfpt(thetadum,rdum))
                 dpzdum1 = (e**3*jr2(rdum)/rfpx(thetadum,rdum)-&
                          e/rf(thetadum,rdum)**2*&
                          (dum2/dum3/&
                          b_nor(thetadum,rdum)*lamdum(j)+&
                          dum3/b_nor(thetadum,rdum)**2))*&
                          (xdum0(l+1)-xdum0(l))    ! analytic delta Pz

                 temp = b_nor(thetadum,rdum)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,rdum,0)
                 temp0 = tmpfpol*e**2 - dum3/temp
!
                 !dpzdum1 = (e**3*jr2(rdum)/rfpr(thetadum,rdum)/e+&
                 !         dum4*(dum2/dum3/&
                 !         b_nor(thetadum,rdum)*lamdum(j)+&
                 !         dum3/b_nor(thetadum,rdum)**2))*&
                 !         (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
!
                 dum = fun_alpha(x1,xdum0(l+1),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum1,x1)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x1,0)
!                particle toroidal canonical momentum normalized to
!                R0**2*omega_c
                 temp2 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl+1

                 dum = fun_alpha(x2,xdum0(l),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum2,x2)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x2,0)
!                particle toroidal canonical momentum 
                 temp1 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl
                 dpzdum2 = abs(temp2 - temp1)       ! delta Pz
!
                ! print *, 'dpzdum=',abs(dpzdum1), dpzdum2
                ! print *, 'E=', en0(iE), 'lambda=',lamdum(j)
!
                 write(6,441) rdum, thetadum, phidum, xdum, &
                              fun_alpha(rdum,xdum,lamdum(j),sigma), &
                              lamdum(j), temp0, en0(iE), dlam(j), &
                              dpzdum2, dE, kk
              enddo
              deallocate(ss)
              deallocate(xdum0)
          enddo
        enddo
      enddo
       deallocate(alpha)
       deallocate(alphadum)
       deallocate(rx)
       deallocate(en0)
       deallocate(lamdum)
       deallocate(aa0)
       deallocate(aa)
       deallocate(dlam)

      ik = 1 ! ik: 0 turn off fort.1900 writing
             ! ik: 1 for turn on fort.1900 writing
      a0 = 0.0
      a1 = 1.0
      tini = 0.0
      tfinal = 50000
      write(1900,733)
      if (ik.eq.1) then
         write(1900,444) kk,a0,a1,tini,tfinal
       else
          kk=5
          write(1900,444) kk,a0,a1,tini,tfinal
       endif
       close(1900)
      end subroutine check_mesh3
!
      subroutine check_mesh4 !  analytic equalibrium
        implicit none
        ! grids set up in (X, alpha, E)
        ! nx: X at alpha = 0, na: alpha at X=Xmin
        ! num: X for alpha as function of X 
        ! num1: X for arc length integral
        integer,parameter :: num=101,na=51,num1=21
        !common/aspect/e
        !rr=R/R0, zz=Z/R0=0
        !da: dalpha=pi/2/na
        !dE: grid of energy
        !xmin=(Rmin-R0)/a,xmax=(Rmax-R0)/a
        !Rmin:high field side major radius
        !Rmax:low filed side major radius
        !R0: magnetic axis major radius
        !a: minor radius
        !real(r8) :: e,lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,vbg,vbgp
        real(r8) :: almin, almax,rdum
        real(r8) :: tmpfpol,temp0,temp,temp1,temp2
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        real(r8) :: dl,  x0, ddx, thetadum, phidum
        real(r8) :: thetadum1,thetadum2, dpzdum1,dpzdum2
        real(r8) :: dum, dum1, dum2, dum3, dum4
        ! s: arc length for X at alpha = 0
        real(r8) ::  s(nrho-1),sdum(num)
        ! lamb/lamdum: mube
        ! xx: grids of [Xmin,Xmax], X=(R-R0)/rminor
        ! xxx: grids of the ith X grid
        real(r8) :: xgrid0(nrho+1),xgrid(nrho),lamb(nrho),&
                    xx(num),xx1(num),xx2(num),xxx(num1),fi(num1)
        ! aa: grids of alpha for X=Xmin, Z=0
        !real(r8) :: lamdum(na),aa0(na+1),aa(na),dlam(na)
        real(r8) :: tini,tfinal,a0,a1
        ! ss: arc length
        ! cscoefx: spline interp coefs between X and s
        real(r8), allocatable :: ss(:),xdum0(:),rx(:),xxk0(:),&
                                 xxk(:),xstart(:),xstart0(:),&
                                 temp_lamdum0(:), temp_xstart0(:)
        real(r8), allocatable :: cscoefsx(:,:),cscoefrx(:,:),&
                                 cscoefrx1(:,:),cscoeflamx(:,:)
        real(r8), allocatable :: en0(:),lamdum0(:),lamdum(:),&
                                 aa0(:),aa(:),dlam(:)
        integer :: n,n1,n2,n3,n4,n5,i,j,k,m,l,nn,kk,kk1,iE
        integer :: nlam,istep,ik,numalp
!
730      format('# x, theta, phi, XR, alpha(pitch angle), &
               mube, Pz, E(keV), dL, dP, dE, kth, sigma')
         open(unit = 6,file = 'mesh_xal_arc.dat')
         write(6,730)
731      format('# XR, alpha(pitch angle), mube, E(keV) sigma')
         open(unit = 7,file = 'mesh_xal.dat')
         write(7,731)
732      format('# Ei(keV)')
         open(unit = 8, file = 'grid_E.dat')
         write(8,732)
733      format('# nprt a0 a1 tini tfinal')

441      format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                1E16.9,1x,I8,1x,1F12.6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
444      format(I8,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
448      format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                  1E16.9,1x,1E16.9,1x,I6)

!
        allocate(alpha(nrho,num),alphadum(numL,num),rx(num))
        allocate(en0(nE))
        allocate(lamdum0(numL))
        numalp = numL - nrho ! total mube grid numL
        allocate(aa0(numalp),aa(numalp)) ! alpha_k 
        
        kk = int(nrho/2)  
        kk1 =  nrho - kk         
        allocate(xxk0(kk),xxk(kk1))       ! X_k
        !write(1003,*) 'numalp:',numalp,'kk:',kk,'kk1:',kk1
        rr = 0.0
        rr1 = PI_D
        zz = 1.0 ! LCFS flux
!
        !adum = (Rmax1 -Rmin1)/2.0
        adum = (rf(rr,zz) - rf(rr1,zz))/2.0
        !xmin = (Rmin1 -rmagx)/rminor
        xmin = (rf(rr1,zz)-1.0)/e + 1e-9
        !print *, 'xmin=',xmin
        !xmax = (Rmax1 - rmagx)/rminor
        xmax = (rf(rr,zz) - 1.0)/e - 1e-9
        !print *, 'xmax=',xmax

        n=size(xgrid0)
        n1=size(xx) ! n1=num
        n2=size(xxx)! n2=num1
     

        allocate(cscoefsx(4,num-1),cscoefrx(4,num-1),cscoefrx1(4,num-1))

!      xtp < xmin
       
!        almax = asin(sqrt(1.0 - bnor*La)) ! max alpha at Xmin
!        almin = asin(sqrt(1.0 - bnor*Lb)) ! min alpha at Xmin



        almax = PIO2
        almin = -PIO2
        n3=size(aa0) ! alpha_k
        n4=size(xxk0) ! X_k for X <=0
        n5=size(xxk)  ! X_k for X >=0
        vbg = 0.0
        vbgp = 0.0
       
!      aa0: pitch angle,i.e. alpha=sigma_paral*sin^(-1)(abs(v_paral/v))
       ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.3, B.4
       ! note grid num n3, n4, n5 relate to mube grid numL
       ! xtp < xmin (B.3) grid of alpha
       call linspace(vbg, almax, n3, aa0) ! above mid-plane at x=xmin
       ! xtp > xmin (B.4) grid of X
       call linspace(xmin, -vbg, n4, xxk0) ! [xmin, 0.0]
       call linspace(vbgp, xmax, n5, xxk)   ! [vbgp, xmax]
!
       ! grids of X= (R-R0)/a for interp
       call linspace(xmin,vbg,n1,xx)
       call linspace(vbgp,xmax,n1,xx1)
       dl = 1e-8 ! error tol
       do j=1,n1
          dum1 = 1.0 + e*xx(j)
          if (xx(j) < 0.0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr1,dum1) ! x=r/a in flux coordinate
            rx(j) = x0                                ! from X
            !write(1003,*) xx(j),x0,abs((rf(rr1,x0)-1.0)/e-xx(j))
          else
            x0 = 0.0
            rx(j) = x0
            !write(1003,*) xx(j),x0,abs((rf(rr1,x0)-1.0)/e-xx(j))
          endif

       enddo

       call inrcsnak(xx,rx,cscoefrx) ! X<=0, given X, obtain x=r/a
                                     ! X = (R-R0)/a
      
!
       do j=1,n1
          dum1 = 1.0 + e*xx1(j)
          if (xx1(j) > 0.0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr,dum1) ! x=r/a in flux coordinate
            rx(j) = x0                               ! from X
            !write(1003,*) xx1(j),x0,abs((rf(rr,x0)-1.0)/e-xx1(j))
          else
            x0 = 0.0
            rx(j) = x0
            !write(1003,*) xx1(j),x0,abs((rf(rr,x0)-1.0)/e-xx1(j))
          endif
       enddo
       !write(1003,*) 'X,r/a end'
       call inrcsnak(xx1,rx,cscoefrx1) ! X>=0, given X, obtain x=r/a
!
    
       ! set up mube grid 
       ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.6
        allocate(xstart0(numL))
        bnor = b_nor(rr1,zz)
       do i=1,numL 
          ! set up mube grid at the axis alpha
          if (i <= n3) then                       
          lamdum0(i) = (1.0- sin(aa0(i))**2)/bnor ! mube corresponds to alpha at x=xmin
                                                 ! above (X,alpha) mid-plane 
          xstart0(i) = xmin
          ! set up mube grid at the axis X <=0
          elseif ( i > n3 .and. i<=(n4+n3)) then ! X <= 0 at alpha = 0 rr1=pi
          rdum = inrcsval(xx,cscoefrx,xxk0(i-n3),0) ! x=r/a
          lamdum0(i) = 1/b_nor(rr1,rdum) 
          xstart0(i) = xxk0(i-n3)
          ! set up mube grid at the axis X >=0
          elseif (i > (n3+n4)) then ! X >= 0 at alpha = 0  rr=0.0
          rdum = inrcsval(xx1,cscoefrx1,xxk(i-n3-n4),0) ! x=r/a
          lamdum0(i) = 1/b_nor(rr,rdum)
          xstart0(i) = xxk(i-n3-n4)
          endif

       enddo
!
      ! delete same value element at X=0
      ! lamdum(n4+n3)=lamdum(n4+n3+1) with xxk0=xxk=vbg=0 
      ! replace mdum(n4+n3) with lamdum(n4+n3+1) 
      ! by move the grids forward
      do i = n3 + n4 + 2, numL 
         lamdum0(i-1) = lamdum0(i)
         xstart0(i-1) = xstart0(i)
        ! write(1003,*) 'X0:',xstart0(i)
      end do
   
      numL = numL - 1
      ! delete same value ellement at (Xmin,alpha=0)
      ! lamdum(1) = lamdum(n3+1) with x=xmin
      ! by move the grid forward
      do i = n3+2, numL
         lamdum0(i-1) = lamdum0(i)
         xstart0(i-1) = xstart0(i)
        ! write(1003,*) 'X0:',xstart0(i)
      end do
       numL = numL - 1
      
      allocate(temp_lamdum0(numL), temp_xstart0(numL))

      ! 将 i <= n3 的元素放在最前面，并且 lamdum0 的元素顺序倒排
      k = 1
      do i = n3, 1, -1
         temp_lamdum0(k) = lamdum0(i)
         temp_xstart0(k) = xstart0(i)
         k = k + 1
      end do

     ! 将 i > n3 且 i <= n3 + n4 的元素紧接着放置
      do i = n3 + 1, n3 + n4
         temp_lamdum0(k) = lamdum0(i)
         temp_xstart0(k) = xstart0(i)
         k = k + 1
      end do

     ! 将 i > n3 + n4 的元素放在最后
      do i = n3 + n4 + 1, numL
        temp_lamdum0(k) = lamdum0(i)
        temp_xstart0(k) = xstart0(i)
        k = k + 1
      end do

     ! 将临时数组复制回原始数组
      lamdum0 = temp_lamdum0
      xstart0 = temp_xstart0
      !do i = 1, numL
      !   write(1003,*) 'x0,lam0:',xstart0(i),lamdum0(i)
      !enddo
      deallocate(temp_lamdum0, temp_xstart0)
!
     ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.7
     ! starting points with resetting lamdum = lamdum + dlam/2
     ! for xstart and mube by using root-finding
      allocate(xstart(numL-1))
      allocate(lamdum(numL-1),dlam(numL-1))     
     ! xtp < xmin  xstart=xmin 
      do i = 1, n3-1
         dlam(i) = lamdum0(i+1) - lamdum0(i)    ! i: [1,n3-1]
         lamdum(i) = lamdum0(i) + dlam(i) /2.0 ! i: [1,n3-1]
         xstart(i) = xstart0(i)
        !write(1003,*) 'xstart,1-lam*B:',xstart(i),1.0-lamdum(i)*b_nor(rr1,zz),i
      enddo
       ! xtp > xmin xstart from mube*bnor=1 
       ! advance an element by moving the n3+1th element to the n3th
      do i = n3+1, n3+n4-2 ! X <= 0 rr1=pi at alpha = 0 
         dlam(i-1) = lamdum0(i+1) - lamdum0(i)     
         lamdum(i-1) = lamdum0(i) + dlam(i-1) /2.0 
         rdum = inrcsval(xx,cscoefrx,xstart0(i),0) ! x=r/a
         x1=0.0
         x2=rdum
         ddx=(x2-x1)/5.0
         x0=(x2+x1)/1.0
         call secant3(dl,x0,ddx,istep,rr1,lamdum(i-1)) ! r/a corresponding to lamdum + dlam/2
         xstart(i-1) = (rf(rr1,x0) - 1.0)/e  ! rf(pi,r/a)
         ! write(1003,*) 'xstart,1-lam*B:',xstart(i-1),1.0-lamdum(i-1)*b_nor(rr1,x0),i-1
         !write(1003,*) 'xstart0:',xstart0(i),i
       enddo
       ! delete i=n3+n4-1 with X=0
       do i = n3+n4, numL-1
          dlam(i-2) = lamdum0(i) - lamdum0(i-1)    
          lamdum(i-2) = lamdum0(i-1) + dlam(i-2) /2.0   
          rdum = inrcsval(xx1,cscoefrx1,xstart0(i),0) ! fixed X=0,rdum=0,root-finding failed
          x1=0.0
          x2=rdum
          ddx=(x2-x1)/5.0
          x0=(x2+x1)/1.0 ! x2=0, root-finding failed
          call secant3(dl,x0,ddx,istep,rr,lamdum(i-2)) ! r/a corresponding to lamdum + dlam/2
          xstart(i-2) = (rf(rr,x0) - 1.0)/e  ! rf(0.0,r/a)
          !write(1003,*) 'xstart0:',xstart0(i),i
      enddo 

      numL = numL - 3 ! reset numL
      !do i = 1, numL
      ! write(1003,*) 'xstart:',xstart(i),'lam:',lamdum(i)
      !enddo
!
       ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.1
       ! uniform En grid
       call linspace(Ea,Eb,nE,en0)
       dE = en0(2)-en0(1)
       do i=1,nE-1
           en0(i) = en0(i) + dE/2.0
           write(8,443) en0(i), dE
       enddo
!
!     (E,mube,X) grid set up
      kk = 0
      do i=1,2  ! co-/ctr- direction
          if (i==1) then
             sigma = 1.0 ! co- direction
          else
             sigma = -1.0 ! ctr- direction
          endif
!
          
        do iE = 1, nE-1    ! energy grids 
          do j= 1,numL     ! pitch/mube grids
              call linspace(xstart(j),xmax,n1,xx2) ! X grids for xtp < xmin and xtp > xmin
              do k=1,num   ! X grids with equal dX
                 ! Given lambda, obtaining alphas at X grids
                 if (xx2(k) >=0) then
                      rdum = inrcsval(xx1,cscoefrx1,xx2(k),0) ! x=r/a
                      !write(1003,*) 'X', xx2(k),1+e*xx2(k)-rf(rr,rdum)
                 else
                      rdum = inrcsval(xx,cscoefrx,xx2(k),0)
                      !write(1003,*) 'X', xx2(k),1+e*xx2(k)-rf(rr1,rdum)
                 endif
                 alphadum(j,k)=fun_alpha(rdum,xx2(k),lamdum(j),sigma) ! pitch angle/ alpha curve/ grid line
                 write(7,442) xx2(k), alphadum(j,k), &
                              lamdum(j), en0(iE), sigma
                 !write(1003,*) 'X,r:',xx2(k),rdum
                 if (k>1) then
                    ! xxx of [Xstart, Xk] grids 
                    call linspace(xx2(1),xx2(k),n2,xxx)
                    dx = abs(xxx(2)-xxx(1))
                    do l=1,n2
                       if (xxx(l) >=0) then
                          rdum = inrcsval(xx1,cscoefrx1,xxx(l),0)
                       else
                          rdum = inrcsval(xx,cscoefrx,xxx(l),0)
                       endif

                       fi(l) = fun_s(rdum,xxx(l),lamdum(j),sigma) ! fi= sqrt(1+(dalpha/dx)^2)
                       !write(1003,*) 'rdum:', rdum, fi(l),xxx(l)
                    enddo
                    
                    fi(1) = -fi(3)+2.0*fi(2) ! fix Inf for dalpha/dx with 1-lam*B=0 
                    
                    call simp(n2,dx,fi,sdum(k))  ! arc length of
                                                 ! range [Xstart,xx2(k)]
                  else
                     sdum(k) = 0.0
                 endif
                 !write(1003,*) 's:', sdum(k),k,xx2(k)
              enddo
!
              ! arc length s corresponding to X
              call inrcsnak(sdum,xx2,cscoefsx) ! given s, obtain X
!             real(nx) should be real(num)?
              ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.10
              nn = max(2,floor(sdum(num)/(xmax-xmin)*real(nrho))+1)
              !
              allocate(ss(nn),xdum0(nn))
              call linspace(vbg,sdum(num),nn,ss)
              !
            ds = ss(2) - ss(1)
            do l = 1, nn
               xdum0(l) = inrcsval(sdum,cscoefsx,ss(l),0) ! X for equal arc len, used to calculate deltax x(xdum_i)-x(xdum_i-1)
               !write(1003,*) 'xdum0=',xdum0(l)
            enddo
              do l = 1, nn-1 ! X grids with equal arc length, should be nn-1
                 kk = kk + 1 ! ion index
                 tmps = ss(l) + ds/2.0
                 ! xdum: Xi at si+ds/2, more accouracy for f(P_z)dPz integral
                 ! note alpha=0,xdum0(1) on X axis
                 ! xdum(1) corresponding to alpha \=0 instead
                 ! accounting for negative and positive sigma
                 ! ejecting alpha=0 points
                 xdum = inrcsval(sdum,cscoefsx,tmps,0)
                 !xdum = xdum0(l) + (xdum0(l+1)-xdum0(l))/2.0 ! another way for set Xi
                 ! rdum: r/a at Xi 
                 if (xdum >=0) then
                   rdum = inrcsval(xx1,cscoefrx1,xdum,0)
                   thetadum = 0.0 ! initial position on X axis
                 else
                   rdum = inrcsval(xx,cscoefrx,xdum,0)
                   thetadum = PI_D ! initial position on X axis
                 endif
                 ! dx=x(Xi+1)-x(Xi)
                 if (xdum0(l+1)>=0) then
                    x1 = inrcsval(xx1,cscoefrx1,xdum0(l+1),0)
                    thetadum1 = 0.0
                 else
                    x1 = inrcsval(xx,cscoefrx,xdum0(l+1),0)
                    thetadum1 = PI_D
                 endif
                ! write(1003,*) xdum0(l+1),x1
                 if (xdum0(l)>=0) then
                    x2 = inrcsval(xx1,cscoefrx1,xdum0(l),0)
                    thetadum2 = 0.0
                 else
                    x2 = inrcsval(xx,cscoefrx,xdum0(l),0)
                    thetadum2 = PI_D
                 endif
                 !write(1003,*) xdum0(l),x2,xdum
                 dxdum = abs(x1-x2)
!
                 phidum = 0.0  ! initial position 
!
                 dum = fun_alpha(rdum,xdum,lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 dum4 = e*(b_norpx(thetadum,rdum)/&
                           rfpx(thetadum,rdum)+&
                           b_norpt(thetadum,rdum)/&
                           rfpt(thetadum,rdum))
                 !A.Bierwage,M.Fitzgerald,Ph.Lauberetal. ComputerPhysicsCommunications275(2022)108305 eq.(36)
                 dpzdum1 = (e**3*jr2(rdum)/rfpx(thetadum,rdum)-&
                          e/rf(thetadum,rdum)**2*&
                          (dum2/dum3/&
                          b_nor(thetadum,rdum)*lamdum(j)+&
                          dum3/b_nor(thetadum,rdum)**2))*&
                          (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
                 temp = b_nor(thetadum,rdum)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,rdum,0)
                 temp0 = tmpfpol*e**2 - dum3/temp  ! orbit_eq normalise for P_phi, and P_phi/e**2 normalise in KAEC
!
                 !dpzdum1 = (e**3*jr2(rdum)/rfpr(thetadum,rdum)/e+&
                 !         dum4*(dum2/dum3/&
                 !         b_nor(thetadum,rdum)*lamdum(j)+&
                 !         dum3/b_nor(thetadum,rdum)**2))*&
                 !         (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
!
                 dum = fun_alpha(x1,xdum0(l+1),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum1,x1)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x1,0)
!                particle toroidal canonical momentum normalized to
!                R0**2*omega_c
                 temp2 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl+1

                 dum = fun_alpha(x2,xdum0(l),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum2,x2)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x2,0)
!                particle toroidal canonical momentum 
                 temp1 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl
                 ! A.Bierwage,M.Fitzgerald,Ph.Lauberetal. ComputerPhysicsCommunications275(2022)108305 eq.(35)
                 dpzdum2 = abs(temp2 - temp1)       ! delta Pz, should be positive ?
!
                ! print *, 'dpzdum=',abs(dpzdum1), dpzdum2
                ! print *, 'E=', en0(iE), 'lambda=',lamdum(j)
!
                 write(6,441) rdum, thetadum, phidum, xdum, &
                              fun_alpha(rdum,xdum,lamdum(j),sigma), &
                              lamdum(j), temp0, en0(iE), dlam(j), &
                              dpzdum2, dE, kk, sigma
              enddo
              deallocate(ss)
              deallocate(xdum0)
          enddo
        enddo
      enddo
      close(6)
       deallocate(alpha)
       deallocate(alphadum)
       deallocate(rx)
       deallocate(en0)
       deallocate(lamdum)
       deallocate(aa0)
       deallocate(aa)
       deallocate(dlam)

      ik = 1 ! ik: 0 turn off fort.1900 writing
             ! ik: 1 for turn on fort.1900 writing
      a0 = 0.0
      a1 = 1.0
      tini = 0.0
      tfinal = 50000
      write(1900,733)
      if (ik.eq.1) then
         write(1900,444) kk,a0,a1,tini,tfinal
       else
          kk=5
          write(1900,444) kk,a0,a1,tini,tfinal
       endif
       close(1900)
      end subroutine check_mesh4
!
!     Bierwage, A. Computer Physics Communications, 5(183) 1107-1123 2012
!     Appendix B. Grid setup procedure
      subroutine check_mesh5 !  analytic equalibrium
        implicit none
        ! grids set up in (X, alpha, E)
        ! nx: X at alpha = 0, na: alpha at X=Xmin
        ! num: X for alpha as function of X 
        ! num1: X for arc length integral
        integer,parameter :: num=101,num1=21
        !common/aspect/e
        !rr=R/R0, zz=Z/R0=0
        !dE: grid of energy
        !xmin=(Rmin-R0)/a,xmax=(Rmax-R0)/a
        !Rmin:high field side major radius
        !Rmax:low filed side major radius
        !R0: magnetic axis major radius
        !a: minor radius
        !real(r8) :: e,lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: lam,x,rr,zz,sigma,sdum0,dx_ref,dx_ref0,da_ref,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,vbg,vbgp
        real(r8) :: almin, almax,rdum,sdum
        real(r8) :: tmpfpol,temp0,temp,temp1,temp2
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        real(r8) :: dl,  x0, ddx, thetadum, phidum
        real(r8) :: thetadum1,thetadum2, dpzdum1,dpzdum2
        real(r8) :: dum, dum1, dum2, dum3, dum4
        ! s: arc length for X at alpha = 0
        real(r8) ::  s(nrho-1)
        ! lamb/lamdum: mube
        ! xx: grids of [Xmin,Xmax], X=(R-R0)/rminor
        ! xxx: grids of the ith X grid
        real(r8) :: xgrid0(nrho+1),xgrid(nrho),lamb(nrho),&
                    xx(num),xx1(num),xx2(num),xxx(num1),fi(num1)
        ! aa: grids of alpha for X=Xmin, Z=0
        !real(r8) :: lamdum(na),aa0(na+1),aa(na),dlam(na)
        real(r8) :: tini,tfinal,a0,a1
        ! ss: arc length
        ! cscoefx: spline interp coefs between X and s
        real(r8), allocatable :: ss(:),xdum0(:),rx(:),xxk0(:),&
                                 xxk(:),xstart(:),xstart0(:),&
                                 temp_lamdum0(:), temp_xstart0(:),&
                                 xxk1(:),xxk2(:)
        real(r8), allocatable :: cscoefsx(:,:),cscoefrx(:,:),&
                                 cscoefrx1(:,:),cscoeflamx(:,:)
        real(r8), allocatable :: en0(:),lamdum0(:),lamdum(:),&
                                 aa0(:),aa(:),dlam(:)
        integer :: n,n1,n2,n3,n4,n5,i,j,k,m,l,nn,kk,kk1,iE,reverse_k
        integer :: nlam,istep,ik,numalp,numalp1,numx_ref0,numx_ref,numL1
!
730      format('# x, theta, phi, XR, alpha(pitch angle), &
               mube, Pz, E(keV), dL, dP, dE, kth, sigma')
         open(unit = 6,file = 'mesh_xal_arc.dat')
         write(6,730)
731      format('# XR, alpha(pitch angle), mube, E(keV) sigma')
         open(unit = 7,file = 'mesh_xal.dat')
         write(7,731)
732      format('# Ei(keV)')
         open(unit = 8, file = 'grid_E.dat')
         write(8,732)
733      format('# nprt a0 a1 tini tfinal')

441      format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                1E16.9,1x,I8,1x,1F12.6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
444      format(I8,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
448      format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                  1E16.9,1x,1E16.9,1x,I6)

!
        allocate(alpha(nrho,num),alphadum(numL,num),rx(num))
        allocate(en0(nE))
        
        numalp = numL - nrho ! grids of alpha above mid-plane
                             ! numL: total mube grids
        !numalp1 = 2*numalp - 1
        allocate(aa0(numalp-1)) ! alpha_k at xmin
        allocate(aa(numalp-1)) ! alpha_k+1/2 at xmin
       
        
        !write(1003,*) 'numalp:',numalp,'kk:',kk,'kk1:',kk1
        rr = 0.0
        rr1 = PI_D
        zz = 1.0 ! LCFS flux
!
        !adum = (Rmax1 -Rmin1)/2.0
        adum = (rf(rr,zz) - rf(rr1,zz))/2.0
        !xmin = (Rmin1 -rmagx)/rminor
        xmin = (rf(rr1,zz)-1.0)/e + 1e-9 ! normalized to minor radius a
        !print *, 'xmin=',xmin
        !xmax = (Rmax1 - rmagx)/rminor
        xmax = (rf(rr,zz) - 1.0)/e - 1e-9 ! normalized to minor radius a
        !print *, 'xmax=',xmax

        n=size(xgrid0)
        n1=size(xx) ! n1=num
     

        allocate(cscoefsx(4,num-1),cscoefrx(4,num-1),cscoefrx1(4,num-1))

!      xtp < xmin
       
!        almax = asin(sqrt(1.0 - bnor*La)) ! max alpha at Xmin
!        almin = asin(sqrt(1.0 - bnor*Lb)) ! min alpha at Xmin



        almax = PIO2
        almin = -PIO2
        n3=size(aa0) ! alpha_k
        vbg = 0.0
        vbgp = 0.0


       
!      aa0: pitch angle,i.e. alpha=sigma_paral*sin^(-1)(abs(v_paral/v))
       ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.3, B.4
       ! note grid num n3, n4, n5 relate to mube grid numL
       
       ! xtp < xmin (B.3) grid of alpha,  [1, numalp-1] for k
       da_ref = almax/(numalp - 1)      
       do k = 1, numalp-1
          ! Descending Order for alpha due to  Ascending Order for mube
          reverse_k = numalp - k 
          aa0(k)=almax*abs(real(reverse_k)*&
                 da_ref/almax)**0.8*&
                 real(reverse_k)/abs(real(reverse_k))  ! alpha_k at xmin
          aa(k)=almax*abs((real(reverse_k)-0.5)*&
                 da_ref/almax)**0.8*&
                 real(reverse_k)/abs(real(reverse_k))  ! alpha_k+1/2 at xmin
                                                       ! for 'mube_k+1/2' 
                                                       ! note mube_k+1/2 is not the centers of the grid mube
                                                       ! use alpha_k+1/2 here and X_k+1/2 below to 
                                                       ! obtain the corresponding 'mube_k+1/2' 
                                                       ! since hard to interp mube_k+1/2 through X_k+1/2 
       end do 
       
       ! xtp > xmin 
       ! grids of X= (R-R0)/a for x=r/a interp 
       call linspace(xmin,vbg,n1,xx)
       call linspace(vbgp,xmax,n1,xx1)
       dl = 1e-8 ! error tol
       ! sp coefs for X<=0
       do j=1,n1
          dum1 = 1.0 + e*xx(j)
          if (xx(j) < 0.0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr1,dum1) ! x=r/a in flux coordinate
            rx(j) = x0                             ! from X
            !write(1003,*) xx(j),x0,abs((rf(rr1,x0)-1.0)/e-xx(j))
          else
            x0 = 0.0
            rx(j) = x0
            !write(1003,*) xx(j),x0,abs((rf(rr1,x0)-1.0)/e-xx(j))
          endif

       enddo

       call inrcsnak(xx,rx,cscoefrx) ! X<=0, given X(xx), obtain x=r/a(rx)
                                     ! with X = (R-R0)/a
      
!
       ! sp coefs for X>=0
       do j=1,n1
          dum1 = 1.0 + e*xx1(j)
          if (xx1(j) > 0.0) then
            x1=0.0
            x2=1.0
            ddx=(x2-x1)/5.0
            x0=(x2+x1)/2.0
            call secant2(dl,x0,ddx,istep,rr,dum1) ! x=r/a in flux coordinate
            rx(j) = x0                            ! from X
            !write(1003,*) xx1(j),x0,abs((rf(rr,x0)-1.0)/e-xx1(j))
          else
            x0 = 0.0
            rx(j) = x0
            !write(1003,*) xx1(j),x0,abs((rf(rr,x0)-1.0)/e-xx1(j))
          endif
       enddo
      
       call inrcsnak(xx1,rx,cscoefrx1) ! X>=0, given X(xx1), obtain x=r/a(rx)
       
!     set up X grids
!A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.4
       numx_ref = nrho
       numx_ref0 = max(3, ceiling(real(1.0 * numx_ref)))
       
      
       if (allocated(xxk1)) then
          deallocate(xxk1)
       endif
       allocate(xxk1(0:numx_ref0-1))
       if (allocated(xxk)) then
          deallocate(xxk)
       endif 
       allocate(xxk(0:numx_ref0-2))
       if (allocated(xxk0)) then
          deallocate(xxk0)
       endif 
       allocate(xxk0(0:numx_ref-1))
       dx_ref0 = (xmax - xmin)/real(numx_ref0 - 1)
       do k = 0,numx_ref0 - 1
          dum1 = abs(real(k)* dx_ref0/(xmax-xmin))**0.8
          xxk1(k) = xmin + dum1*(xmax-xmin) ! (B.4) X_k at alpha = 0
          if (k < numx_ref0 - 1) then
             dum1 = abs((real(k)+0.5)* dx_ref0/(xmax-xmin))**0.8 
             xxk(k) = xmin + dum1*(xmax-xmin) ! X_k+1/2 for alpha line
                                              ! note total X_k+1/2 grids = total X_k grids - 1
          endif
       enddo 
       
!
    
       ! set up mube grid for delta_mube
       ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.6
        numL1 = numalp+numx_ref0 - 1 ! total grids of (B.6)
        allocate(xstart0(numL1))
        allocate(lamdum0(numL1))
        ! first expression (B.6) 
        bnor = b_nor(rr1,zz)
        !! set up mube grid at the axis alpha
        do i=1,numalp - 1                              
          lamdum0(i) = (1.0- sin(aa0(i))**2)/bnor ! mube corresponds to alpha at x=xmin
                                                  ! mube with Ascending Order
          xstart0(i) = xmin
        enddo
        ! second expression (B.6)
        !! set up mube grid at the axis X
        do i = numalp, numL1
           if (xxk1(i-numalp) <= 0) then  !X <= 0 at alpha = 0 rr1=pi
              if (i==numalp) then  ! at xmin
                 rdum = 1.0
              else
                 rdum = inrcsval(xx,cscoefrx,xxk1(i-numalp),0) ! x=r/a
              endif
              lamdum0(i) = 1/b_nor(rr1,rdum) 
              xstart0(i) = xxk1(i-numalp)
           else                          ! X >= 0 at alpha = 0  rr=0.0
              if (i==numL1) then  ! at xmax
                  rdum = 1.0
              else
                  rdum = inrcsval(xx1,cscoefrx1,xxk1(i-numalp),0) ! x=r/a
              endif
              lamdum0(i) = 1/b_nor(rr,rdum)
              xstart0(i) = xxk1(i-numalp)      
           endif
        enddo
!

!
     ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.7
     ! obtain delta_ mube from (B.7) and mube_i+1/2, X_i+1/2
      allocate(xstart(numL1-1))
      allocate(lamdum(numL1-1),dlam(numL1-1))     
     ! xtp < xmin  xstart=xmin  at the axis of alpha
     ! delta_mube is not const. in general
      bnor = b_nor(rr1,zz)
      do i = 1, numalp-1
            if (i < numalp -1) then
               dlam(i) = lamdum0(i+1) - lamdum0(i)
            else
               dlam(i) = 1/bnor - lamdum0(i) ! mube(alpha = 0) - mube(alpha_numalp-1) 
            endif    
            lamdum(i) = (1.0- sin(aa(i))**2)/bnor  ! mube_i+1/2 from alpha_i+1/2
            xstart(i) = xstart0(i)   ! starting point for X_i+1/2, x=xmin
            !write(1003,*) 'xstart,1-lam*B:',xstart(i),1.0-lamdum(i)*b_nor(rr1,zz),i
      enddo
      ! xtp > xmin xstart from mube*bnor=1 at the axis of X
      ! alpha = 0
      do i = numalp, numL1-1 ! total grids should be -1 due to X_k+1/2 grids
            dlam(i) = lamdum0(i+1) - lamdum0(i) 
            xstart(i) = xxk(i-numalp)  ! starting x point,X_k+1/2, alpha = 0
            if (xxk(i-numalp) <= 0) then  !X <= 0 at alpha = 0 rr1=pi
              rdum = inrcsval(xx,cscoefrx,xxk(i-numalp),0) ! x=r/a
              lamdum(i) = 1/b_nor(rr1,rdum) ! mube_i+1/2 from X_i+1/2 
                                            ! not use mube_i+dmube/2
           else                          ! X >= 0 at alpha = 0  rr=0.0
              rdum = inrcsval(xx1,cscoefrx1,xxk(i-numalp),0) ! x=r/a
              lamdum(i) = 1/b_nor(rr,rdum)   ! mube_i+1/2 from X_i+1/2  
                                             ! not use mube_i+dmube/2           
           endif                 
      enddo

      numL = numL1 - 1 ! reset numL
!
      
       ! A. Bierwage et al. / Computer Physics Communications 183 (2012) 1107–1123 B.1
       ! uniform En grid
       
       call linspace(Ea,Eb,nE,en0)
       dE = en0(2)-en0(1)
       do i=1,nE-1
          ! A.Bierwage et al. Computer Physics Communications 275 (2022) 108305
          ! 4. CoM mesh and drift orbit types
           en0(i) = en0(i) + dE/2.0
           write(8,443) en0(i), dE
       enddo
      
!
    
    
!     set up (E,mube,X) grids along alpha lines
      kk = 0
       do i=1,2  ! co-/ctr- direction
           if (i==1) then
              sigma = 1.0 ! co- direction
           else
              sigma = -1.0 ! ctr- direction
           endif
! !
          
         do iE = 1, nE-1    ! energy grids 
           do j= 1,numL     ! pitch/mube grids
           
               
               !call linspace(xstart(j),xmax,num,xx2) !debug
               call linspace(xstart(j),xmax,n1,xx2) ! X grids for xtp < xmin or xtp > xmin
               
              
              
!               ! Given lambda, obtaining X grids along alpha lines, (E,mube,X)
!               ! X grids with equal dX
               do k=1,n1   
                
                  if (xx2(k) >=0) then
                       if (abs(xmax-xx2(k))<1e-6) then
                          rdum = 1.0
                      else
                          rdum = inrcsval(xx1,cscoefrx1,xx2(k),0) ! x=r/a
                      !write(1003,*) 'X', xx2(k),1+e*xx2(k)-rf(rr,rdum)
                      endif
                 else
                       if (abs(xmin-xx2(k))<1e-6) then
                          rdum = 1.0
                       else
                          rdum = inrcsval(xx,cscoefrx,xx2(k),0)
                      !write(1003,*) 'X', xx2(k),1+e*xx2(k)-rf(rr1,rdum)
                       endif
                 endif
                 alphadum(j,k)=fun_alpha(rdum,xx2(k),lamdum(j),sigma) ! pitch angle/ alpha curve/ grid line
                 write(7,442) xx2(k), alphadum(j,k), &
                              lamdum(j), en0(iE), sigma
               enddo ! X grids with equal dX
                 
              !
              ! A.Bierwage et al. Computer Physics Communications 275 (2022) 108305 
              ! Given lambda, obtaining X grids along alpha lines, B.8-B.10, 
              ! (E,mube,X) and dPz, dL, dE
              ! A.Bierwage et al. Computer Physics Communications 275 (2022) 108305 
              ! obtain grids of alpha lines, N_D=kk1 in B.10
              dx_ref = (xmax - xstart(j))/(numx_ref - 1)                  
              do l = 0,numx_ref - 1         
                 dum1 = abs(real(l)* dx_ref/(xmax-xstart(j)))**1.0
                 xxk0(l) = xstart(j) + dum1*(xmax-xstart(j)) ! X_k (B.9) 
              enddo 
             
!               ! obtain the length of alpha line, D(mube)=sdum in B.9
              sdum = 0.0
              do l = 1, numx_ref - 1
                 dum1 = xxk0(l) - xxk0(l-1) ! dX_k (B.9)
                 dum1 = dum1 **2 
                 if (xxk0(l) >=0) then
                    if (abs(xmax-xxk0(l))<1e-6) then
                        rdum = 1.0
                    else
                        rdum = inrcsval(xx1,cscoefrx1,xxk0(l),0) ! x=r/a
                    
                    endif
                 else
                    if (abs(xmin-xxk0(l))<1e-6) then
                        rdum = 1.0
                    else
                        rdum = inrcsval(xx,cscoefrx,xxk0(l),0)
                     
                    endif
                 endif
                 dum2 = fun_alpha(rdum,xxk0(l),lamdum(j),sigma) ! alpha(X_k) (B.9)
                 
                 if (xxk0(l-1) >=0) then
                    if (abs(xmax-xxk0(l-1))<1e-6) then
                       rdum = 1.0
                    else
                       rdum = inrcsval(xx1,cscoefrx1,xxk0(l-1),0) ! x=r/a
                    
                    endif
                 else
                    if (abs(xmin-xxk0(l-1))<1e-6) then
                       rdum = 1.0
                    else
                       rdum = inrcsval(xx,cscoefrx,xxk0(l-1),0)
                      
                    endif
                 endif
                 dum3 = fun_alpha(rdum,xxk0(l-1),lamdum(j),sigma)  ! alpha(X_k-1) (B.9)
                 sdum = sdum + sqrt(dum1 + (dx_ref/da_ref*(dum2-dum3))**2) ! D(mube) (B.9)
              enddo
              kk1 = max(2, ceiling(sdum/(xmax-xmin)*real(numx_ref))) ! (B.10)
              
              
              ! set up grids based on the length of alpha line
              ! to obtain points Fig. B.12
              ! avoid reallocate
              if (allocated(xdum0)) then
                  deallocate(xdum0)
              endif
              allocate(xdum0(0:kk1-1)) 
              dx = (xmax - xstart(j)) / real(kk1 - 1)  ! reset dX
              do l = 0, kk1 - 1         
                 xdum0(l) = xstart(j) + dx*real(l) ! points Fig. B.12
                                                     ! for uniform X_i
                                                     ! based on the length of alpha line
              enddo 
!   
             ! xdum: Xn at n=i+1/2, more accouracy for f(P_z)dPz integral
             ! note alpha=0, xdum0(0) on X axis
             ! xdum(1) corresponding to alpha \=0 instead
             ! accounting for negative and positive sigma
             ! ejecting xdum for alpha=0 points
                               
             do l = 1, kk1-1    ! loop from 1 instead of 0       
                kk = kk + 1 ! ion index        
                xdum = xdum0(l) - dx/2.0 ! X_i-1/2 along alpha line start from i=1 
                                           ! at alpha axis, x=xmin?
                if (xdum >=0) then
                   thetadum = 0.0 ! initial position on X axis
                   if (abs(xmax-xdum)<1e-6) then
                          rdum = 1.0
                   else
                          rdum = inrcsval(xx1,cscoefrx1,xdum,0) ! x=r/a
                      
                   endif
                else
                    thetadum = PI_D ! initial position on X axis
                    if (abs(xmin-xdum)<1e-6) then
                          rdum = 1.0
                    else
                          rdum = inrcsval(xx,cscoefrx,xdum,0)
                      
                    endif
                endif
                 
                 ! obtain x(X_i),x(X_i-1) for dx=x(X_i)-x(X_i-1)
                 if (xdum0(l)>=0) then
                   if (abs(xmax-xdum0(l))<1e-6) then
                          x1 = 1.0
                   else
                          x1 = inrcsval(xx1,cscoefrx1,xdum0(l),0)
                   endif
                    thetadum1 = 0.0
                 else
                    if (abs(xmin-xdum0(l))<1e-6) then
                          x1 = 1.0
                   else
                          x1 = inrcsval(xx,cscoefrx,xdum0(l),0)
                   endif
                    
                    thetadum1 = PI_D
                 endif
               
                 if (xdum0(l-1)>=0) then
                   if (abs(xmax-xdum0(l-1))<1e-6) then
                          x2 = 1.0
                   else
                          x2 = inrcsval(xx1,cscoefrx1,xdum0(l-1),0)
                   endif
                    thetadum2 = 0.0
                 else
                   if (abs(xmin-xdum0(l-1))<1e-6) then
                          x2 = 1.0
                   else
                          x2 = inrcsval(xx,cscoefrx,xdum0(l-1),0)
                   endif
                    thetadum2 = PI_D
                 endif
                 
                 dxdum = abs(x1-x2)
!
                 phidum = 0.0  ! initial position 
!
                 dum = fun_alpha(rdum,xdum,lamdum(j),sigma)
                 dum1 = sin(dum) ! vparal/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1 ! vparal
                 dum4 = e*(b_norpx(thetadum,rdum)/&
                           rfpx(thetadum,rdum)+&
                           b_norpt(thetadum,rdum)/&
                           rfpt(thetadum,rdum))
                 ! A.Bierwage et al. Computer Physics Communications 275 (2022) 108305 eq.(36)
                 dpzdum1 = (e**3*jr2(rdum)/rfpx(thetadum,rdum)-&
                          e/rf(thetadum,rdum)**2*&
                          (dum2/dum3/&
                          b_nor(thetadum,rdum)*lamdum(j)+&
                          dum3/b_nor(thetadum,rdum)**2))*&
                          (xdum0(l)-xdum0(l-1))    ! analytic delta Pz
                 temp = b_nor(thetadum,rdum)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,rdum,0) ! pol
                 ! P_phi at X_i-1/2 along alpha line
                 temp0 = tmpfpol*e**2 - dum3/temp  
                 ! orbit_eq normalise for P_phi, and P_phi/e**2 normalise in KAEC
                 ! particle toroidal canonical momentum normalized to
!                ！R0**2*omega_c in orbit_eq
                 ! a**2*omega_c in KAEC
!
                 !dpzdum1 = (e**3*jr2(rdum)/rfpr(thetadum,rdum)/e+&
                 !         dum4*(dum2/dum3/&
                 !         b_nor(thetadum,rdum)*lamdum(j)+&
                 !         dum3/b_nor(thetadum,rdum)**2))*&
                 !         (xdum0(l+1)-xdum0(l))    ! analytic delta Pz
!
                 ! P_phi at Xl
                 dum = fun_alpha(x1,xdum0(l),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum1,x1)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x1,0)

                 temp2 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl
                 !
                 ! P_phi at Xl-1
                 dum = fun_alpha(x2,xdum0(l-1),lamdum(j),sigma)
                 dum1 = sin(dum) ! vpara/v
                 dum2 = rhoh_bar**2*en0(iE)/ekev
                 dum3 = sqrt(2.0*dum2)*dum1
                 temp = b_nor(thetadum2,x2)
                 tmpfpol = inrcsval(gridx2,cscoefxpsi,x2,0)
!                particle toroidal canonical momentum 
                 temp1 = tmpfpol*e**2 - dum3/temp   ! P_phi at Xl-1
                 ! A.Bierwage et al. Computer Physics Communications 275 (2022) 108305 eq.(35)
                 dpzdum2 = abs(temp2 - temp1)       ! dPz, should be positive ?
!
                ! print *, 'dpzdum=',abs(dpzdum1), dpzdum2
                ! print *, 'E=', en0(iE), 'lambda=',lamdum(j)
!
                 write(6,441) rdum, thetadum, phidum, xdum, &
                              fun_alpha(rdum,xdum,lamdum(j),sigma), &
                              lamdum(j), temp0, en0(iE), (dlam(j)), &
                              dpzdum2, dE, kk, sigma
             enddo
             
           enddo      ! pitch/mube grids
         enddo        ! energy grids 
       enddo          ! co-/ctr- direction
      close(6)
       deallocate(alpha)
       deallocate(alphadum)
       deallocate(rx)
       deallocate(en0)
       deallocate(lamdum)
       deallocate(aa0)
       deallocate(aa)
       deallocate(dlam)
       !
       deallocate(xxk)
       deallocate(xxk0)
       deallocate(xxk1)
       deallocate(cscoefsx)
       deallocate(cscoefrx)
       deallocate(cscoefrx1)
       deallocate(xstart0)
       deallocate(xstart)
       deallocate(lamdum0)
       
       

      ik = 1 ! ik: 0 turn off fort.1900 writing
             ! ik: 1 for turn on fort.1900 writing
      a0 = 0.0
      a1 = 1.0
      tini = 0.0
      tfinal = 50000
      write(1900,733)
      if (ik.eq.1) then
         write(1900,444) kk,a0,a1,tini,tfinal
       else
          kk=5
          write(1900,444) kk,a0,a1,tini,tfinal
       endif
       close(1900)
      end subroutine check_mesh5
!


!
      subroutine check_flux1 ! given q=q(x=r/a), obtain psi normailzed to a**2*B(0)
        implicit none        ! range of [0,pw] with pw normalized to a**2*B(0)
        integer,parameter :: num=101,num1=101
        !common/aspect/e
        !real(r8) :: e,lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,temp,vbg
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        real(r8) ::  sdum(num)
        ! xx: grids of [Xmin,Xmax], X=(R-R0)/rminor
        ! xxx: grids of the ith X grid 
        real(r8) ::  xx(num),xxx(num1),fi(num1)
        ! cscoefxpsi: spline interp coefs between x and psi
        integer :: n,n1,n2,n3,i,j,k,m,l,nn,kk,iE
        integer :: nlam
!
730      format('# x, psi')
         open(unit = 60,file = 'psix.dat')
         write(60,730)

731      format('# x, jr2')
         open(unit = 61,file = 'xjr2.dat')
         write(61,731)

441      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,&
               1F12.6,1x,1F12.6,1x,I6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
!
        !xmin = 0.0
        !xmax = 0.5
        !call linspace(xmin,xmax,num,xx)
!
         if (allocated(cscoefxpsi)) then
            deallocate(cscoefxpsi)
          endif
        allocate(cscoefxpsi(4,num-1))
       if (allocated(gridx2)) then
            deallocate(gridx2)
        endif
        allocate(gridx2(num))

       if (allocated(cscoefxjr2)) then
            deallocate(cscoefxjr2)
       endif
       ! allocate(cscoefxjr2(4,num-1))       

         if (allocated(cscoefxnh)) then
            deallocate(cscoefxnh)
          endif
        allocate(cscoefxnh(4,num-1))

       if (allocated(gridpsi2)) then
            deallocate(gridpsi2)
        endif
        allocate(gridpsi2(num))

       if (allocated(cscoefxjr2)) then
            deallocate(cscoefxjr2)
       endif
              
        n1=size(gridx2) ! n1 = num
        n2=size(xxx) ! n2 = num1

        xmin = 0.0
        xmax = 1.0
       call linspace(xmin,xmax,n1,gridx2)
!
       da = TWOPI_D/real(n2-1)
       dxdum = 1.0/real(n1-1)
!        
        !do k=1,n1
        !   rr1 = dxdum*(k-1)
        !  if (k>1) then
        !    do l=1,n2
        !      temp = da*(l-1)
        !      fi(l)=jr20(temp,rr1)
              !fi(l)=jr2(rr1)
        !    enddo
        !      call simp(n2,da,fi,sdum(k))
        !  else
        !     sdum(k) = 0.0
        !  endif
        !     write(61,443) rr1, sdum(k)/TWOPI_D
             !write(61,443) rr1, delx1(rr1)
        !enddo
        sdum = 0.0
!
        do k=1,num
           if (k>1) then
           ! xxx of grids [Xmin, Xk]
              call linspace(gridx2(1),gridx2(k),n2,xxx)
              dx = abs(xxx(2)-xxx(1))
              do l=1,num1
                 fi(l) = jr2(xxx(l))/qfun(xxx(l))
              enddo
               call simp(n2,dx,fi,sdum(k))
             else
                sdum(k) = 0.0
             endif
             write(60,443) gridx2(k), sdum(k)
         enddo
!
         ! calculate coef of psi on uniform x grid points
         call inrcsnak(gridx2,sdum,cscoefxpsi)
         pw = sdum(num) ! normaized to a**2*B(0)
         do k=1,num
            gridpsi2(k) = sdum(k)/pw
            xxx(k) = den_h(gridpsi2(k))
         enddo
         call inrcsnak(gridpsi2,xxx,cscoefxnh)
 
      end subroutine check_flux1

!
      subroutine check_flux2 ! given q=q(psi), obtain psi normaized 
        implicit none        ! to a**2*B(0) too
        integer,parameter :: num=101,num1=101
        real(r8) :: lam,x,rr,zz,sigma,sdum0,da,dE,ds
        real(r8) :: xmin,xmax,adum,bnor,dx,xdum,temp,vbg
        real(r8) :: tmps, dxdum, rr1, rr2, pol1, pol2, x1, x2
        real(r8) ::  sdum(num)
        real(r8) ::  xx(num),xxx(num1),fi(num1)
        integer :: n,n1,n2,n3,i,j,k,m,l,nn,kk,iE
        integer :: nlam
!
730      format('# x psi')
         open(unit = 60,file = 'psix.dat')
         write(60,730)

731      format('# x, jr2')
         open(unit = 61,file = 'xjr2.dat')
         write(61,731)

441      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,&
               1F12.6,1x,1F12.6,1x,I6)
442      format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443      format(1F12.6,1x,1F12.6)
!
         if (allocated(cscoefxpsi)) then
            deallocate(cscoefxpsi)
          endif
        allocate(cscoefxpsi(4,num-1))
       if (allocated(gridx2)) then
            deallocate(gridx2)
        endif
        allocate(gridx2(num))
        
        
              
        n1=size(gridx2) ! n1 = num
        n2=size(xxx) ! n2 = num1

        xmin = 0.0
        xmax = 1.0
!
       da = TWOPI_D/real(n2-1)
       dxdum = 1.0/real(n1-1)
!        
        do k=1,n1
           rr1 = dxdum*(k-1)
          if (k>1) then
            do l=1,n2
              temp = da*(l-1)
              fi(l)=jr20(temp,rr1)
            enddo
              call simp(n2,da,fi,sdum(k))
          else
             sdum(k) = 0.0
          endif
             write(61,443) rr1, sdum(k)/TWOPI_D
        enddo
        sdum = 0.0

       call linspace(xmin,xmax,n1,sdum) ! uniform psi grids with range
                                        ! of [0,1]
!        
!       q=1.6667+0.5*psi+0.8333*psi**2, psi=psi_p/psi_wall range of [0,1]
!       dpsi_bar/dx=jr2/q
!       int_0^\psi q(psi)dpsi*dpsi_bar/dpsi=int_0^x jr2(x)dx 
!        x = sqrt((3.3334*psi+0.5*psi**2+1.6666/3.0*psi**3)*pw)

        !pw = 0.5/(1.6667+0.25+0.8333/3.0)
        pw = 1.5/6.5834 ! normaized to a**2*B(0)

        do k=1,num  
           xdum = sdum(k) ! psi       
           temp = (3.3334*xdum+0.5*xdum**2+1.6666/3.0*xdum**3)*pw
           gridx2(k) = sqrt(temp) ! x=r/a
            write(60,443) gridx2(k), sdum(k)
         enddo
!
         ! calculate coef of psi on non-uniform x grid points
         sdum = sdum*pw
         call inrcsnak(gridx2,sdum,cscoefxpsi)
        ! cscoefxpsi: spline interp coefs between x and psi

        !do k=1,num
        !   temp = inrcsval(gridx2,cscoefxpsi,sdum(k),0) ! psi  
        !   write(4000,*) sdum(k),temp
        !enddo
      end subroutine check_flux2
!

      subroutine check3 ! record flux surfaces for straight field lines coordinates
        implicit none
        real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
        integer,parameter :: n=201,na=361
        !common/aspect/e!aspect rate:e=a/R_0
        real(r8) :: R(na),Z(na)
        integer :: i,j
        !real(r8) :: dt,dx,dum1,dum2,e
        real(r8) :: dt,dx,dum1,dum2
        real(r8) :: x(n),t(na)
        dt=two_pi/(na-1)
        dx=1.0/(n-1)

        do i=1,n
           x(i)=(i-1)*dx
           if (x(i)==0) then
              x(i)=1e-03
           endif
        enddo

        do i=1,na
           t(i)=(i-1)*dt
        enddo

        do i=1,n
          do j=1,na
          R(j)=rf(t(j),x(i))
          Z(j)=zf(t(j),x(i))
         if (mod(i,10)==0) then
          write(40,*) R(j),Z(j)
         endif
          enddo
       enddo

        do i=1,n
          do j=1,na
          R(j)=rf(t(j),x(i))
          Z(j)=zf(t(j),x(i))
           if (mod(j,12)==0) then
              write(40,*) R(j),Z(j)
           endif
          enddo
       enddo
!
      end subroutine check3
!
       subroutine check_xy  ! interp R(x,y),Z(x,y) 
        implicit none
        !common/aspect/e!aspect rate:e=a/R_0
        real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
        integer,parameter :: nx=201,ny=201,na=1801
        integer :: i,j,k,l,m,n,istep
        real(r8),  allocatable :: R_old(:,:), Z_old(:,:), theta1(:,:)
        real(r8),  allocatable :: theta(:), dum4(:), dum5(:)
        !real(r8) :: e,dy,dx,dx1,dt
        real(r8) :: dy,dx,dx1,dt
        real(r8) :: a,b,dl,ddx,x0
        real(r8) :: rdum, tdum, dum, dum1, dum2, dum3
        allocate(gridx(nx),gridy(ny),knotx(nx+kr),knoty(ny+kz))
        allocate(R_old(nx,ny),Z_old(nx,ny))
        allocate(bscoefxy2R(nx,ny),bscoefxy2Z(nx,ny))
        allocate(theta(na),dum4(na),dum5(na))
        allocate(theta1(nxefit,nyefit), bscoeftheta(nxefit,nyefit))
        dt = two_pi/real(na-1)
        do i=1,na
           theta(i) = (i-1)*dt
        end do
!
        dx = 3.0/real(nx-1)
        dy = 3.0/real(ny-1)
        do i=1,nx
                gridx(i) = -1.5 + dx*(i-1) + 5e-6
                if (gridx(i)==0) then
                    gridx(i)=5e-05
                endif
                !write(23,*) 'x=', gridx(i),i
        end do
        do j=1,ny
                gridy(j) = -1.5 + dy*(j-1)+5e-6
                if (gridy(j)==0) then
                    gridy(j)=5e-05
                endif
               ! write(23,*) gridy(j)
        end do

        do i = 1, nx
                do j = 1, ny
                       rdum = sqrt(gridx(i)**2+gridy(j)**2)
                       if (gridx(i) > 0 .and. gridy(j) > 0) then
                          tdum = atan(gridy(j)/gridx(i))
                       elseif (gridx(i) < 0 .and. gridy(j) > 0) then
                          tdum = atan(gridy(j)/gridx(i)) + pi
                       elseif (gridx(i) < 0 .and. gridy(j) < 0) then
                          tdum = atan(gridy(j)/gridx(i)) + pi
                       else
                          tdum = atan(gridy(j)/gridx(i)) + two_pi
                       end if
                         
                        R_old(i,j) = rf(tdum,rdum)
                        Z_old(i,j) = zf(tdum,rdum)

                end do
        end do
!
        call cdbbsnak(gridx,kr,knotx) !生成4阶r节点
        call cdbbsnak(gridy,kz,knoty) !生成4阶z节点
        call cdbbscoef2d(gridx,gridy,R_old,knotx,knoty,kr,kz,bscoefxy2R) ! R
        call cdbbscoef2d(gridx,gridy,Z_old,knotx,knoty,kr,kz,bscoefxy2Z) ! Z
       !print *,"Rmin,Rmax,Zmin,Zmax",Rmin1,Rmax1,Zmin1,Zmax1
       do i=1,nxefit
           do j=1,nyefit   
              dum = cdbbsval2d(knotr,knotz,kr,kz,bscoefpsi,gridr(i),gridz(j),0,0)
              dum1 = gridr(i)*rmag
              dum2 = gridz(j)*rmag
              dum3 = sqrt(dum/pw)! r/a
              !print *, dum1, dum2, dum3
              if ((dum1 >= Rmin1 .and. dum1 <= Rmax1) .and. &
                 (dum2 >= Zmin1 .and. dum2 <= Zmax1).and. dum3 <=1.0) then
                ! print *,"R,Z=",dum1,dum2
                ! print *, 'r/a=', dum3,i,j
                 dum1 = dum1 - rmag
                 dl = 1.0e-6
                 if (dum1>=0 .and. dum2 >=0) then
                     a = 0.0
                     b = pi/2.0
                     ddx = (b-a)/10.0
                     x0 = (a+b)/2.0
                     call secant(dl,x0,ddx,istep,dum3,i,j)
! 
                     !WRITE (23,"(I4,4F16.8)") ISTEP,X0,DDX,dum1+rmag,dum2
                     theta1(i,j) = X0
!
                 elseif (dum1 < 0 .and. dum2 > 0) then
                     a = pi/2.0
                     b = pi
                     ddx = (b-a)/10.0
                     x0 = (a+b)/2.0
                     call secant(dl,x0,ddx,istep,dum3,i,j)
!
                     !WRITE (23,"(I4,4F16.8)") ISTEP,X0,DDX,dum1+rmag,dum2
                     theta1(i,j) = X0
!
                 else if (dum1 < 0 .and. dum2 < 0) then
                     a = pi
                     b = 3.0*pi/2.0
                     ddx = (b-a)/10.0
                     x0 = (a+b)/2.0
                     call secant(dl,x0,ddx,istep,dum3,i,j)
!
                     !WRITE (23,"(I4,4F16.8)") ISTEP,X0,DDX,dum1+rmag,dum2
                     theta1(i,j) = X0
!
                 else 
                     a = 3.0*pi/2.0
                     b = two_pi
!
                     ddx = (b-a)/10.0
                     x0 = (a+b)/2.0
                     call secant(dl,x0,ddx,istep,dum3,i,j)
!
                     !WRITE (23,"(I4,4F16.8)") ISTEP,X0,DDX,dum1+rmag,dum2
                     theta1(i,j) = X0
!
                 end if
                 ! search theta for grid (R,Z)
                  do k=1,na
                     dum4(k) = (rf(theta(k),dum3) -gridr(i))/abs(gridr(i))
                     dum5(k) = (e*zf(theta(k),dum3) -gridz(j))/abs(gridz(j))
                  end do
              else
                     theta1(i,j) = 0.0
                     !WRITE (23,"(I4,4F16.8)") ISTEP,theta1(i,j),DDX,dum1,dum2
!
              end if
!

           end do
       end do
       call cdbbscoef2d(gridr,gridz,theta1,knotr,knotz,kr,kz,bscoeftheta) ! theta
        dx1 = 3/100.0
        do i= 1, 100
           do j= 1, 100
              tdum = -1.5 + i*dx1 + 5e-6
              rdum = -1.5 + j*dx1 + 5e-6
              dum = cdbbsval2d(knotx,knoty,kr,kz,bscoefxy2R,tdum,rdum,1,0)
              dum1 = cdbbsval2d(knotx,knoty,kr,kz,bscoefxy2R,tdum,rdum,0,1)
              dum2 = cdbbsval2d(knotx,knoty,kr,kz,bscoefxy2Z,tdum,rdum,1,0)
              dum3 = cdbbsval2d(knotx,knoty,kr,kz,bscoefxy2Z,tdum,rdum,0,1)
              !write(23,106) dum, dum1, dum2, dum3
           end do
         end do
         106 FORMAT (4ES25.16)
           ! rdum = 0.3988
            rdum = 0.3882
            tdum = 1.873476 
            a = 0.172117
            !tdum = 1.891255  
            !a = 0.175989
         do i = 1, na
             x0 = (rf(theta(i),rdum)*rmag - tdum)!/abs(tdum)
              b =    (e*zf(theta(i),rdum)*rmag - a)!/abs(a)
             !write(23,*) theta(i),x0,b
         end do 
       end subroutine check_xy
!
!       subroutine C_LKCT6 ! load phase space grids and data
!         implicit none
!         integer :: i, j, k, p, ir, jr, im, im1, l, kk, nx
!         integer, parameter :: nE=1, nL=1, nprt0=10, num=3000, nprt=8
!         integer, parameter :: nomode=3, nopmode=3, nog = 11
!         real(8), allocatable :: clkct6d_r(:,:,:,:,:,:),clkct6d_i(:,:,:,:,:,:)
!         real(8) :: dum, dum1, dum2, dum3, dum4, dum5, edum
!         integer :: ii(nprt), jj(nprt), k2(nprt)
!         real(8), allocatable :: ikct4d(:,:,:,:)
 
!         allocate(clkct6d_r(nprt0,nog,nog,nomode,nomode,nopmode))
!         allocate(clkct6d_i(nprt0,nog,nog,nomode,nomode,nopmode))
!         allocate(ikct4d(nog,nog,nomode,nomode))

!         open(89,file='grid_E.dat',status='unknown')
!         read(89,*)
!         open(90,file='mesh_xal_arc.dat',status='unknown')
!         read(90,*)
!         open(91,file='phasefreq.dat',status='unknown')
!         read(91,*)
!         open(92,file='Y6D.dat',status='unknown')
!         read(92,*)

!         do l=1, nprt0
!            ii(l) = 0
!            jj(l) = 0
!            k2(l) = 0
!            do ir =1, nog
!               do jr=1, nog
!                  do im=1, nomode
!                     do im1=1, nomode
!                       do p=1, nopmode

!                          clkct6d_r(l,ir,jr,im,im1,p) = 0.0
!                          clkct6d_i(l,ir,jr,im,im1,p) = 0.0


!                       enddo
!                     enddo
!                  enddo
!               enddo
!            enddo
!         enddo

!         do i =1, num
!            read(92,*) dum, dum1, dum2, dum3, dum4, kk, p, ir, jr, im, im1
!            clkct6d_r(kk, ir, jr, im, im1, p) = dum3
!            clkct6d_i(kk, ir, jr, im, im1, p) = dum4
!         enddo

!         do k = 1, nL
!            read(88,*) mubedum
!            do j = 1, nE
!               read(89,*) edum
!               nx = 0
!               do i = 1, nprt
!                  read(90,*) dum, dum1, dum2, dum3, dum4, dum5, kk
!                  if (abs(edum-dum1)<=1e-4.and.abs(dum2-mubedum)<=1e-4) then
!                     nx = nx +1 ! x grid changed with ions
!                     ii(kk) = nx  ! kth vs grids
!                     jj(kk) = j
!                     k2(kk) = k
!                  end if
!                enddo
!                nx0(k,j) = nx
!             enddo
!          enddo

!            do ir =1, nog
!               do jr=1, nog
!                  do im=1, nomode
!                     do im1=1, nomode
!                       do p=1, nopmode
!                          do l = 1, nprt0
!                             dum = clkct6d_r(l,ir,jr,im,im1,p)
!                             if (dum.eq.0) then
!                                ikct4d(ir,jr,im,im1) = 0.0
!                             else   
!                               do k = 1, nL
!                                  do j = 1, nE
!                                     allocate(fi(nx0(k,j)))
!                                     do nx=1, nx0(k,j)
!                                           if (ii(l)==nx.and.jj(l)=j.and.k2(l)=k) then
!                                              fi(nx) = dum
!                                           endif
!                                     enddo
!                                     do i = 1, nx-1
!                                        dum = (fi(i) + fi(i+1))/2.0
!                                        dum1 = x3d(i+1,j,k) - x3d(i,j,k)
!                                        ikct1(k,j) = ikct1(k,j) + dum*dum1
!                                     enddo      
!                                     deallocate(fi)
!                                  enddo
!                                  call simp(nE,dedum,ikct1(k,1:nE),sdum) 
!                                  ikct2(k) = sdum  
            
!                               enddo

!                                  call simp(nL,dldum,ikct2(1:nL),sdum)
!                                  ikct4d(ir,jr,im,im1) = sdum 
!
!                             endif
!                          enddo
!
!                       enddo
!                     enddo
!                  enddo
!               enddo
!            enddo
         


!       end subroutine C_LKCT6

!
       ! SUBROUTINE SIMP(N,H,FI,S)
!
! Subroutine for integration over f(x) with the Simpson rule.  FI:
! integrand f(x); H: interval; S: integral.  Copyright (c) Tao Pang 1997.
!
         !use omp_lib
        ! IMPLICIT NONE
        ! INTEGER :: N
        ! INTEGER :: I
        ! REAL(r8) :: H
        ! REAL(r8) :: S0,S1,S2
        ! REAL(r8) :: S
        ! REAL(r8), DIMENSION (N) :: FI
!
        ! S  = 0.0
        ! S0 = 0.0
        ! S1 = 0.0
        ! S2 = 0.0
        !!$omp parallel
        ! print *,'cpucore:',omp_get_num_procs(), omp_get_thread_num()
        !!$omp do firstprivate(N,FI) reduction(+:S0,S1,S2)

         !DO I = 2, N-1, 2
         !   S1 = S1+FI(I-1)
         !   S0 = S0+FI(I)
         !   S2 = S2+FI(I+1)
         !END DO
         !!$omp end do
         !!$omp end parallel
         !S = H*(S1+4.0*S0+S2)/3.0
!
! If N is even, add the last slice separately
!
         !IF (MOD(N,2).EQ.0) S = S &
         !     +H*(5.0*FI(N)+8.0*FI(N-1)-FI(N-2))/12.0
       !END SUBROUTINE SIMP
!
        !SUBROUTINE SIMPC(N,H,FI,S)
!
! Subroutine for integration over f(x) of complex with the Simpson rule.  FI:
! integrand f(x); H: interval; S: integral.  Copyright (c) Tao Pang 1997.
!
         !IMPLICIT NONE
         !INTEGER :: N
         !INTEGER :: I
         !REAL(r8) :: H
         !COMPLEX(8) :: S0,S1,S2
         !COMPLEX(8) :: S
         !COMPLEX(8), DIMENSION (N) :: FI
!
         !S  = (0.0,0.0)
         !S0 = (0.0,0.0)
         !S1 = (0.0,0.0)
         !S2 = (0.0,0.0)
         !DO I = 2, N-1, 2
         !   S1 = S1+FI(I-1)
         !   S0 = S0+FI(I)
         !   S2 = S2+FI(I+1)
         !END DO
         !S = H*(S1+4.0*S0+S2)/3.0
!
! If N is even, add the last slice separately
!
         !IF (MOD(N,2).EQ.0) S = S &
         !     +H*(5.0*FI(N)+8.0*FI(N-1)-FI(N-2))/12.0
       !END SUBROUTINE SIMPC

!
       SUBROUTINE SECANT (DL,X0,DX,ISTEP,DUM,I,J)
!
       ! Subroutine for the root of f(x)=0 with the secant method.
       ! Copyright (c) Tao Pang 1997.
!
          IMPLICIT NONE
          INTEGER, INTENT (INOUT) :: ISTEP
          INTEGER, INTENT (IN) :: I,J
          REAL(r8), INTENT (INOUT) :: X0,DX
          REAL(r8) :: X1,X2,D
          REAL(r8), INTENT (IN) :: DL,DUM
!
          ISTEP = 0
          X1 = X0+DX
          DO WHILE (ABS(DX).GT.DL)
             D  = FX(X1,DUM,I,J)-FX(X0,DUM,I,J)
             X2 = X1-FX(X1,DUM,I,J)*(X1-X0)/D
             X0 = X1
             X1 = X2
             DX = X1-X0
             ISTEP = ISTEP+1
          END DO
      END SUBROUTINE SECANT
!
      SUBROUTINE SECANT1 (DL,X0,DX,ISTEP,DUM,R,Z)
!
       ! Subroutine for the root of f(x)=0 with the secant method.
       ! Copyright (c) Tao Pang 1997.
!
          IMPLICIT NONE
          INTEGER, INTENT (INOUT) :: ISTEP
          REAL(r8), INTENT (INOUT) :: X0,DX
          REAL(r8) :: X1,X2,D
          REAL(r8), INTENT (IN) :: DL,DUM,R,Z
!
          ISTEP = 0
          X1 = X0+DX
          DO WHILE (ABS(DX).GT.DL)
             D  = FX1(X1,DUM,R,Z)-FX1(X0,DUM,R,Z)
             X2 = X1-FX1(X1,DUM,R,Z)*(X1-X0)/D
             X0 = X1
             X1 = X2
             DX = X1-X0
             ISTEP = ISTEP+1
             IF (ISTEP >= 100) THEN
                EXIT
             END IF
          END DO
      END SUBROUTINE SECANT1
!
      SUBROUTINE SECANT2 (DL,X0,DX,ISTEP,DUM,R)
!
       ! Subroutine for the root of f(x)=0 with the secant method.
       ! Copyright (c) Tao Pang 1997.
!
          IMPLICIT NONE
          INTEGER, INTENT (INOUT) :: ISTEP
          REAL(r8), INTENT (INOUT) :: X0,DX
          REAL(r8) :: X1,X2,D
          REAL(r8), INTENT (IN) :: DL,DUM,R
!
          ISTEP = 0
          X1 = X0+DX
          DO WHILE (ABS(DX).GT.DL)
             D  = FX2(X1,DUM,R)-FX2(X0,DUM,R)
             X2 = X1-FX2(X1,DUM,R)*(X1-X0)/D
             X0 = X1
             X1 = X2
             DX = X1-X0
             ISTEP = ISTEP+1
             IF (ISTEP >= 100) THEN
                EXIT
             END IF
          END DO
      END SUBROUTINE SECANT2
!
      SUBROUTINE SECANT3 (DL,X0,DX,ISTEP,DUM,LAM)
!
       ! Subroutine for the root of f(x)=0 with the secant method.
       ! Copyright (c) Tao Pang 1997.
!
          IMPLICIT NONE
          INTEGER, INTENT (INOUT) :: ISTEP
          REAL(r8), INTENT (INOUT) :: X0,DX
          REAL(r8) :: X1,X2,D
          REAL(r8), INTENT (IN) :: DL,DUM,LAM
!
          ISTEP = 0
          X1 = X0+DX
          DO WHILE (ABS(DX).GT.DL)
             D  = FX3(X1,DUM,LAM)-FX3(X0,DUM,LAM)
             X2 = X1-FX3(X1,DUM,LAM)*(X1-X0)/D
             X0 = X1
             X1 = X2
             DX = X1-X0
             ISTEP = ISTEP+1
             IF (ISTEP >= 100) THEN
                EXIT
             END IF
          END DO
      END SUBROUTINE SECANT3



!      FUNCTION FX(X,DUM,I,J) 
       FUNCTION FX(X,DUM,I,J)
          IMPLICIT NONE
          !common/aspect/e
          !REAL(r8) :: FX,e
          REAL(r8) :: FX
          INTEGER, INTENT (IN) :: I,J
          REAL(r8), INTENT (IN) :: X,DUM
          !REAL(r8), INTENT (IN) :: X
!
          !FX = EXP(X)*LOG(X)-X*X-DUM-I-J-e
         FX = (rf(X,DUM) -gridr(I))/abs(gridr(I))-&
             (e*zf(X,DUM) -gridz(J))/abs(gridz(J))
       END FUNCTION 
!
       FUNCTION FX1(X,DUM,R,Z)
          IMPLICIT NONE
          !common/aspect/e
          !REAL(r8) :: FX1,e
          REAL(r8) :: FX1
          REAL(r8), INTENT (IN) :: X,DUM,R,Z
!
!         FX1 = (rf(X,DUM) - R)/abs(R)-&
!               (e*zf(X,DUM) - Z)/abs(Z)
         !FX1 = rf(X,DUM) - R
         FX1 = e*zf(X,DUM) - Z
       END FUNCTION
!
       FUNCTION FX2(X,DUM,R)
          IMPLICIT NONE
          !common/aspect/e
          !REAL(r8) :: FX2,e
          REAL(r8) :: FX2
          REAL(r8), INTENT (IN) :: X,DUM,R

          FX2 = rf(DUM,X) - R

       END FUNCTION

!
       FUNCTION FX3(X,DUM,LAM)
          IMPLICIT NONE
          REAL(r8) :: FX3
          REAL(r8), INTENT (IN) :: X,DUM,LAM

          FX3 =1.0- b_nor(DUM,X)*LAM 

       END FUNCTION

!
      SUBROUTINE linspace(from, too, n,array)
        IMPLICIT NONE
        integer :: n,i
        real(r8) :: from, too
        real(r8),dimension(n) :: array
        REAL(r8) :: rangee

        n = size(array)
        !print *, n
        rangee= too - from
        if (n==0) return
        if (n==1) then
           array(1) = from
           return
        end if

        do i=1,n
           array(i) = from + rangee*real(i-1)/real(n-1)
        end do
      END SUBROUTINE linspace


!

      end module orbit_eq_mod
