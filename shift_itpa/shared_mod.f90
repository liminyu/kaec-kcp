module vars
  use inrtype, only: sp
  implicit none
  integer :: nprt,nomode,nopmode,ntor,midx
  integer numL,nE,nt,pa,pb, nrho,iname0
  integer, allocatable :: mode(:),pmode(:)
  real(sp) :: rhoh_bar,vh,rmaj,omeg0,ekev,zprt,bkg,prot,vA_bar
  real(sp) :: La,Lb,Ea,Eb,E0,Ed,Ec,L0,Ld,rr0,rrd,prot0
  real(sp) :: rhoh, ne0, v_A0
  real(sp) :: Ec32,E32,tE,tL,Cn
  real(sp), allocatable :: en(:), mube(:),ptch(:),ptcha(:),dxa(:),dE0(:)
  real(sp), allocatable :: pz0(:),dLa(:),dpz0(:),gridx2(:),gridpsi2(:)
  real(sp), allocatable :: R0(:), Z0(:), phi0(:),sgnb(:),XR0(:)
  real(sp), allocatable :: xx0(:), theta0(:),endum(:),&
                            fedum(:),fldum(:), lambdum(:)
  real(sp), allocatable :: tdum(:),Rdum(:),Zdum(:),&
                           phidum(:),psidum(:),thetadum(:),&
                           psirdum(:),vparaldum(:)
   real(sp), allocatable :: tdum1(:), Rdum1(:),Zdum1(:),&
                             phidum1(:),phi1dum1(:),&
                             psidum1(:),thetadum1(:),&
                             psirdum1(:)
   real(sp), allocatable :: f1(:,:),f2(:,:),f3(:,:),&
                             f4(:,:),f5(:,:),f6(:,:)
   real(sp), allocatable :: cscoefxpsi(:,:),cscoefxjr2(:,:),&
                            cscoefxnh(:,:)
   real(sp) :: pw
   
   ! pw: psi wall normaized to B(0)*a**2
   real(sp), allocatable :: w_theta(:), w_phi(:),&
                            w_theta1(:),w_phi1(:),&
                            tau_theta1(:),psir_bar(:)
   complex(sp), allocatable :: ylmp(:,:,:,:) ! ylmp(pmode,idxion,mmode,lfem)
   complex(sp), allocatable :: fylmp(:,:,:,:)
   complex(sp), allocatable :: fylmpp(:,:,:,:)
   integer, allocatable :: name0(:)
   integer ikct
   real(sp) res
end module vars

module vars_k
 implicit none
 integer, parameter, private :: r8 = selected_real_kind(12,100)
 complex(r8),allocatable :: aint_KCT4D(:,:,:,:)
 complex(r8),allocatable :: temp_aint_KCT4D(:,:,:,:)
 complex(r8),allocatable :: aint_KCT4Dp(:,:,:,:)
 complex(r8),allocatable :: temp_aint_KCT4Dp(:,:,:,:)
! integer,allocatable :: imdum(:),im1dum(:),&
!            ildum(:),il1dum(:),&
!            ipdum(:)
!  integer test1_number
!  real, allocatable :: local_a(:),total_a(:)
end module vars_k

module vars_e
 implicit none
 integer, parameter, private :: r8 = selected_real_kind(12,100)
 real(r8) :: e ! !aspect rate:e=a/R_0
 real(r8) :: etai,eta_i
 complex(r8) :: omg,lam_omg !omega_mode, resisitive eta_i, eigen_value
 complex(kind=4) :: rho
 real(r8) :: q0,q1,cq1,xini,xend !qfundat
 real(r8) :: delta ! step of raidal grids xgrid1
 real(r8) :: b0,b0_h ! beta_eff=P/B^2 at core 
                     ! with respect to bulk plasmas and fast ions
 real(r8), allocatable :: xg(:) ! xgrida
 real(r8), allocatable :: del(:) !shafranov shift
 real(r8), allocatable :: cden(:),ctem(:) ! dendat, tfundat
 real(r8) :: dalpha,dedge,dedge1 ! denf
 real(r8) :: gam ! adiabatic_index
 complex(r8),allocatable :: aa(:,:),bb(:,:),dd(:,:) ! matab
 real(r8), allocatable :: C1tt(:,:) ! inertia_tt
 real(r8), allocatable :: C1tr(:,:) ! inertia_tr
 real(r8), allocatable :: C1rr(:,:) ! inertia_rr
 real(r8), allocatable :: C2rr(:,:) ! curvature_rr
 real(r8), allocatable :: C2tr(:,:) ! curvature_tr
 real(r8), allocatable :: C2tt(:,:) ! curvature_tt
 real(r8), allocatable :: C3rt(:,:) ! current_rt
 real(r8), allocatable :: C3tt(:,:) ! current_tt
 real(r8), allocatable :: C4rt(:,:) ! pressure_rt
 real(r8), allocatable :: C4tt(:,:) ! pressure_tt
 real(r8), allocatable :: C4rr(:,:) ! pressure_rr
 real(r8), allocatable :: C4tr(:,:) ! pressure_tr
 integer nog ! nogr radial grid
end module vars_e

module shared_mod
  use vars_e
  use vars
  use inrtype
  use splines

  implicit none
  integer, parameter, private :: r8 = selected_real_kind(12,100)

   contains
!===================================================================================
      function gfun1(x,m)
      complex(r8) :: gfun1
      integer :: m
      real(r8) :: x, xi0,h0
        xi0=abs(0.1D0/ak1(x,m))/sqrt(tfun(x))
        if(xi0.gt.10.0D0)xi0=10.0D0
        h0=-0.5*xi0*exp(-xi0**2)*sqrt(3.14)
!        h0=-0.00D0
!        gfun=(0.2D0)**2*cmplx(1.0D0,h0)*tfun(x)
        gfun1=ak1(x,m)**2*cmplx(0.88D0,h0)*tfun(x)
!        gfun=(ak1(x,m)**2*cmplx(0.5D0,h0)+0.04*den(x))*tfun(x)
!resistive         gfun=cmplx(0.0D0,-1.0D0)*ak1(x,m)**2
!        gfun=cmplx(0.0D0,-1.0D0)*exp(-(x-1.0)**2/1.0)
!        gfun=cmplx(0.16,0.0D0)
        if(x.gt.1.0D0)gfun1=cmplx(0.0D0,0.0D0)
!        return
      end function gfun1
!  
!  
 function gfun(x,m)
      complex(r8) :: gfun,zio,zout,fzout,zin,fone
      integer :: m
      real(r8) :: x, xi0,h0,etax
        xi0=abs(0.3D0/ak1(x,m))/sqrt(tfun(x))
        etax=etai*den(x)/tfun(x)**1.5
        zio=xi0*sqrt(omg)*cmplx(1.0,etax)
        zin=xi0*sqrt(omg)*cmplx(0.0,etax)
!        if(xi0.gt.10.0D0)xi0=10.0D0
!        zio=cmplx(xi0,0.0)
        call zeta1(zio,zout)
        
        fone=cmplx(1.0,0.0)
        fzout=(fone+zin*zout)/(fone+zio*zout)
        if(m.eq.4.and.xi0.gt.0.97)then
        write(61,*)"z,fz=",zio,zout,fzout
        endif
!        h0=-0.5*xi0*exp(-xi0**2)*1.77245
         h0=-0.00D0
!        gfun=(0.4D0)**2*cmplx(1.0D0,-etai)*tfun(x)
!        gfun=ak1(x,m)**2*(cmplx(0.375D0,h0)+0.5*fzout)*tfun(x)
!        gfun=ak1(x,m)**2*(cmplx(0.88D0,h0)*tfun(x)+cmplx(0.0D0,-1.0D0)*eta/tfun(x)**1.5D0) 
!        gfun=(ak1(x,m)**2*cmplx(0.5D0,h0)+0.04*den(x))*tfun(x)
!resistive         gfun=cmplx(0.0D0,-1.0D0)*ak1(x,m)**2
          gfun=cmplx(1.0D0,0.0)*lam_omg
!        gfun=cmplx(0.0D0,-1.0D0)*exp(-(x-1.0)**2/1.0)
!        gfun=cmplx(0.16,0.0D0)
        if(x.gt.1.0D0)gfun=cmplx(0.0D0,0.0D0)
!        gfun1=gfun
!        return
 end function gfun
!

      subroutine zeta(z,zetaoz)
      complex(r8) :: z,zetaoz,dzetaz,term,fmult,terme,an1,bn1,zsquar,hold,temp1,temp2,ddzeta,dddzet
      real(r8) :: imagte,imagmu,imagse,imagsu,realte,realmu,realsu,realse,x,y,fn,erro
      complex :: zsquar1,zsquar3
      integer :: n
      erro=1.e-15
      zsquar=z*z
      zsquar1=-zsquar
      x=real(z)
      y=aimag(z)
      fn=real(zsquar)
      if (y.gt.0.) go to 99
      if (abs(fn).lt.174..and.abs(aimag(zsquar)).lt.5.e4) go to 98
      if (fn.gt.0.) go to 97
      write(16,1) z
1     format (" argument wp of subroutine zeta has too large a negative imaginary part, wp = ",1pe14.7," + ",e14.7," i")
97    hold=(0.,0.)
      go to 99
98    hold=(0.,1.77245385090551603)*cexp(zsquar1)
99    if (x*x+y*y.gt.16.) go to 200
      if (abs(y).ge.1.) go to 300
      realte=-2.*x
      imagte=-2.*y
      realmu=.5*(imagte*imagte-realte*realte)
      imagmu=-imagte*realte
      realsu=realte
      imagsu=imagte
      if (x.eq.0..and.y.eq.0.) go to 103
      fn=3.
100   realse=realte
      imagse=imagte
      realte=(realse*realmu-imagse*imagmu)/fn
      imagte=(realse*imagmu+imagse*realmu)/fn
      realse=realsu
      imagse=imagsu
      realsu=realsu+realte
      imagsu=imagsu+imagte
      fn=fn+2.
      if (abs(realse-realsu).gt.erro.or.abs(imagse-imagsu).gt.erro) go to 100
103   x=realsu
      fn=imagsu
      if (y.gt.0.)hold=(0.,1.77245385090551603)*cexp(-zsquar1)
      zetaoz=cmplx(x,fn)+hold
      go to 401
200   fn=5.
      dddzet=6.
      term=dddzet
      fmult=.5/zsquar
201   terme=term
      term=term*fmult*fn*(fn-1.)/(fn-3.)
      zetaoz=term/terme
      if (abs(real(zetaoz))+abs(aimag(zetaoz)).gt.1.) go to 250
      zetaoz=dddzet
      dddzet=dddzet+term
      fn=fn+2.
      zsquar3=zetaoz-dddzet
      if (cabs(zsquar3).gt.erro) go to 201
250   dddzet=dddzet/(zsquar*zsquar)
      if (y.gt.0.) go to 260
      fn=1.
      if (y.lt.0.) fn=2.
      dddzet=dddzet-4.*fn*hold*z*(2.*zsquar-3.)
260   ddzeta=-(4.+(zsquar-.5)*dddzet)/(z*(2.*zsquar-3.))
      dzetaz=(2.-z*ddzeta)/(2.*zsquar-1.)
      zetaoz=-(1.+.5*dzetaz)/z
      go to 401
300   if (y.lt.0.) z=conjg(z)
      terme=(1.,0.)
      term=(0.,0.)
      dzetaz=term
      fmult=terme
      n=0
      an1=z
      bn1=-z*z+.5
301   temp1=bn1*term+an1*terme
      temp2=bn1*fmult+an1*dzetaz
      zetaoz=temp1/temp2
      dzetaz=(zetaoz-term/fmult)/zetaoz
      if (abs(real(dzetaz)).lt.erro.and.abs(aimag(dzetaz)).lt.erro)go to 302
      bn1=bn1+2.
      n=n+1
      an1=-.5*float(n*(n+n-1))
      terme=term
      dzetaz=fmult
      term=temp1
      fmult=temp2
      if (n.lt.30) go to 301
302   if (y.ge.0.) go to 401
      zetaoz=conjg(zetaoz)+2.*hold
401   continue
      end subroutine zeta

      subroutine zeta1(z,zout)
      complex(r8) :: z,zout,z1,zout2,zout1
      complex :: zsquare
      z1=conjg(z)
      zsquare=z*z
      call zeta(z1,zout1)
      zout2=cexp(-zsquare)*2.0*sqrt(3.1415926535)*cmplx(0.0,1.0)
      zout=conjg(zout1)+zout2
      end subroutine zeta1
  
!
        !x:r/a,t:theta,e:a/R_0
        function ddel(x)!Delta^prime(r)
          implicit none
          real(r8) x,al
       
          real(r8) :: ddel
          if(x==0.0D0) x=1.0D-6
          if(e==0.0D0) then
          ddel=0.0D0
          else
          al=-2.0D0*qfun(x)**2*betap(x)/e
          ddel=(e*x+al)/4.0
          endif
 !      
        endfunction ddel
!
          
        function ddel1(x)
          implicit none
          real(r8) x,al
          real(r8)::ddel1
          if(x==0.0D0) x=1.0D-6
          if (e==0.0D0) then
          ddel1=0.0D0
          else
          al=-2*qfun(x)**2*betap(x)
          ddel1=(e**2*x+al)/4.0
          endif
        endfunction ddel1
!
      !  subroutine delx !Shafranov shift :Delta(r)
      !    implicit none
      !    real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
      !    integer i,j,k !j:index of theta(j*theta),m:maxium of j
      !    
      !    integer,parameter ::n=101,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
       
      !    real(r8) x,dt,s,t,ta,tb,dx,y
      !    real(r8) :: fi(2010)
      !    if (allocated(del)) then
      !      deallocate(del)
      !    endif
      !    allocate(del(n1))
          !!          real(r8) ::f_C1tt
!        dt=5.0D-4
      !    dx=delta/float(nq)
        !           dt=two_pi/real(n)
        !           ta=0.0
        !           tb=two_pi
        !
        !  do k=0,nq*nog
       !      x=dx*k+xini
           !          write(30,*)'x=',x
      !       if (x==0.0) then
      !          x=1.0D-6
      !       endif
      !       dt=(1.0D0-x)/float(n-1)
!              dt=x/float(n-1)
      !       do i=1,n
      !          t=(i-1)*dt+x
!                 t=(i-1)*dt
      !          fi(i)=ddel1(t)
      !       enddo
             !
      !       call simp(n,dt,fi,s)
                   !
      !       del(k)=-s
      !       write(30,*) x, del(k)
      !    enddo
      !  endsubroutine delx
!
        function delx1(x)!analytic Delta(x) Shafranov shift
          implicit none
          real(r8) :: x, dum, dum1
          real(r8) :: delx1
          dum = e**2*x**2/8+&
                b0/2*(0.2*q0**2*x**2/4+&
                      1.9*q0**2*x**2/2+&
                      0.064*q0*x**6/6+&
                      0.608*q0*x**4/4+&
                      0.00512*x**8/8+&
                      0.04864*x**6/6)
          !dum1 = e**2/8+b0/2*(0.2*q0**2/4+&
          !           1.9*q0**2/2+0.064*q0/6+&
          !           0.608*q0/4+0.00512/8+&
          !           0.04864/6)
          !delx1= dum - dum1
          delx1 = dum
        endfunction delx1      


        function eps1(x)!epsilon
          implicit none
          real(r8) x
          real(r8) :: eps1
          eps1=e*x
        endfunction eps1
!
        function eta(x)
          implicit none
          real(r8) x
!          common/aspect/e
          real(r8) ::eta
          eta=(ddel(x)+eps1(x))/2.0D0
        endfunction eta
!
        function theta_s(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::theta_s
          theta_s=t+(eps1(x)*sin(t)+&
                     eps1(x)**2/4*sin(2.0*t))/&
                     (1.0-eps1(x)**2/2)
        endfunction theta_s        
!
        function thetaspx(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::thetaspx
          thetaspx=e*((1.0-eps1(x)**2/2.0)*(sin(t)+&
                     eps1(x)/2.0*sin(2.0*t))+eps1(x)**2*&
                      (sin(t)+eps1(x)/4.0*sin(2.0*t)))/&
                      (1.0-eps1(x)**2/2)**2
        endfunction thetaspx
!
        function thetaspt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::thetaspt
          thetaspt=rf(t,x)
        endfunction thetaspt
!
        function thetaspxt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::thetaspxt
          thetaspxt=e*((1.0-eps1(x)**2/2.0)*(cos(t)+&
                     eps1(x)*cos(2.0*t))+eps1(x)**2*&
                      (cos(t)+eps1(x)/2.0*cos(2.0*t)))/&
                      (1.0-eps1(x)**2/2)**2
        endfunction thetaspxt
!

        function rf(t,x) !!R/R0 of straight field line coordinates
          implicit none
          integer,parameter ::nq=4,nog1=800,n1=nq*nog1
          integer n
          real(r8) t,x,dx
          real(r8) ::rf
!
          if(x==0.0D0) x=1.0D-6
          !dx=delta/float(nq)
!        x=0.5D0
          !n=x/dx+0.01D0
!
!          rf=1.0+e*x*cos(t)-delx1(x)+eps1(x)*eta(x)*(cos(2.0*t)-1.0)
          rf=1.0+e*x*cos(t)-(e*x)**2/8.0+eps1(x)*eta(x)*(cos(2.0*t)-1.0)
!           rf = 1.0 + eps1(x)*cos(theta_s(t,x))
        endfunction rf
!
        function zf(t,x) !!Z/R0 of straight field line coordinates
          implicit none
          real(r8) t,x
          real(r8) ::zf
!
          zf=eps1(x)*sin(t)+eps1(x)*eta(x)*sin(2.0*t)
!          zf = eps1(x)*sin(theta_s(t,x))
        endfunction zf

!
        function rf2(t,x)
          implicit none
          real(r8) t,x
!          common/aspect/e
          real(r8) ::rf2
          rf2=rf(t,x)**2
        endfunction rf2
!
!
        function rfpr(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::rfpr
!
          !if (e==0.0D0) then
          !rfpr=cos(t)
          !else
          !rfpr=cos(t)-ddel(x)+(eta_r(x)+0.625D0*e*x-0.25D0*qfun(x)*&
          !     x*betap(x)*dqfun(x)/e-0.125D0*qfun(x)**2&
          !     *x*betapp(x)/e)*(cos(2.0*t)-1.0D0)
          !endif
!
          rfpr=e*(cos(t)-ddel(x)+2.0*eta(x)*(cos(2.0*t)-1.0)) ! 20240831
!          rfpr = e*cos(theta_s(t,x))-eps1(x)*sin(theta_s(t,x))*thetaspx(t,x)
        endfunction rfpr
!
       function rfpx(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::rfpx
!
!          rfpx = e*cos(theta_s(t,x))-eps1(x)*sin(theta_s(t,x))*thetaspx(t,x)
          rfpx=e*(cos(t)-ddel(x)+2.0*eta(x)*(cos(2.0*t)-1.0)) ! 20240831
        endfunction rfpx
!
        function rfpt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::rfpt
!
!          rfpt=-x*sin(t)-2.0*x*eta_r(x)*sin(2.0*t)
          rfpt=-eps1(x)*sin(t)-2.0*eps1(x)*eta(x)*sin(2.0*t) ! 20240831
!           rfpt =-eps1(x)*sin(theta_s(t,x))*thetaspt(t,x)
        endfunction rfpt
!
!

        function jf(t,x)!jacobian
          implicit none
          real(r8) x,t
          real(r8) ::jf
          jf=x*rf2(t,x)
          !jf=x*(1.0+2.0*eps1(x)*cos(t))
        endfunction jf
!
        function jfpx(t,x)
          implicit none
          real(r8) x,t
          real(r8) ::jfpx
          jfpx=rf2(t,x)+2.0*x*rf(t,x)*rfpr(t,x)
          !jfpx=1.0+4.0*eps1(x)*cos(t)
        endfunction jfpx

!
        function jfpt(t,x)!jacobian
          implicit none
          real(r8) x,t
          real(r8) ::jfpt
          jfpt=2.0*x*rf(t,x)*rfpt(t,x)
          !jfpt=-x*(2.0*eps1(x)*sin(t))
        endfunction jfpt

!
        function grr(t,x)!metrix coefficients
          implicit none
          real(r8) t,x
          real(r8) ::grr
          !grr=1.0
          grr=1.0+2.0*ddel(x)*cos(t)
        endfunction grr
!
        function grt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::grt
          if(x.le.1.0D-3) x=1.0D-3
          if(x.ge.1.0D0) x=1.0-1e-06

          !grt=-thetaspx(t,x)/thetaspt(t,x)
           grt=-(eps1(x)+2.0*ddel(x))*sin(t)/x
        endfunction grt
!
        function gzz(t,x)
          implicit none
          real(r8) t,x
          real(r8) :: gzz
!          common/aspect/e
          if(x.le.1.0D-3) x=1.0D-3
          if(x.ge.1.0D0) x=1.0-1e-06

          gzz=1.0/rf(t,x)/rf(t,x)
          !gzz=1.0-2.0*eps1(x)*cos(t)
        endfunction gzz

!
        function gtt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::gtt
          !gtt=1.0/x**2/rf2(t,x)+grt(t,x)**2
          gtt=(1.0-2.0*(eps1(x)+ddel(x))*cos(t))/x**2
        endfunction gtt
!
        function f_C1tt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C1tt
          f_C1tt=jf(t,x)*rf2(t,x)*gtt(t,x)
        endfunction f_C1tt
!
        function f_C1tr(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C1tr
          f_C1tr=jf(t,x)*rf2(t,x)*grt(t,x)
        endfunction f_C1tr
!
        function f_C1rr(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C1rr
          f_C1rr=jf(t,x)*rf2(t,x)*grr(t,x)
        endfunction f_C1rr
!
        subroutine C1_tt!get fourier coupling coefficients
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C1tt)) then
            deallocate(C1tt)
          endif
          allocate(C1tt(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=dx*k+xini
          if (x==0.0) then
           x=1.0D-6
          endif
          do j=0,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C1tt(t,x)*cos(float(j)*t)
             enddo
!
             call simp(n,dt,fi,s)
             if (j==0)then
                C1tt(j,k)=s/pi
             else
                C1tt(j,k)=2*s/pi
             endif
          enddo
       enddo
        endsubroutine C1_tt
!
!
        subroutine C1_tr
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
!
          if (allocated(C1tr)) then
            deallocate(C1tr)
          endif
          allocate(C1tr(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=dx*k+xini
          if (x==0.0)then
            x=1.0D-6
          endif
          do j=1,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C1tr(t,x)*sin(float(j)*t)
             enddo
             call simp(n,dt,fi,s)
             C1tr(j,k)=2*s/pi
          
          enddo
       enddo
     endsubroutine
!
     subroutine C1_rr
       implicit none
       real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C1rr)) then
            deallocate(C1rr)
          endif
          allocate(C1rr(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=dx*k+xini
          if (x==0.0)then
           x=1.0D-6
          endif
       do j=0,m
          do i=1,n
             t=(i-1)*dt
             fi(i)=f_C1rr(t,x)*cos(float(j)*t)
          enddo
!
          call simp(n,dt,fi,s)
          if (j==0)then
             C1rr(j,k)=s/pi
          else
             C1rr(j,k)=2.0*s/pi
          endif
       enddo
       enddo
     endsubroutine C1_rr
!
!
     function gamma1(x,k,m,j)
       implicit none
       integer m,j,k,jp,i,n
       integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
       real(r8) x,x1
       real(r8) ::gamma1
!
       x1=x-xini
       n=x1*float(nq)/delta+0.01D0

! 
       jp=abs(j)
       if (jp<=m1)then
!
          if (jp==0) then
             gamma1=C1tt(jp,n)*float(k*m)
             return
          elseif (k-m+jp==0)then
             gamma1=C1tt(jp,n)*float(k*m)/2.0
          elseif (k-m-jp==0)then
             gamma1=C1tt(jp,n)*float(k*m)/2.0
          endif
       else
          gamma1=0.0
       endif
!
     endfunction gamma1
!
!
     function delta1(x,k,m,j)
       implicit none
      
       real(r8) ::delta1
!
       integer m,j,k,jp,i,n
       integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
       real(r8) x,x1
!
       x1=x-xini
       n=x1*float(nq)/delta+0.01D0      
!
       jp=abs(j)
       if(jp<=m1)then
!
!
          if (jp==0) then
             delta1=0.0
             return
          elseif (k-m+jp==0)then
             delta1=float(k)*C1tr(jp,n)/2.0
          elseif(k-m-jp==0)then
             delta1=-float(k)*C1tr(jp,n)/2.0
          endif
       else
          delta1=0.0
       endif
!
     endfunction delta1
!
     function theta1(x,k,m,j)
       implicit none
       
       real(r8) ::theta1
!
       integer m,j,k,jp,i,n
       integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
       real(r8) x,x1
!
       x1=x-xini
       n=x1*float(nq)/delta+0.01D0      
!
       jp=abs(j)
       if (jp<=m1) then
!
!
          if (jp==0) then
             theta1=0.0
             return
          elseif (k-m+jp==0)then
             theta1=float(m)*C1tr(jp,n)/2.0
          elseif(k-m-jp==0)then
             theta1=-float(m)*C1tr(jp,n)/2.0
          endif
       else
          theta1=0.0
       endif
!
     endfunction theta1
!
!
     function lammbda1(x,k,m,j)
       implicit none
       
       real(r8) ::lammbda1
!
       integer m,j,k,jp,i,n
       
       integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
       real(r8) x,x1
!
       x1=x-xini
       n=x1*float(nq)/delta+0.01D0      
!
       jp=abs(j)
       if (jp<=m1) then
!
!
          if (jp==0) then
             lammbda1=C1rr(jp,n)
          elseif (k-m+jp==0)then
             lammbda1=C1rr(jp,n)/2.0
          elseif (k-m-jp==0)then
             lammbda1=C1rr(jp,n)/2.0
          endif
       else
          lammbda1=0.0
       endif
!
     endfunction lammbda1
!
      subroutine check4
        implicit none
        real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
        integer,parameter :: n=101,na=361
        real(r8) :: R(na),Z(na)
        integer :: i,j
        real(r8) :: dt,dx,dum
        real(r8) :: x(n),t(na)
        dt=two_pi/(na-1)
        dx=(1.0)/(n-1)

        do i=1,n
           x(i)=(i-1)*dx
           if (x(i)==0) then
              x(i)=5.1e-03
           endif
        enddo

        do i=1,na
           t(i)=(i-1)*dt
        enddo
       do i=1,n
          do j=1,na
          if (mod(i,10)==0) then
          write(40,*) rf(t(j),x(i)),zf(t(j),x(i))
          endif
          enddo
       enddo
!
       do i=1,n
          do j=1,na
          if (mod(j,12)==0) then
          write(40,*) rf(t(j),x(i)),zf(t(j),x(i))
          endif
          enddo
       enddo

!
        do i=1,n
           write(41,441) x(i),qfun(x(i)),dqfun(x(i))
        enddo

441   format(1E12.6,1x,1E12.6,1x,1E12.6)

       do i=1,n
           !dum = inrcsval(gridx2,cscoefxpsi,x(i),0)/pw
           !write(42,441) dum,den(dum),denp(dum)
           !write(42,441) x(i),den(x(i)),denp(x(i))
           dum = inrcsval(gridx2,cscoefxpsi,x(i),0)/pw
           write(42,441) dum,den_h(dum),denp_h(dum)
        enddo
        close(40)
        close(41)
        close(42)
     end subroutine check4
!
        function F(x)
          implicit none
          real(r8) x
          real(r8) ::F
          F=jr2(x)
        endfunction F
!
        subroutine C2_rr
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
!
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C2rr)) then
            deallocate(C2rr)
          endif
          allocate(C2rr(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=dx*k+xini
          if (x==0.0)then
            x=1.0D-6
          endif
!
          do j=0,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=grr(t,x)*cos(float(j)*t)
             enddo
!
             call simp(n,dt,fi,s)
             if (j==0)then
                C2rr(j,k)=s/pi
             else
                C2rr(j,k)=2.0*s/pi
             endif
          enddo
       enddo
!
        endsubroutine C2_rr
!
        subroutine C2_tr
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
!
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C2tr)) then
            deallocate(C2tr)
          endif
          allocate(C2tr(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nog*nq
          x=dx*k+xini
          if(x==0.0)then
           x=1.0D-6
          endif
!
          do j=1,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=grt(t,x)*sin(float(j)*t)
             enddo
!
             call simp(n,dt,fi,s)
             C2tr(j,k)=2*s/pi
        
          enddo
       enddo
        endsubroutine C2_tr
!
!
        subroutine C2_tt
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
!
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C2tt)) then
            deallocate(C2tt)
          endif
          allocate(C2tt(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=dx*k+xini
          if(x==0.0)then
           x=1.0D-6
          endif
!
          do j=0,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=gtt(t,x)*cos(float(j)*t)
             enddo
             call simp(n,dt,fi,s)
             if (j==0)then
                C2tt(j,k)=s/pi
             else
                C2tt(j,k)=2.0*s/pi
             endif
          enddo
       enddo
        endsubroutine C2_tt
!
        function gamma2(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::gamma2
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0        
! 
          jp=abs(j)
          if (jp<=m1) then
! 
!
             if (j==0) then
                gamma2=C2tt(jp,n)*float(k*m)
             elseif (k-m+jp==0)then
                gamma2=C2tt(jp,n)*float(k*m)/2.0
             elseif (k-m-jp==0)then
                gamma2=C2tt(jp,n)*float(k*m)/2.0
             endif
          else
             gamma2=0.0
          endif
!
        endfunction gamma2
!
!
        function delta2(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::delta2
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0

          jp=abs(j)
!
          if (jp<=m1) then
!
!
             if (j==0) then
                delta2=0.0
             elseif (k-m+jp==0)then
                delta2=C2tr(jp,n)*float(k)/2.0
             elseif (k-m-jp==0)then
                delta2=-C2tr(jp,n)*float(k)/2.0
             else
                delta2=0.0
             endif
          else
             delta2=0.0
          endif
!
        endfunction delta2

!
        function theta2(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::theta2
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
!
          jp=abs(j)
          if (jp<=m1) then
!
!
             if (j==0) then
                theta2=0.0
                return
             elseif (k-m+jp==0)then
                theta2=C2tr(jp,n)*float(m)/2.0
             elseif (k-m-jp==0)then
                theta2=-C2tr(jp,n)*float(m)/2.0
             endif
          else
             theta2=0.0
          endif

!
        endfunction theta2
!
!
        function lammbda2(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::lammbda2
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
!
          jp=abs(j)
          if(jp<=m1) then
!
!
             if (j==0) then
                lammbda2=C2rr(jp,n)
             elseif (k-m+jp==0)then
                lammbda2=C2rr(jp,n)/2.0
             elseif (k-m-jp==0)then
                lammbda2=C2rr(jp,n)/2.0
             endif
          else
             lammbda2=0.0
          endif
!
        endfunction lammbda2
!

!
!
        function chip(x)
          implicit none
          real(r8) x
          real(r8) ::chip
          chip=jr2(x)/qfun(x)
        endfunction chip
!
!
        function jr2(x)
          implicit none
          real(r8) x
          real(r8) :: jr2
          jr2=x
        endfunction jr2
!
        function jr20(t,x)
          implicit none
          real(r8) x,t
          real(r8) :: jr20
          jr20=jf(t,x)/rf2(t,x)
        endfunction jr20

!
        function jr2p(x)
          implicit none
          real(r8) x
          real(r8) ::jr2p
          jr2p=1.0
        endfunction jr2p
!
        function jr20p(t,x)!jacobian
          implicit none
          real(r8) x,t
          real(r8) ::jr20p
          jr20p=-2.0/rf(t,x)**3*jf(t,x)*&
                e*rfpr(t,x)+jfpx(t,x)/rf2(t,x)

        endfunction jr20p

!
        function grrpx(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::grrpx
          !grrpx=0.0
          grrpx=e*cos(t)/2.0
        endfunction grrpx
!
        function grrpt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::grrpt
          !grrpt=0.0
          grrpt = -2.0*ddel(x)*cos(t)
        endfunction grrpt

!
        function gtrpt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::gtrpt
          if(x.le.1.0D-3) x=1.0D-3
          if(x.ge.1.0D0) x=1.0-1e-06
!
          !gtrpt=-thetaspxt(t,x)/thetaspt(t,x)-&
          !       rfpt(t,x)/thetaspt(t,x)*grt(t,x)
          gtrpt=-(eps1(x)+2.0*ddel(x))*cos(t)/x
        endfunction gtrpt
!
        function gtrpptx(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::gtrpptx
!
          
          gtrpptx=-1.5*e*cos(t)/x+(eps1(x)+2.0*ddel(x))*cos(t)/x**2
        endfunction gtrpptx
!
        function chipp(x)
          implicit none
          real(r8) x
          real(r8) ::chipp
          chipp=jr2p(x)/qfun(x)-jr2(x)*dqfun(x)/qfun(x)**2
        endfunction chipp
!
        function jr2pp(x)
          implicit none
          real(r8)x
          real(r8)::jr2pp
          jr2pp=0.0
        endfunction jr2pp
!
        function grrppx(t,x)
        implicit none
        real(r8) t,x
        real(r8)::grrppx
        grrppx=0.0
        endfunction
!
        function chippp(x)
          implicit none
          real(r8) x
          real(r8)::chippp
          chippp=jr2pp(x)/qfun(x)-jr2p(x)*dqfun(x)/qfun(x)**2- &
               jr2p(x)*dqfun(x)/qfun(x)**2-jr2(x)*ddqfun(x)/qfun(x)**2+ &
               2*jr2(x)*dqfun(x)**2/qfun(x)**3
        endfunction chippp
!
        
        function f_C3rt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C3rt
          f_C3rt=(jr2p(x)*grr(t,x)*chip(x)+jr2(x)*grrpx(t,x)*chip(x)+ &
           jr2(x)*grr(t,x)*chipp(x))/jr2(x)+gtrpt(t,x)*chip(x)
		  
!
        endfunction f_C3rt
!
        function f_C3tt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C3tt
!
          f_C3tt=-jr2p(x)*(jr2p(x)*grr(t,x)*chip(x)+jr2(x)*grrpx(t,x)*chip(x)+ &
               jr2(x)*grr(t,x)*chipp(x))/jr2(x)**2+(jr2pp(x)*grr(t,x)*chip(x)+ &
               jr2p(x)*grrpx(t,x)*chip(x)+jr2p(x)*grr(t,x)*chipp(x)+ &
               jr2p(x)*grrpx(t,x)*chip(x)+jr2(x)*grrppx(t,x)*chip(x)+ &
               jr2(x)*grrpx(t,x)*chipp(x)+jr2p(x)*grr(t,x)*chipp(x)+ &
               jr2(x)*grrpx(t,x)*chipp(x)+jr2(x)*grr(t,x)*chippp(x))/jr2(x)+ &
               gtrpt(t,x)*chipp(x)+gtrpptx(t,x)*chip(x)
		  
!
        endfunction f_C3tt
!
        subroutine C3_rt
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C3rt)) then
            deallocate(C3rt)
          endif
          allocate(C3rt(0:m,0:n1))
!
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
         do k=0,nq*nog
            x=k*dx+xini
            if (x==0)then
             x=1.0D-6
            endif
!
          do j=0,m
!
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C3rt(t,x)*cos(float(j)*t)
             enddo
!
             call simp(n,dt,fi,s)
             if (j==0)then
                C3rt(j,k)=s/pi
             else
                C3rt(j,k)=2.0*s/pi
             endif
          enddo
       enddo
!
        endsubroutine C3_rt
!
        subroutine C3_tt
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C3tt)) then
            deallocate(C3tt)
          endif
          allocate(C3tt(0:m,0:n1))
!
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
        do k=0,nog*nq
           x=k*dx+xini
           if(x==0.0) then
           x=1.0D-6
           endif
!
          do j=0,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C3tt(t,x)*cos(float(j)*t)
             enddo
             call simp(n,dt,fi,s)
             if (j==0)then
                C3tt(j,k)=s/pi
             else
                C3tt(j,k)=2.0*s/pi
             endif
          enddo
       enddo
!
        endsubroutine C3_tt
!
        function gamma3(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::gamma3
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
!
          jp=abs(j)
          if(jp<=m1) then
!
!
             if (j==0) then
                gamma3=float(k)*C3tt(jp,n)
             elseif (k-m+jp==0)then
                gamma3=C3tt(jp,n)*float(k)/2.0
             elseif (k-m-jp==0)then
                gamma3=C3tt(jp,n)*float(k)/2.0
             endif
          else
             gamma3=0.0
          endif
!
        endfunction gamma3
!
!
!
!
        function theta3(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::theta3
!
          x1=x-xini
          n=x1*float(nq)/delta+0.01D0
!
          jp=abs(j)
          if(jp<=m1) then
!
!
             if (j==0) then
                 theta3=C3rt(jp,n)
                return 
             elseif (k-m+jp==0)then
                theta3=C3rt(jp,n)/2.0
             elseif (k-m-jp==0)then
                theta3=C3rt(jp,n)/2.0
             endif
          else
             theta3=0.0
          endif
!
        endfunction theta3
!
!
        function beta(x) ! normalized pressure, P/B^2
          implicit none
          real(r8) x, dum
          real(r8) :: beta
          dum = inrcsval(gridx2,cscoefxpsi,x,0) ! psi
          dum = dum/pw ! for q=q(r)
          !dum = dum/1.0 ! for q=q(psi)
          !beta=b0*(1-x**2)!beta(0)=0.2%
!           beta=0.0
!          beta=b0*(1.0 - 0.95*dum - 0.05*dum**2)
          beta=b0*(1.0D0-x**2)
          ! beta=b0*(1-x**2)**3
!           beta=b0
!        beta=b0*(6.0141D0*x**5-19.613D0*x**4+22.183D0*x**3-9.6863D0*x**2+0.11273D0*x+0.99488D0)
!          beta=b0*(0.97996D0-10.79691D0*x**2+63.52315D0*x**4-213.40D0*x**6+&
!                   419.77D0*x**8-4.742163D0*x**10+284.694D0*x**12-70.22167D0*x**14)
!        beta=b0*(2.452D0*exp(-5.396D0*x)-1.46D0*exp(-9.566D0*x))
!        beta=b0*(19.95D0*x**6-66.51D0*x**5+78.57D0*x**4-&
!                 37.29D0*x**3+5.035D0*x**2-0.7637D0*x+1.005D0)
!        beta=b0*(-9.425D0*x**5+18.19D0*x**4-8.598D0*x**3-0.8569D0*x**2-0.3393D0*x+0.9966D0)
!       beta=b0*(1.0D0-1.0D0*x**2)
!       beta=b0*(1.0 - 0.95*x**2- 0.05*x**4) ! ITPA n=6
        endfunction beta
!
        function betap(x)
          implicit none
          real(r8) x, dum, dum1
          real(r8) :: betap
          dum = inrcsval(gridx2,cscoefxpsi,x,0) ! psi
          dum1 = inrcsval(gridx2,cscoefxpsi,x,1) ! dpsi/dx
          dum = dum/pw ! for q=q(r)
          dum1 = dum1/pw
          !dum = dum/1.0 ! for q=q(psi)
          !betap=-2*b0*x
!           betap=0.0
!          betap=b0*(-0.95-0.1*dum)*dum1
!       betap=b0*(30.0705D0*x**4-78.452D0*x**3+66.549D0*x**2-19.3726D0*x+0.11273D0)          
!       betap=b0*(-21.59382D0*x+254.0926D0*x**3-1280.4D0*x**5+3358.16D0*x**7-47.42163D0*x**9+3416.328D0*x**11-983.10338D0*x**13)
!        betap=b0*(-13.231D0*exp(-5.396D0*x)+13.96636D0*exp(-9.566D0*x))
!        betap=b0*(119.7D0*x**5-332.55D0*x**4+314.28D0*x**3-111.87D0*x**2+10.07D0*x-0.7637D0)
!        betap=b0*(-47.125D0*x**4+72.76D0*x**3-25.794D0*x**2-1.7138D0*x-0.3393D0)
        betap=b0*(-2.0D0*x)
!        betap=b0*( - 2.0*0.95*x- 4.0*0.05*x**3) ! ITPA n=6
!         betap=-6.0D0*b0*(1-x**2)**2*x
        endfunction betap
!
        function betapp(x)
          implicit none
          real(r8) x, dum, dum1, dum2
          real(r8) :: betapp
!          dum = inrcsval(gridx2,cscoefxpsi,x,0) ! psi
!          dum1 = inrcsval(gridx2,cscoefxpsi,x,1) ! dpsi/dx
!          dum2 = inrcsval(gridx2,cscoefxpsi,x,2) ! dpsi^2/dx^2
!          dum = dum/pw ! for q=q(r)
!          dum1 = dum1/pw
!          dum2 = dum2/pw
          !dum = dum/1.0 ! for q=q(psi)
          betapp=-2*b0
!           betap=0.0
!          betapp=b0*((-0.1)*dum1+&
!                 (-0.95-0.1*dum)*dum2)
!           betapp=b0*( - 2.0*0.95- 3.0*4.0*0.05*x**2) ! ITPA n=6
!         betapp=24.0D0*b0*(1-x**2)*x**2-6.0D0*b0*(1-x**2)**2
        endfunction betapp
!
        function kappa_r(t,x)
          implicit none
          real(r8) x,t
          real(r8) ::kappa_r
!          kappa_r=-rfpr(t,x)/rf(t,x)-e**2*x/(qfun(x)*rf(t,x))**2
!           kappa_r=-rfpr(t,x)*(1-e*cos(t))-e**2*x*(1-2*e*cos(t))/(qfun(x))**2
!           kappa_r=-e*rfpr(t,x)/rf(t,x)-e**2*x/(qfun(x)*rf(t,x))**2
          kappa_r=b_norpx(t,x)/b_nor(t,x)
        endfunction kappa_r
!
        function kappa_t(t,x)
          implicit none 
          real(r8) t,x
          real(r8) ::kappa_t
!          kappa_t=-rfpt(t,x)/rf(t,x)
!           kappa_t=-rfpt(t,x)*(1-e*cos(t))
!            kappa_t=-e*rfpt(t,x)/rf(t,x)
          kappa_t = b_norpt(t,x)/b_nor(t,x)
        endfunction kappa_t
!
      function b_nor(t,x) ! B/B(0)
       implicit none
       real(r8) ::  x, t
       real(r8) ::b_nor

       b_nor = sqrt(1.0+(e*jr2(x)/qfun(x))**2*grr(t,x))/rf(t,x)

      end function b_nor

!
!
      function b_norpt(t,x)
       implicit none
       real(r8) ::  x, t
       real(r8) ::b_norpt
       real(r8) :: rr,rr1,gr,gr1
       real(r8) :: dum,dum1,qq
       rr=rf(t,x)
       rr1=rfpt(t,x)
       gr=grr(t,x)
       gr1=grrpt(t,x)
       dum=e*jr2(x)
       qq=qfun(x)
       dum1=sqrt(1.0+dum**2/qq**2*gr)
       b_norpt = -1.0/rr/rr*rr1*dum1+&
                 0.5/rr/dum1*dum**2/qq**2*gr1
      end function b_norpt

!
      function b_norpx(t,x)
       implicit none
       real(r8) ::  x, t
       real(r8) ::b_norpx
       real(r8) :: rr,rr1,rr2,dum,dum1,qq2
       real(r8) :: dum2, dum3
       real(r8) :: gr,gr1
       rr=rf(t,x)
       rr1=rfpr(t,x)
       rr2=rf(t,x)**2
       dum=(e*jr2(x))**2
       qq2=qfun(x)**2
       gr = grr(t,x)
       gr1 = grrpx(t,x)
       dum1=sqrt(1.0+dum/qq2*gr)
       dum2 = 2.0*jr2(x)/qfun(x)*&
              (jr2p(x)/qfun(x)-jr2(x)/qq2*dqfun(x))
       dum3 = jr2(x)**2/qq2
       b_norpx = -1.0/rr2*rr1*dum1+&
                  0.5/rr/dum1*e**2*(gr*dum2+gr1*dum3)
      end function b_norpx

!
      ! readme_x1 Eq.(36) for b_norpt
      ! b_norpt_bar_rfpt = b_norpt/rfpt
      function b_norpt_bar_rfpt(t,x)
       implicit none
       real(r8) ::  x, t
       real(r8) ::b_norpt_bar_rfpt
       real(r8) :: rr,dum,dum1,qq,gr
       rr=rf(t,x)
       dum=e*jr2(x)
       qq=qfun(x)
       gr = grr(t,x)
       dum1=sqrt(1.0+dum**2/qq**2*gr)
       ! case of grrpt = 0 
       b_norpt_bar_rfpt = -1.0/rr/rr*dum1
        
      end function b_norpt_bar_rfpt
!
      function b_star(t,x,v4) ! b**(theta,r,v||)
       implicit none
       real(r8) :: t,x,v4,b2,g
       real(r8) :: g1,g2,g3,g4
       real(r8) :: j1,j2,j3,j4,j5,j6
       real(r8) :: rr,rr1,rr2,qq,qq1
       real(r8) :: b_star,bb
       bb=b_nor(t,x)
       b2=b_nor(t,x)**2
       g=gzz(t,x)
       g1=grrpx(t,x)
       g2=grr(t,x)
       g3=gtrpt(t,x)
       g4=grt(t,x)
       j1=jf(t,x)
       j2=jfpx(t,x)
       j3=jfpt(t,x)
       j4=jr2(x)
       j5=jr2p(x)
       rr=rf(t,x)
       rr1=rfpr(t,x)
       rr2=rfpt(t,x)
       qq=qfun(x)
       qq1=dqfun(x)
       b_star = bb + v4*g*j4/b2/qq*&
                (g1+j2/j1*g2+g3+j3/j1*g4)-&
                v4/b2/qq/rr*j4**g*(g2*rr1+g4*rr2)+&
                v4/b2*(j5/qq-j4*qq1/qq**2)*g*g2
      end function b_star
!
        function realfuny(mube,x,t,theta,phi1,wtheta,nq,j,m,p)
          implicit none
          real(r8) :: mube, x, t,theta,phi1,wtheta
          integer j, m, p, nq
          real(r8) :: realfuny, dum1, dum2
          real(r8) :: xa,xb
          real(r8) :: eps
          !dum1 = mube*b_nor(theta,x)+2.0*(1.0-mube*b_nor(theta,x))
           dum1 = 2.0/b_nor(theta,x)-mube
          dum2 = 1.0/b_nor(theta,x)/jf(theta,x)
          realfuny = dum1*dum2*(-kappa_t(theta,x)*femp(j,m,x)*&
                cos(float(nq)*phi1-float(m)*theta-float(p)*wtheta*t)+&
                float(m)*kappa_r(theta,x)*fem(j,m,x)*&
                sin(float(nq)*phi1-float(m)*theta-float(p)*wtheta*t))

!
         endfunction realfuny

!
        function imagfuny(mube,x,t,theta,phi1,wtheta,nq,j,m,p)
          implicit none
          integer n
          real(r8) :: mube, x, t,theta,phi1,wtheta
          integer j, m, p, nq
          real(r8) :: imagfuny, dum1, dum2
          real(r8) :: xa,xb
          real(r8) :: eps
          !dum1 = mube*b_nor(theta,x)+2.0*(1.0-mube*b_nor(theta,x)) !b~=R/R0
           dum1 = 2.0/b_nor(theta,x)-mube
          dum2 = 1.0/b_nor(theta,x)/jf(theta,x)
          imagfuny = dum1*dum2*(-kappa_t(theta,x)*femp(j,m,x)*&
                 sin(float(nq)*phi1-float(m)*theta-float(p)*wtheta*t)-&
                 float(m)*kappa_r(theta,x)*fem(j,m,x)*&
                cos(float(nq)*phi1-float(m)*theta-float(p)*wtheta*t))
!
         endfunction imagfuny

!
        function f_C4rt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C4rt
!
          f_C4rt=rf2(t,x)*kappa_t(t,x)*(2*gam*beta(x)*kappa_r(t,x)-betap(x))
        endfunction f_C4rt
!
        function f_C4tt(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C4tt
          !
          f_C4tt=(2*gam*beta(x)*kappa_r(t,x)-betap(x))*kappa_r(t,x)*rf2(t,x)
        endfunction f_C4tt
!
        function f_C4rr(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C4rr
!
          f_C4rr=2*gam*beta(x)*kappa_t(t,x)*kappa_t(t,x)*rf2(t,x)
        endfunction f_C4rr
!
        function f_C4tr(t,x)
          implicit none
          real(r8) t,x
          real(r8) ::f_C4tr
          !
          f_C4tr=2*gam*beta(x)*kappa_t(t,x)*kappa_r(t,x)*rf2(t,x)
        endfunction f_C4tr
!
        subroutine C4_rt
          implicit none
          real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          real(r8) x,dt,s,t,dx
          real(r8) :: fi(n)
          if (allocated(C4rt)) then
            deallocate(C4rt)
          endif
          allocate(C4rt(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=k*dx+xini
          if(x==0.0)then
          x=1.0D-6
          endif
          do j=1,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C4rt(t,x)*sin(float(j)*t)
             enddo
             call simp(n,dt,fi,s)
             C4rt(j,k)=2.0*s/pi
          enddo
       enddo
     endsubroutine C4_rt
!
        subroutine C4_tt
          implicit none
          real,parameter ::pi=3.14159265358979,two_pi=6.28318530717959d0
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C4tt)) then
            deallocate(C4tt)
          endif
          allocate(C4tt(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
        do k=0,nq*nog
           x=k*dx+xini
           if(x==0.0)then
           x=1.0D-6
           endif
          do j=0,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C4tt(t,x)*cos(j*t)
             enddo
!
             call simp(n,dt,fi,s)
             if (j==0)then
                C4tt(j,k)=s/pi
             else
                C4tt(j,k)=2.0*s/pi
             endif
          enddo
       enddo
          !
        endsubroutine C4_tt
!
!
        subroutine C4_rr
          implicit none
          real,parameter :: two_pi=6.28318530717959d0,pi=3.14159265358979
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          real(r8) x,dt,s,t,ta,tb,dx
          real(r8) :: fi(n)
          if (allocated(C4rr)) then
            deallocate(C4rr)
          endif
          allocate(C4rr(0:m,0:n1))
          dt=pi/float(n-1)
          dx=delta/float(nq)
!
       do k=0,nq*nog
          x=k*dx+xini
          if(x==0.0)then
          x=1.0D-6
          endif
          do j=0,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C4rr(t,x)*cos(float(j)*t)
             enddo
!
             call simp(n,dt,fi,s)
             if (j==0)then
                C4rr(j,k)=s/pi
             else
                C4rr(j,k)=2.0*s/pi
             endif
          enddo
       enddo
        endsubroutine C4_rr
!
!
        subroutine C4_tr
          implicit none
          real,parameter :: two_pi=6.28318530717959d0,pi=3.14159265358979
          integer,parameter ::n=15,m=3,nq=4,nog1=800,n1=nq*nog1 !n:integral steps
          integer i,j,k !j:index of theta(j*theta),m:maxium of j
          real(r8) x,dt,s,t,dx
          real(r8) :: fi(n)
          if (allocated(C4tr)) then
            deallocate(C4tr)
          endif
          allocate(C4tr(0:m,0:n1))
!
          dt=pi/float(n-1)
          dx=delta/float(nq)
       do k=0,nq*nog
          x=k*dx+xini
          if(x==0.0)then 
           x=1.0D-6
          endif
          do j=1,m
             do i=1,n
                t=(i-1)*dt
                fi(i)=f_C4tr(t,x)*sin(float(j)*t)
             enddo
!
             call simp(n,dt,fi,s)
             C4tr(j,k)=2.0*s/pi
          enddo
       enddo
        endsubroutine C4_tr
!
!
        function gamma4(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::gamma4
!
      
          jp=abs(j)
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
          if(jp<=m1)then
!
! 
             if (j==0) then
                gamma4=float(k*m)*C4tt(jp,n)
                return
             elseif (k-m+jp==0)then
                gamma4=C4tt(jp,n)*float(k*m)/2.0
             elseif (k-m-jp==0)then
                gamma4=C4tt(jp,n)*float(k*m)/2.0
             endif
          else
             gamma4=0.0
          endif
!
        endfunction gamma4
!
!
        function delta4(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::delta4
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
          jp=abs(j)
          if(jp<=m1)then
!
!
             if (j==0) then
                delta4=0.0
             elseif (k-m+jp==0)then
                delta4=C4tr(jp,n)*float(k)/2.0
             elseif (k-m-jp==0)then
                delta4=-C4tr(jp,n)*float(k)/2.0
             endif
          else
             delta4=0.0
          endif
!
        endfunction delta4
!
!
        function theta4(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::theta4
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
          jp=abs(j)
!
          if(jp<=m1)then
!
! 
             if (j==0) then
                theta4=0.0
             elseif (k-m+jp==0)then
                theta4=C4rt(jp,n)*float(m)/2.0
             elseif (k-m-jp==0)then
                theta4=-C4rt(jp,n)*float(m)/2.0
             endif
          else
             theta4=0.0
          endif
        endfunction theta4
!
!
        function lammbda4(x,k,m,j)
          implicit none
          integer,parameter ::m1=3,nq=4,nog1=800,n1=nq*nog1 !m1:max number of fourier coefficients
          integer m,j,k,jp,n
          real(r8) x,x1
          real(r8) ::lammbda4
!
          x1=x-xini
          n=float(nq)*x1/delta+0.01D0
          jp=abs(j)
!
          if(jp<=m1)then
!
! 

             if (j==0) then
                lammbda4=C4rr(jp,n)
                return
             elseif (k-m+jp==0)then
                lammbda4=C4rr(jp,n)/2.0
             elseif (k-m-jp==0)then
                lammbda4=C4rr(jp,n)/2.0
             endif
          else
             lammbda4=0.0
          endif
        endfunction lammbda4
!
      function GetFileN(ifileUnit)
      implicit none
      integer :: GetFileN
      integer,intent(IN) :: ifileUnit
      integer :: ios
      character(len=1) :: cDummy
      GetFileN = 0
      Rewind(ifileUnit)
      do
         read(ifileUnit,*,iostat=ios) cDummy
         if(ios/=0) exit
         GetFileN=GetFileN+1
      enddo
    end function

!
        subroutine set_couple_coefficients
          call C1_tt
          call C1_tr
          call C1_rr
          call C2_rr
          call C2_tr
          call C2_tt
          call C3_rt
          call C3_tt
          call C4_rt
          call C4_tt
          call C4_rr
          call C4_tr
         ! call KCT4D
        endsubroutine set_couple_coefficients
!
!  

      subroutine irange(ng,nofmode,j,i1,i2)
      integer ng,nofmode,j,i1,i2
      i1=((j-1)/nofmode)*nofmode+1-3*nofmode
      i2=((j-1)/nofmode)*nofmode+4*nofmode
      if(i1.lt.1)i1=1
      if(i2.gt.ng)i2=ng
!      return
      end subroutine irange

      subroutine xrange(j,i,xa,xb,nx)
      implicit none
      integer n,nq,i,j,nx,ny
      parameter(n=800,nq=4)
      real(r8) :: xa,xb
      real(r8) :: eps
      eps=1.0D-6
      eps=xini+eps
      !if(j.eq.1)then
      if(j.eq.0) then
      xa=eps
      xb=xg(j+2)
      else if(j.eq.(nog-1))then
      xa=xg(j-2)
      xb=xg(nog)
      else
      xa=xg(j-2)
      xb=xg(j+2)
     ! if(i.lt.j)xb=xg(i+2)
     ! if(i.gt.j)xa=xg(i-2)
     ! if(xa.lt.eps)xa=eps
      endif

      if(i.lt.j)xb=xg(i+2)
      if(i.gt.j)xa=xg(i-2)
      if(xa.lt.eps)xa=eps

      ny=(xb-xa)/delta+0.01D0
      nx=ny*nq
      return
      end subroutine xrange

!
      subroutine xrange1(j,i,xa,xb,nx)
      implicit none
      integer n,nq,i,j,nx,ny
      parameter(n=800,nq=4)
      real(r8) :: xa,xb
      real(r8) :: eps
      eps=1.0D-6
      if(j.eq.0) then
      xa=eps
      xb=xg(j+2)
      else if(j.eq.(nog-1))then
      xa=xg(j-2)
      xb=xg(nog)
      else
      xa=xg(j-2)
      xb=xg(j+2)
      endif

      if(i.lt.j)xa=xg(i-2)
      if(i.gt.j)xb=xg(i+2)
      if(xa.lt.eps)xa=eps
      if(i.eq.(nog-1))xb=xg(nog)


      ny=(xb-xa)/delta+0.01D0
      nx=ny*nq
      return
      end subroutine xrange1

!

      subroutine setgrid
      implicit none
      integer n,i
      parameter(n=800)
      real(r8) :: h
!  
      if (allocated(xg)) then
         deallocate(xg)
      endif
      allocate(xg(-3:n+3))
!     h=xend/float(nog)
      h= (xend-xini)/float(nog)
      delta=h
      do i=-3,nog+3
      xg(i)=h*float(i)+xini
      enddo
!      return
      end subroutine setgrid

!
      function fem(i,m,x)
      implicit none
      integer n,i,m
      parameter(n=800)
      real(r8) :: fem,x
      if(i.eq.1)then
      fem=felst(m,x)
      else if(i.eq.(nog-1))then
      fem=felend(m,x)
      else
      fem=fel(i,x)
      endif
!      return
      end function fem


      function femp(i,m,x)
      implicit none
      integer n,i,m
      parameter(n=800)
      real(r8) :: x,femp
      if(i.eq.1)then
      femp=felstp(m,x)
      else if(i.eq.(nog-1))then
      femp=felendp(m,x)
      else
      femp=felp(i,x)
      endif
!      return
      end function femp


      function fempp(i,m,x)
      implicit none
      integer n,i,m
      parameter(n=800)
      real(r8) :: x,fempp
      if(i.eq.1)then
      fempp=felstpp(m,x)
      else if(i.eq.(nog-1))then
      fempp=felendpp(m,x)
      else
      fempp=felpp(i,x)
      endif
!      return
      end function fempp

!  
      function fel(i,x)
      implicit none
      integer n
      parameter(n=800)
      integer i,i0,im1,im2,ip1,ip2
      real(r8) :: fel,x,xbar
      i0=i
      im1=i-1
      im2=i-2
      ip1=i+1
      ip2=i+2

      if(x.gt.xg(im2).and.x.le.xg(im1))then
      xbar=(x-xg(im2))/(xg(im1)-xg(im2))
      fel=xbar**3
!  
      else if(x.gt.xg(im1).and.x.le.xg(i0))then
!  
      xbar=(x-xg(i0))/(xg(i0)-xg(im1))
      fel=4.0D0-xbar**2*(6.0D0+3.0D0*xbar)
!  
      else if(x.gt.xg(i0).and.x.le.xg(ip1))then
      xbar=(x-xg(i0))/(xg(ip1)-xg(i0))
      fel=4.0D0-xbar**2*(6.0D0-3.0D0*xbar)

      else if(x.gt.xg(ip1).and.x.le.xg(ip2))then
      xbar=(x-xg(ip2))/(xg(ip2)-xg(ip1))
      fel=-xbar**3
      else
      fel=0.0D0
      endif
!
!      return
      end function fel
!  
      function felend(m,x)
!  
      implicit none
      integer n,m
      parameter(n=800)
      real(r8) :: felend,x,cm1,cn,deltax,augden

      augden=1.0D0
      deltax=delta*(denp(augden)/den(augden)+1.0D0)
      cm1=(3.0D0+deltax)/(3.0D0-deltax)
      cn = -0.5D0*deltax/(3.0D0-deltax)
      felend=cm1*fel(nog-1,x)+cn*fel(nog,x)-fel(nog+1,x)
!      return
      end function felend

      function felendp(m,x)
      implicit none
      integer n,m
      parameter(n=800)
      real(r8) :: felendp,x,cm1,cn,deltax,augden
      augden=1.0D0
      deltax=delta*(denp(augden)/den(augden)+1.0D0)
      cm1=(3.0D0+deltax)/(3.0D0-deltax)
      cn = -0.5D0*deltax/(3.0D0-deltax)
      felendp=cm1*felp(nog-1,x)+cn*felp(nog,x)-felp(nog+1,x)
!      return
      end function felendp

      function felendpp(m,x)
      implicit none
      integer n,m
      parameter(n=800)
      real(r8) :: felendpp,x,cm1,cn,deltax,augden
      augden=1.0D0
      deltax=delta*(denp(augden)/den(augden)+1.0D0)
      cm1=(3.0D0+deltax)/(3.0D0-deltax)
      cn = -0.5D0*deltax/(3.0D0-deltax)
      felendpp=cm1*felpp(nog-1,x)+cn*felpp(nog,x)-felpp(nog+1,x)
!      return
      end function felendpp
!  
      function felst(m,x)

      implicit none
      integer m
      real(r8) :: felst,x

      if(m.ne.1)then
      felst=fel(1,x)+fel(-1,x)-0.5D0*fel(0,x)
      else
      felst=fel(1,x)-fel(-1,x)
      endif
!      return
      end function felst

      function felstp(m,x)
      implicit none
      integer m
      real(r8) :: felstp,x
      if(m.ne.1)then
      felstp=felp(1,x)+felp(-1,x)-0.5D0*felp(0,x)
      else
      felstp=felp(1,x)-felp(-1,x)
      endif
!      return
      end function felstp

      function felstpp(m,x)
      implicit none
      integer m
      real(r8) :: felstpp,x
      if(m.ne.1)then
      felstpp=felpp(1,x)+felpp(-1,x)-0.5D0*felpp(0,x)
      else
      felstpp=felpp(1,x)-felpp(-1,x)
      endif
!      return
      end function felstpp

      function felp(i,x)
      implicit none
      integer n
      parameter(n=800)
      integer i,i0,im1,im2,ip1,ip2
      real(r8) :: felp,x,xbar
      i0=i
      im1=i-1
      im2=i-2
      ip1=i+1
      ip2=i+2

      if(x.gt.xg(im2).and.x.le.xg(im1))then
      xbar=(x-xg(im2))/(xg(im1)-xg(im2))
      felp=3.0D0*xbar**2/(xg(im1)-xg(im2))
!  
      else if(x.gt.xg(im1).and.x.le.xg(i0))then
!  
      xbar=(x-xg(i0))/(xg(i0)-xg(im1))
      felp=-(12.0D0*xbar+9.0D0*xbar**2)/(xg(i0)-xg(im1))
!  
      else if(x.gt.xg(i0).and.x.le.xg(ip1))then
      xbar=(x-xg(i0))/(xg(ip1)-xg(i0))
      felp=-(12.0D0*xbar-9.0D0*xbar**2)/(xg(ip1)-xg(i0))

      else if(x.gt.xg(ip1).and.x.le.xg(ip2))then
      xbar=(x-xg(ip2))/(xg(ip2)-xg(ip1))
      felp=-3.0D0*xbar**2/(xg(ip2)-xg(ip1))
      else
      felp=0.0D0
      endif
!  
!      return 
      end function felp

      function felpp(i,x)
      implicit none
      integer n
      parameter(n=800)
      integer i,i0,im1,im2,ip1,ip2
      real(r8) :: felpp,x,xbar
      i0=i
      im1=i-1
      im2=i-2
      ip1=i+1
      ip2=i+2

      if(x.gt.xg(im2).and.x.le.xg(im1))then
      xbar=(x-xg(im2))/(xg(im1)-xg(im2))
      felpp=6.0D0*xbar/(xg(im1)-xg(im2))**2
!  
      else if(x.gt.xg(im1).and.x.le.xg(i0))then
!  
      xbar=(x-xg(i0))/(xg(i0)-xg(im1))
      felpp=-(12.0D0+18.0D0*xbar)/(xg(i0)-xg(im1))**2
!  
      else if(x.gt.xg(i0).and.x.le.xg(ip1))then
      xbar=(x-xg(i0))/(xg(ip1)-xg(i0))
      felpp=-(12.0D0-18.0D0*xbar)/(xg(ip1)-xg(i0))**2

      else if(x.gt.xg(ip1).and.x.le.xg(ip2))then
      xbar=(x-xg(ip2))/(xg(ip2)-xg(ip1))
      felpp=-6.0D0*xbar/(xg(ip2)-xg(ip1))**2
      else
      felpp=0.0D0
      endif
!  
!      return
      end function felpp


       function ak2(x,m) 
       implicit none
       integer nq,m
       real(r8) :: ak2,x,akp
       akp=float(ntor)-float(m)/qfun(x)
       ak2=akp**2
!       return
       end function ak2
!
       function akpp(x,m)
         implicit none
         integer nq,m
         real(r8) :: akpp,x
         akpp=float(m)*dqfun(x)/qfun(x)**2
!       return
       end function akpp
!
       function akfemp(i,m,x)
         integer i,m
         real(r8) x
         real(r8) ::akfemp
         akfemp=akpp(x,m)*fem(i,m,x)+femp(i,m,x)*ak1(x,m)
       endfunction akfemp
       
!
       function akfem(i,m,x)
         implicit none
         integer i,m
         real(r8) x
         real(r8) ::akfem
         akfem=ak1(x,m)*fem(i,m,x)
       endfunction akfem
!
       function ak1(x,m) 
       implicit none
       integer nq,m
       real(r8) :: ak1,x,akp
       ak1=float(ntor)-float(m)/qfun(x)
!       return
       end function ak1

       function ak2p(x,m)
       implicit none
       integer nq,m
       real(r8) :: ak2p,akpp,x,akp
       akp=float(ntor)-float(m)/qfun(x)
       akpp=float(m)*dqfun(x)/qfun(x)**2
       ak2p=2.0D0*akp*akpp
!       return
       end function ak2p
!

       function den(x)
       implicit none
       integer i,nfit
       parameter(nfit=5)
       real(r8) :: x,xx,den,x1
       den=cden(1)
       x1=x
       if(x.gt.1.0D0)x1=1.0D0
       xx=x1**2
       do i=1,nfit-1
       den=den+cden(i+1)*xx**i
       enddo
       !den=den-dedge*exp((x1-1.0)/dalpha)
       den=1.0D0
       !den=(1.00001-0.8*x**1.6)**1.8
       end function den

       function den_h(psi)
       implicit none
       real(r8) :: psi,den_h
       real(r8) :: c0,c1,c2,c3
       real(r8) :: dum,dum1
       c0 = 0.49123
       c1 = 0.298228
       c2 = 0.198739
       c3 = 0.521298
       if(psi.gt.1.0D0)psi=1.0D0
       dum = psi ! psi_n
       dum1 = tanh((sqrt(dum)-c0)/c2)
       den_h= c3*exp(-c2*dum1/c1)
       end function den_h

       function den1(x)
       implicit none
       integer i,nfit
       parameter(nfit=5) 
       real(r8) :: x,xx,den,den1,x1
       den1=cden(1)
       x1=x
       if(x.gt.1.0D0)x1=1.0D0
       xx=x1**2
       do i=1,nfit-1
       den1=den1+cden(i+1)*xx**i
       enddo
       den1=den1-dedge1*exp((x1-1.0)/dalpha)
!       den1=1.0D0
       end function den1

       function denp(x)
       implicit none
       integer i,nfit
       parameter(nfit=5)
       real(r8) :: x,xx,denp
       xx=x**2
       do i=1,nfit-1
       denp=denp+(2*i)*cden(i+1)*x**(2*i-1)
       enddo
       !denp=denp-dedge*exp((x-1.0)/dalpha)/dalpha
       if(x.gt.1.0D0)denp=0.0D0
        denp=0.0D0    
        !denp=-2.304*(1.00001-0.8*x**1.6)**0.8*x**0.6
       end function denp
!  
       function denp_h(psi)
       implicit none
       real(r8) :: psi,denp_h
       real(r8) :: c0,c1,c2,c3
       real(r8) :: dum,dum1
       c0 = 0.49123
       c1 = 0.298228
       c2 = 0.198739
       c3 = 0.521298
       if(psi.gt.1.0D0)psi=1.0D0
       if(psi.eq.0.0D0)psi=1e-6
       dum = psi ! psi_n
       dum1 = tanh((sqrt(dum)-c0)/c2)
       denp_h = -c3/2/c1/sqrt(dum)*exp(-c2*dum1/c1)*&
              (1.0 - dum1**2)
       !denp_h = inrcsval(gridpsi2,cscoefxnh,dum,1) ! psi  
       end function denp_h  
!
       function tfun(x)
       implicit none
         integer i,nfit
         parameter(nfit=5)
         real(r8) :: x,xx,tfun
        tfun=ctem(1)
        xx=x**2
        if(x.gt.1.0D0)xx=1.0D0
        do i=1,nfit-1
        tfun=tfun+ctem(i+1)*xx**i
        enddo
        tfun=1.0D0        
        end function tfun
!  

       function qfun(x)
       implicit none
         real(r8) :: x,qfun,dq,cq2,dum,cq3,psi
         integer isp
         isp = 0 ! 0 for analytic, 1 for interp
         dq=-9.0D0
         cq2=4.0D0
         cq3=4.2D0
         psi=(cq2-cq3+q0)/(dq+cq2-2*(cq3-q0))
         !if (isp==1) then
         !  dum = inrcsval(gridx2,cscoefxpsi,x,0) ! psi    
         !  dum = dum/pw ! psi_n
         !  dq=q1-q0
         !  cq2=1.0D0-cq1
         !  qfun=q0+dq*(cq1*dum+cq2*dum**2)
         !else
         !  dum=x**2
           !qfun = q0 + 0.55D0*x**2
         !   qfun = 1.6667+0.65695*dum+1.3388*dum**2-&
         !          0.868*dum**3+9.0268e-2*dum**4+&
         !          0.2001*dum**5-8.4828e-2*dum**6
         ! endif
          qfun = q0 +0.16*x**2 ! ITPA n=6 TAE
         !qfun=q0+x**2*(cq3-q0+(cq2-cq3+q0)*(1-psi)*(x**2-1)/(x**2-psi))
         if(x.gt.1.0D0)qfun=q0+dq*(cq1+cq2)

!	return
       end function qfun
!  
       function dqfun(x)
       implicit none
         real(r8) :: x,dqfun,dq,cq2,dum,dum1,cq3,psi
         integer isp
         isp = 0 ! 0 for anlytic, 1 for interp
         dq=-9.0D0
          cq2=4.0D0
          cq3=4.2D0
          psi=(cq2-cq3+q0)/(dq+cq2-2*(cq3-q0))
!  
         !if (isp==1) then
         !dum = inrcsval(gridx2,cscoefxpsi,x,0) ! psi
         !dum = dum/pw ! psi_n
         !dum1 = inrcsval(gridx2,cscoefxpsi,x,1) ! dpsi/dx    
         !dum1 = dum1/pw ! dpsi_n/dx
         !dq=q1-q0
         !cq2=1.0D0-cq1
         !dqfun=dq*(1.0D0*cq1+2.0D0*cq2*dum)*dum1
         !else
         !dum = x
         !dqfun=2.7029e-4+1.29*dum+0.52801*dum**2+&
         !      1.1417*dum**3+15.515*dum**4-&
         !      33.74*dum**5+24.784*dum**6-&
         !      6.3498*dum**7
         !dqfun=1.1D0*x
         !endif
          dqfun=0.32D0*x ! ITPA n=6 TAE
         ! dqfun=-2*x*(q0-cq3+((psi-1)*(cq2-cq3+q0)*(x**2-1))/(x**2-psi))&
         !      -x**2*((2*x*(psi-1)*(cq2-cq3+q0))/(x**2-psi)&
         !      -(2*x*(psi-1)*(x**2-1)*(cq2-cq3+q0))/(x**2-psi)**2)
         if(x.gt.1.0D0)dqfun=0.0D0

!	return
       end function dqfun
!
      function ddqfun(x)
        implicit none
        real(r8) :: x,ddqfun,dq,cq2,dum,dum1,dum2,cq3,psi
        integer isp
        isp = 0 ! 0 for anlytic, 1 for interp
        dq=-9.0D0
        cq2=4.0D0
        cq3=4.2D0
        psi=(cq2-cq3+q0)/(dq+cq2-2*(cq3-q0))
!  
        !if (isp==1)then
        !dum = inrcsval(gridx2,cscoefxpsi,x,0) ! psi
        !dum = dum/pw ! psi_n
        !dum1 = inrcsval(gridx2,cscoefxpsi,x,1) ! dpsi/dx  
        !dum1 = dum1/pw ! dpsi_n/dx
        !dum2 = inrcsval(gridx2,cscoefxpsi,x,2) ! dpsi^2/dx^2  
        !dum2 = dum2/pw ! dpsi_n^2/dx^2
        !dq=q1-q0
        !cq2=1.0D0-cq1
        !ddqfun=dq*(2.0D0*cq2)*dum1**2+&
        !       dq*(1.0D0*cq1+2.0D0*cq2*dum)*dum2
        ! else
        ! dum = x
        ! ddqfun = 1.3292-0.50144*dum+23.101*dum**2-&
        !          53.792*dum**3+200.22*dum**4-&
        !          525.38*dum**5+662.01*dum**6-&
        !          394.65*dum**7+91.049*dum**8
        !ddqfun=1.1D0
        !endif
         ddqfun=0.32D0
        !ddqfun=2*cq3-2*q0&
        !         -4*x*((2*x*(psi-1)*(cq2-cq3+q0))/(x**2-psi)&
        !        -(2*x*(psi-1)*(x**2-1)*(cq2-cq3+q0))/(x**2-psi)**2)&
        !        -x**2*((2*(psi-1)*(cq2-cq3+q0))/(x**2-psi)&
        !        -(2*(psi-1)*(x**2-1)*(cq2-cq3+q0))/(x**2-psi)**2&
        !        -(8*x**2*(psi-1)*(cq2-cq3+q0))/(x**2-psi)**2&
        !        +(8*x**2*(psi-1)*(x**2-1)*(cq2-cq3+q0))/(x**2-psi)**3)&
        !        -(2*(psi-1)*(x**2-1)*(cq2-cq3+q0))/(x**2-psi)
        if(x.gt.1.0D0)ddqfun=0.0D0

!	return
       end function ddqfun

!
       SUBROUTINE SIMP(N,H,FI,S)
!
 ! Subroutine for integration over f(x) with the Simpson rule.  FI:
! integrand f(x); H: interval; S: integral.  Copyright (c) Tao Pang 1997.
!
         IMPLICIT NONE
         INTEGER :: N
         INTEGER :: I
         REAL(r8) :: H
         REAL(r8) :: S0,S1,S2
         REAL(r8) :: S
         REAL(r8), DIMENSION (N) :: FI
!
         S  = 0.0
         S0 = 0.0
         S1 = 0.0
         S2 = 0.0
         DO I = 2, N-1, 2
            S1 = S1+FI(I-1)
            S0 = S0+FI(I)
            S2 = S2+FI(I+1)
         END DO
         S = H*(S1+4.0*S0+S2)/3.0
!
! If N is even, add the last slice separately
!
         IF (MOD(N,2).EQ.0) S = S &
              +H*(5.0*FI(N)+8.0*FI(N-1)-FI(N-2))/12.0
       END SUBROUTINE SIMP


!
       subroutine profiledat
         integer nfit,i,j
         parameter(nfit=5)
         integer nE1,numL1
         real(r8) dum,dedum,dldum
        if (allocated(cden)) then
             deallocate(cden)
         endif
         if (allocated(ctem)) then
             deallocate(ctem)
         endif
         allocate(cden(nfit))
         allocate(ctem(nfit))
         q0=1.35D0
         q1=4.6D0
         cq1=0.5D0
         !e=0.125D0
!          
         prot0 = 1.0 ! the mass in proton unit 1: H 2: D 3: T (BuP)
         prot = 2.0  ! the mass in proton unit 1: H 2: D 3: T (EP)
         zprt=1.0    ! the charge in proton unit
         !ekev=100    ! particle energy in kev for normalization
                     ! or temperature of maxwellian distribution
!
         !ne0 = 1.08e+13 !cm^-3 EAST
         ! ne0 = 1.0e+13 ! n=3 TAE, GMEC,Liu2024
         ne0 = 2.0e+13 ! n=6 TAE, A. Könies et al., Nuclear Fusion 58, 126027 (2018)
         !ne0 = 1.2e+13 ! HL-2A
!
         eta_i = 5.0e-5
         read(2100,*)q0,q1,cq1,xini,xend,dalpha,dedge,&
                   dedge1,e,gam,b0,b0_h,etai,lam_omg
         print *, 'q0=',q0,'q1=',q1,'cq1=',cq1,'xini=',xini,&
                  'xend=',xend,'dalpha=',dalpha,'dedge=',dedge,&
                  'dedge1=',dedge1,'e=',e,'gam=',gam,&
                  'b0=',b0,'b0_h=',b0_h,'etai=',etai,'lam=',lam_omg
!        load r0 and delta r (a) for density profile of fast ions, nrho
!        for X grid of sub check_mesh2
         read(2100,*) rr0,rrd, nrho
         print*, 'rr0=',rr0,'rrd=',rrd,'nrho=',nrho
!        load min and max energy (E_i0), grid of energy, injection energy, delta
!        E, energy width and
!        critical energy (E_i0) for slowing down distribution
!        nE for energy grid of sub check_mesh2
         read(2100,*) Ea,Eb,nE,E0,Ed,Ec,ekev
         print*,'Ea=',Ea,'Eb=',Eb,'nE=',nE,&
                'E0=',E0,'Ed=',Ed,'Ec=',Ec,'ekev=',ekev
!        load min and max Lammbda, its grid, L0 and delta L for distribution
         read(2100,*) La,Lb,numL,L0,Ld
         print*,'La=',La,'Lb=',Lb,'numL=',numL,&
                'L0=',L0,'Ld=',Ld
!        load min and max p-index of bounce harmonics
        ! read(2100,*) pa,pb,nopmode,nomode
         read(2100,*) pa,pb,nopmode,ntor ! ntor: toroidal mode number
         print*,'ntor=',ntor
!        load finite element grid
         read(2100,*) nog
 
        do i=1,nfit
        read(2100,*)cden(i)
        enddo
!
        do i=1,nfit
        read(2100,*)ctem(i)
        enddo
        close(2100)

        read(3100,*) lam_omg
        close(3100)

         omg = sqrt(lam_omg)
         nE1 = 101
         numL1 = 101
         if (allocated(endum)) then
            deallocate(endum)
         endif
         if (allocated(fedum)) then
            deallocate(fedum)
         endif
         if (allocated(fldum)) then
            deallocate(fldum)
         endif
         if (allocated(lambdum)) then
            deallocate(lambdum)
         endif
         allocate(endum(nE1),fedum(nE1))
         allocate(fldum(numL1),lambdum(numL1))
         dedum = (10*Eb - Ea)/float(nE1-1)
         dldum = (Lb - La)/float(numL1-1)
         Ec32=(Ec/ekev)**1.5
         do j=1,numL1
            lambdum(j) = La + (j-1)*dldum
            do i=1,nE1
               endum(i)=Ea + (i-1)*dedum
               E32=(endum(i)/ekev)**1.5+Ec32
               tE=(endum(i)-E0)/Ed
               tL=(lambdum(j)-L0)/Ld
               tL=-1.0*tL*tL
               fedum(i)=4.0*sqrt(PI_D)*sqrt(endum(i)/ekev)*&
                        erfc(tE)*exp(tL)/E32/&
                        sqrt(abs(1 - lambdum(j)))
            enddo
            dum = dedum/ekev
          call simp(nE1,dum,fedum,fldum(j))
          enddo
          call simp(numL1,dldum,fldum,Cn)
! Cn for isotropic slow down 
           do i=1,nE1
               endum(i)=Ea + (i-1)*dedum
               E32=(endum(i)/ekev)**1.5+Ec32
               tE=(endum(i)-E0)/Ed
               fedum(i)=2.0*PI_D*sqrt(endum(i)/ekev)*&
                        erfc(tE)/E32
            enddo
            dum = dedum/ekev
          call simp(nE1,dum,fedum,Cn)

         deallocate(endum)
         deallocate(fedum)
         deallocate(fldum)
         deallocate(lambdum)
         print *, 'Cn=',Cn
!
           open(8,file='result.out',status='unknown')
           vh=9.79*1e5/sqrt(prot)*sqrt(ekev*1e3)*sqrt(2.0)   ! EP velocity
                                                             ! with energy ekev
                                                             ! cm/sec
           omeg0 = 9.58e3*zprt/prot*bkg*1e3   ! rad/sec
           v_A0 = 2.18e+11/sqrt(prot0)/sqrt(ne0)*bkg*1e3 ! Alfven velocity
                                                        ! at axis cm/sec
           vA_bar = v_A0/vh ! normalized alfven velocity
!
           rhoh = vh/omeg0/rmaj/e
           ! normalized to minor radius
           rhoh_bar = 1.02e2*sqrt(prot)/zprt*sqrt(ekev*1e3)/(bkg*1e3)/rmaj
           !gyroradius normalized to major radius

!
           write(*,*) 'e=',real(e),'omega0=',omeg0, 'v_A0=',v_A0,&
                      'rhoh_bar:',real(rhoh_bar),'rhoh=',rhoh,&
                       'vA_bar:', vA_bar,'bkg=',bkg



        end subroutine profiledat

        subroutine profcoefdat
         integer nx0,k0,i !
         parameter(nx0=91)
         parameter(k0=5)
         common/knot_prof/knotx(nx0+k0)
         common/coeftfun/coefte(nx0),coefti(nx0)
         common/coefdenfun/coefden(nx0)
         common/coefbetafun/coefbeta(nx0)
         common/coefq/coefqfun(nx0)
         common/coefnbfun/coefnb(nx0)
         common/profiles_grids/nx,ksp
         integer :: nx,ksp
         real(r8) :: knotx,coefte,coefti
         real(r8) :: coefqfun,coefden,coefbeta,coefnb
!
         read(530,*)
         read(530,*)
         read(530,*) ksp,nx
         do i=1,nx+ksp
            read(530,*) knotx(i)
         enddo
!
!
         read(540,*)
         do i=1,nx
            read(540,*) coefte(i),coefti(i),coefnb(i)
         enddo
!
         read(550,*)
         do i=1,nx
            read(550,*) coefden(i),coefbeta(i),coefqfun(i)
         enddo
!
      end subroutine profcoefdat

      subroutine tensordat   ! coefficients for R and Z Fourier components
         integer nx0,nf0,nc0,k0,i,j !
         parameter(nx0=91)
         parameter(nf0=20)
         parameter(nc0=8)
         parameter(k0=5)
         common/knot_sp/knot(nx0+k0)
         common/jacob/coefc(nx0)
         common/rfdat/coefaf(nf0,nx0),coefaf0(nx0)
         common/zfdat/coefbf(nf0,nx0)
         common/tensor_grids/nx,ksp,nc,nf
         integer :: nx,ksp,nc,nf
         real(r8) :: knot,coefc,coefaf,coefaf0,coefbf
!
         read(53,*)
         read(53,*)
         read(53,*) ksp,nx,nc,nf
         do i=1,nx+ksp
            read(53,*) knot(i)
         enddo
!
        read(56,*)
       do i=1,nx
          read(56,*) coefc(i)
       enddo
!

!
       read(57,*)
       do i=1,nx
          read(57,*) coefaf0(i)
       enddo
!
        read(58,*)
       do j=1,nf
          do i=1,nx
             read(58,*) coefaf(j,i),coefbf(j,i)
          enddo
        enddo
!
       end subroutine tensordat

end module shared_mod
