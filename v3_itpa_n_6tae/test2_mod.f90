          module test2_mod
           !use mpi
           contains

             subroutine main
             !use omp_lib
             use mpi
             use vars
             use vars1
             use global
             use vars_k
             use vars_e
             use inrtype
             use splines
             use shared_mod
             use orbit_eq_mod
             implicit none
             integer, parameter :: nstep2=51,nstep3=201, nseca=4, ngyro=4
             !integer, parameter :: nq = 3 ! toroidal mode number
             integer :: ierr, numprocs, proc_num, points_per_proc
             real(8) :: start_time, end_time  
             integer :: i,j,m,l,k,p,modstep,nt1,istep,jn,ndist,nn,iE
             integer :: il,il1,ip,n,nprt1,itrap,nx,iname
             integer :: ifem, mfem, local_ifow, local_ifow1
             integer :: istart, iend, array_size, iarray
             integer ::  nprt_proc
             integer :: im, im1,ma, mb, idm
             integer :: st2,et2,count_rate
             integer(8) :: ik_proc,ik
             !common/resistive/etai,omg
             !common/nogr/nog
             !common/aspect/e!aspect rate:e=a/R_0
             !common/eigen_value/lam
             !real(8) :: etai,e,lam,omg
             !integer :: nog,nogp
             integer :: nogp
             real(8) :: x,h,x1,x2,wtheta,wphi,ttdum,dtime,dtdum0,dtdum1,s
             real(8) :: dum,dum0,dum1,dum2,dum3,dum4,dum5,&
                        dum6,dum7,dum8,dum9,dum10,dum11,xrdum
             real(8) :: dume,ttdum1,wtheta1,wphi1,wtheta2
             real(8) :: dela,a0,a1,b00
             real(8) :: pol0,pol,tf,pol00
             real(8) :: b,c,dl,ddx,x0,x00,v00
             real(8) :: a,w,vv,temp,sgn,temp1,temp2
             real(8) :: vstart(nvar), v(nmax), dv(nmax),vout(nmax)
             real(8) :: v1(nmax),dv1(nmax)
             real(8) :: tt(nvar),pp(nvar),thth(nvar)
             real(8) :: xx(nseca),yy(nseca)
             real(8) :: ttt,xxx,zzz,ppp,aaa,vvv,hhh
             real(8) :: emin, emax, phimin, phimax
             real(8), allocatable :: en0(:)  
             complex(8) :: zsquar, aint0
             real(8) :: fi(nstep2),fi2(nstep2),fi1(nstep3)
             real(8) :: fi3(nstep2),fi4(nstep2),xxi(nstep2)
             complex(8) :: gi(nstep2),gi1(nstep3)
             real(8) :: sdum,sdum1,sdum2
             complex(8) :: sdum3,sdum4,sdum5,fdum
             real(8) :: runt2,fac
             real(8) :: xa,xb,xw
             real(8) :: ra,rb
             real(8) :: eps
             complex(8) :: omg1
             real(8) :: Cn_proc,Cn1
             real(8) :: xla,xlb,xia,xib
             real(8) :: d1,d2,delta_r,delta_theta
             real(8) :: elapsed_time
             !real(8) :: delta_theta_j(ngyro),delta_r_j(ngyro)
             integer :: file_write,fow_info
             integer :: kct_info,orb_info,gyro_info
             complex(8) :: dum1_cplx,dum2_cplx
             complex(8) :: dum3_cplx,dum4_cplx
             complex(8) :: sdum_cplx


             call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
             call mpi_comm_rank(MPI_COMM_WORLD, proc_num, ierr)

           

            ! Ask the user for the number of points
           if (proc_num == 0) then
               print *, "Using ",numprocs," processors"
               call system_clock(st2,count_rate)
               call CPU_TIME(ctime0)

           endif

           ! 在子程序开始时插入 MPI_Barrier
            call mpi_barrier(MPI_COMM_WORLD, ierr)

             ! nprt : number of particles for computing
             ! x : time variable, h : time step
             ! x1:  starting 'time', x2: finial 'time'
             !  
             
             ! set computing parameters, particle number,wave parameters
             ! from fort.19
             ! nprt : same to mesh_xal_arc.dat if ndist = 2
           !if (proc_num == 0) then
!            equilibrium for device configuration
             call check_eq ! using gfile, bkg, rmaj input
!            parameters for device and particles
             call profiledat
             call setgrid
!            flux surface recording
             !call check3
             call check_flux1
             !call check_flux2
          ! endif


             read(1900,*)
             read(1900,*) nprt, a0, a1, x1, x2
             close(1900)
!
             allocate(R0(nprt),Z0(nprt),phi0(nprt),&
                      XR0(nprt),name0(nprt),dE0(nprt))
             allocate(xx0(nprt),theta0(nprt))
             allocate(en(nprt),ptch(nprt),ptcha(nprt),&
                      mube(nprt),sgnb(nprt),dxa(nprt))
             allocate(pz0(nprt),dLa(nprt),dpz0(nprt))
             allocate(w_theta(nprt),w_theta1(nprt),&
                      w_phi1(nprt),tau_theta1(nprt),&
                      psir_bar(nprt))
             w_theta=0.0D0
             w_theta1=0.0D0
             w_phi1=0.0D0
             tau_theta1=0.0D0
             psir_bar=1.0D-6
             ! note turn on file_write, record fort.6*** with freq info
             ! fort.7*** with orbit integral YY*, fort.8*** with dist
             ! info
             ! To record traj data, turn off kct_info and turn on
             ! orb_info, run code with single core due to traj data
             ! recording
             file_write = 0 ! 1 : write data 0:  not write 
             kct_info = 1 ! 1: obtain kct_4D term 0: not obtain
             orb_info = 0 ! 1: record traj1,traj2 0: not record
             fow_info = 1 ! 1: induce fow 0: not induce
             gyro_info = 1 ! 1: gyro-average 0: not gyro-average
             omg1 = omg*vA_bar ! omega in orbit code unit
             omg1 = omg1 + cmplx(0.0,etai)
             !omg1 = cmplx(real(omg),etai)
!            mode numbers
             nomode = 5
             nopmode = 28
             nogp = nog + 1
             if (allocated(mode)) then
                deallocate(mode)
             endif
            if (allocated(pmode)) then
               deallocate(pmode)
            endif
             allocate(mode(nomode),pmode(nopmode))
             if (allocated(aint_KCT4D)) then
                deallocate(aint_KCT4D)
             endif
             allocate(aint_KCT4D(0:nog,0:nog,nomode,nomode))
             allocate(temp_aint_KCT4D(0:nog,0:nog,nomode,nomode))
             if (allocated(aint_KCT4Dp)) then
                deallocate(aint_KCT4Dp)
             endif
             allocate(aint_KCT4Dp(0:nog,0:nog,nomode,nomode))
             !allocate(temp_aint_KCT4Dp(0:nog,0:nog,nomode,nomode))
             !allocate(ylmp(nopmode,nprt,nomode,0:nog))
             !ylmp=(0.0D0,0.0D0)
             aint_KCT4D=(0.0D0,0.0D0)
             temp_aint_KCT4D=(0.0D0,0.0D0)
             aint_KCT4Dp=(0.0D0,0.0D0)
             !temp_aint_KCT4Dp=(0.0D0,0.0D0)

!            harmonics setting
             mode(1) = 8
             mode(2) = 9
             mode(3) = 10
             mode(4) = 11
             mode(5) = 12
            ! mode(6) = 9
            !mode(7) = 10
             !mode(8) = 8
             !mode(9) = 9
             !mode(10) = 10

!            poloidal mode index shift value
             midx = mode(1) - 1

             pmode(1) = 13
             pmode(2) = 12
             pmode(3) = 11
             pmode(4) = 10
             pmode(5) = 9
             pmode(6) = 8
             pmode(7) = 7
             pmode(8) = 6
             pmode(9) = 5
             pmode(10) = 4
             pmode(11) = 3
             pmode(12) = 2
             pmode(13) = 1
             pmode(14) = 0
             pmode(15) = -1
             pmode(16) = -2
             pmode(17) = -3
             pmode(18) = -4
             pmode(19) = -5
             pmode(20) = -6
             pmode(21) = -7
             pmode(22) = -8
             pmode(23) = -9
             pmode(24) = -10
             pmode(25) = -11
             pmode(26) = -12
             pmode(27) = -13
             pmode(28) = -14

             print *, 'delta=',delta,'dtheta=',dtheta
             print *, 'pmode=', pmode
!            reset mode used
             nomode = 5
             nopmode = 28
!            h = (x2-x1)/nstep
             h = 1.0 ! dt: time step with unit of gyro-frequency
             nstep=int((x2-x1)/h) ! number of time step
             modstep = 10  ! step to record 
!
              
             nt1= int(nstep/modstep)
             allocate(tdum(nt1),Rdum(nt1),Zdum(nt1),psidum(nt1),&
                     phidum(nt1),thetadum(nt1),&
                     psirdum(nt1),vparaldum(nt1))
!
             allocate(f1(nprt,na),f2(nprt,na),f3(nprt,na),&
                      f4(nprt,na),f5(nprt,na),f6(nprt,na))

!
!            initial conditions of particles with fort.3 
!            R0: position in R with cm unit
!            Z0: position in Z with cm unit
!            phi0: toroidal position with rad unit
!            en: energy with kev unit
!            ndist: 1 for test case, 2 for mesh grids case

!            input initial data of particles
             ndist = 2
             if (ndist==1) then
                read(3,*)
                do j = 1,nprt
                   read(3,*) xx0(j), theta0(j), phi0(j), en(j), ptch(j),name0(j)
!                   normalise
                   en(j) = en(j)/ekev
                enddo
             endif
            close(3)
!
           

             if (ndist == 2) then
                open(89,file="mesh_xal_arc.dat",status='unknown')
                read(89,*)
                do j = 1,nprt
                   read(89,*) xx0(j), theta0(j), phi0(j),&
                              XR0(j), ptcha(j), mube(j), pz0(j),&
                              en(j),dLa(j),dpz0(j),&
                              dE0(j),name0(j),sgnb(j)

                   R0(j) = 1.0 + e*XR0(j) ! R = 1.0 + a/R0*X
                   Z0(j) = 0.0
                   ptch(j) = sin(ptcha(j)) ! vpara/v
                   en(j) = en(j)/ekev
                enddo
             end if
!
            close(89)
             Ec = Ec/ekev
             Ed = Ed/ekev
             E0 = E0/ekev
             
             ! reset phi    
             phimin = 0.0
             phimax = TWOPI_D
             call linspace(phimin,phimax,nprt,phi0)

     
             do j=1,na
                do i=1,nprt
!               co: -1 of sgn  counter: +1 of sgn for cylindrical 
!               co: direction of current
                    sgn =  1.0
!                 normalise
                  temp = b_nor(theta0(i),xx0(i))
                  en(i) = rhoh_bar**2*en(i)
                  if (ndist == 1) then
                      mube(i) = (1.0 - ptch(i)**2)/temp ! pitch: mu*B0/E
                  endif
                  f1(i,j) = xx0(i)
                  f2(i,j) = theta0(i)
                  f3(i,j) = phi0(i)
                  f4(i,j) = sgn*sqrt(2.0*en(i))*ptch(i) ! v_||
                  f5(i,j) = en(i) 
                  f6(i,j) = mube(i)
                enddo
             enddo

             dela = (a1-a0)/na
             x = 0.0

             open(unit = 1,file = 'traj1.dat')
             write(1,70) 
  70        format('# trajectory, t/omega_c,X,Z,phi,E,Pz')
             open(unit = 2, file = 'traj2.dat')
             write(2,71)
  71        format('# trajectory, t/omega_c,pol,theta,r/a,q,vpll')
             open(unit = 3, file = 'lost.dat')
             write(3,72)
  72        format('# trajectory, X0,Z0,phi0,alpha0,E,Pz,kth')
  73        format('# splines trajectory, t/omega_c,X,Z,phi,theta,&
                   r/a,p,itrap')
            open(unit = 5,file = 'phasefreq.dat')
             write(5,74)
  74        format('# x, E, mube, tau/(R0/vh), omega_theta,&
                  omega_phi, kth')

          
            open(8,file='orb_mod_tim.out',status='unknown')
            write(8,*) 'rmaj=',rmaj
  75        format('# XR0, r/a=x, dx/dpsi, kth')
            !open(unit = 9,file = 'dxp.dat')
            !write(9,75)
  76        format('# dpsidx, I1, I2, dx, dE,  kth')
            open(unit = 10,file = 'I12.dat')
            write(10,76)
           

  77        format('# realyy*, imagyy*, im,im1,il,il1,p')
            open(unit = 12,file = 'Y6D.dat')
            write(12,77)

            open(unit = 13, file = 'Y6D_1.dat')
            write(13,77)
78          format('# nprt1 nlim/rows_of_Y6D_1.dat vA/vh')


441        format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
442        format(1F16.6,1x,1F12.6,1x,1F12.6,1x,1F12.6,1x,1F12.6)
443        format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                  1E16.9,1x,I8,I6)
444        format(1F16.6,1x,1E16.8,1x,1E16.8,1x,I8)
445        format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,I8)
446        format(1E16.9,1x,1E16.9,1x,I6)
447        format(1E16.9,1x,1E16.9,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4)
448        format(1E16.9,1x,1E16.9,1x,1E16.9,1x,1E16.9,1x,&
                  1E16.9,1x,1E16.9,1x,I6)



            ! Broadcast to all procs; everybody gets the value of nprt from
            ! proc 0
            call mpi_bcast(nprt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)


!           time advance
           l = 0   ! l: particle name0 choosen to run, l=0: all particles to run
           if (l==0) then
             do j=1,na
                a = a1-(j-1)*dela
                nprt1 = 0 ! num of confined ions 20230313 yulm 
                xw = 1.0
                ik = 0
                !pw = inrcsval(gridx2,cscoefxpsi,xw,0)     
                print *, 'pw=',pw

               ! Determine how many points to handle with each proc
                points_per_proc = (nprt + numprocs - 1)/numprocs
                if (proc_num == 0) then   ! Only one proc should print to avoid clutter
                   print *, "points_per_proc = ", points_per_proc
                end if


                nprt_proc = 0
                ik_proc= 0
                Cn_proc = 0.0
              ! Determine start and end index for this proc's points
               istart = proc_num * points_per_proc + 1
               iend = min((proc_num + 1)*points_per_proc, nprt)

              ! Diagnostic: tell the user which points will be handled by which proc
         print '("Process ",i6," will take i = ",i8," through i = ",i8)', &
                 proc_num, istart, iend
                 
!
               array_size = iend - istart + 1
               allocate(ylmp(nopmode,array_size,nomode,0:nog))
               allocate(fylmp(nopmode,array_size,nomode,0:nog))
               !allocate(fylmpp(nopmode,array_size,nomode,0:nog))
                ylmp=(0.0D0,0.0D0)
                fylmp=(0.0D0,0.0D0)
               !fylmpp=(0.0D0,0.0D0)
!
               iname0 = 405
               do l=istart, iend
                    vstart(1) = f1(l,j) ! x
                    vstart(2) = f2(l,j) ! theta
                    vstart(3) = f3(l,j) ! phi: toroidal angle
                    vstart(4) = f4(l,j) ! v_||: v parallel
                    vstart(5) = f5(l,j) ! E: energy
                    vstart(6) = f6(l,j) ! mube: Lambda
                    do i=1,nvar
                       v(i) = vstart(i)
                    enddo
!
!
                    nprt_proc = nprt_proc + 1 
!
                   iname = name0(l) 
                   !print* , 'iname1=', iname, 'pro:', proc_num
                   iarray = iname - istart + 1
                   !print* , 'iname1=', iname, 'pro:', proc_num,'iarray:',iarray
                    ! initial positions record
                    temp = b_nor(v(2),v(1))
                    tmpfpol = inrcsval(gridx2,cscoefxpsi,v(1),0)
!                   particle energy with kev unit
                    temp1 = (v(4)**2/2.0 + v(5)*v(6)*temp)/&
                              rhoh_bar/rhoh_bar*ekev
!                   particle toroidal canonical momentum  
                    temp2 = tmpfpol*e**2 - v(4)/temp   ! P_phi 


                    x=x1  ! reset initial time for each ion 
                    m=1  ! initial index for sp interp
                    i=0
                    itrap = 0 ! passing itrap = 0, trapped itrap = 1 
                    zzz=0.0
                    p=0
                    xw = 1.0
                    !pw = inrcsval(gridx2,cscoefxpsi,xw,0)
                    ! initial positions record for sp interp
                    tdum(m)=x
                    Rdum(m)=rf(v(2),v(1))
                    Zdum(m)=zf(v(2),v(1))
                    phidum(m)=v(3)
                    thetadum(m)=v(2)
                    psirdum(m)=v(1)
                    vparaldum(m)=v(4)
!                advance of particle with RK4 scheme 
                 do k = 1,nstep
                       call derivs1(x,v,dv)
                       call rk4(v,dv,nvar,x,h,v)
!                   lost condition
                    if (v(1)>xw) then
                        nprt_proc = nprt_proc - 1 ! neglect loss ions 20230313
                        !if (file_write ==1) then
                        if (orb_info ==1) then 
                           write(9000+proc_num,448)  XR0(l), Z0(l), &
                                     phi0(l), ptcha(l),&
                                     en(l), pz0(l),name0(l)
                         endif
                        exit  ! particle loss to last flux surface
                    endif

!                   singular/crossing axis point
                    if (v(1)/xw<=0.002) then
                        nprt_proc =nprt_proc -1 
                        exit  ! particle cross axis
                    endif


!                   particle data to record
                     if (mod(k,modstep)==0) then
                        m=m+1

                        temp = b_nor(v(2),v(1))
                        tmpfpol = inrcsval(gridx2,cscoefxpsi,v(1),0)
                        !tmpfpol = (10.0*log(v(1)**2 +21.0/11.0))/11.0-&
                        !       (10.0*log(21.0/11.0))/11.0

!                       particle energy with kev unit
                        temp1 = (v(4)**2/2.0 + v(5)*v(6)*temp)/&
                                rhoh_bar/rhoh_bar*ekev  
!                       particle toroidal canonical momentum normalized
!                       to R0**2*omega_c
                         temp2 = tmpfpol*e*e - &
                                v(4)/b_nor(v(2),v(1))   ! P_phi 
                       
                        if (orb_info ==1) then
                          if (iname == iname0) then
                          !write(1,441) x, rf(v(2),v(1)), zf(v(2),v(1)), &
                          !            v(3), temp1, temp2 ! turn off if nprt particles run suggested
                          !write(2,441) x, tmpfpol, v(2), v(1), qfun(v(1)), v(4)
                            write(6000+proc_num,448) x, rf(v(2),v(1)),zf(v(2),v(1)), &
                                          v(3), temp1, temp2, iname ! turn off if nprt particles run suggested
                            write(7000+proc_num,448) x, tmpfpol, v(2),v(1), qfun(v(1)),&
                                           v(4),iname
                          endif
                        endif
                        if ( sign(1.0,v(4)*vvv)<0.0) then
                            itrap = 1
                        endif
!
!                      poloidal frequency and toroidal 
!                      frequency finding as recording
!                      duration of crossing mid plane 
                       tdum(m)=x  ! t
                       Rdum(m)=rf(v(2),v(1))
                       Zdum(m)=zf(v(2),v(1))  
                       phidum(m)=v(3) ! phi
                       thetadum(m)=v(2)
                       psirdum(m)=v(1) ! x
                       vparaldum(m)=v(4) ! v||

!                      transit time finding
!                      one transit data to record
                       if (m>1) then
                              if (zzz>0 .and. Zdum(m)<0) then   ! cross the mid plane by change sign of z
                                 i=i+1 ! number for circulating poloidal
                                       ! cross section
                            ! interporate time linearly  for z=0
                                 tt(i) = ttt*(-Zdum(m))/(zzz-Zdum(m)) + &
                                         x*(-zzz)/(Zdum(m)-zzz)
                            ! interporate phi linearly for z=0
                                 pp(i) = ppp*(-Zdum(m))/(zzz-Zdum(m)) + &
                                         v(3)*(-zzz)/(Zdum(m)-zzz)
                            ! interporate theta linearly for z=0
                                 thth(i) = hhh*(-Zdum(m))/(zzz-Zdum(m)) + &
                                         v(2)*(-zzz)/(Zdum(m)-zzz)
                            
                              end if
!
                              !if (v(4)*vvv<0.0.or.(v(6)>=1.0)) then ! condition of trapped for v|| reverse
                              !if (v(4)*vvv<0.0) then
                              !   itrap = 1
                              !endif

!
                       end if
!

!
!                            data of last time to keep                      
                             ttt = x
                             xxx = rf(v(2),v(1))
                             zzz = zf(v(2),v(1))
                             hhh = v(2)
                             ppp = v(3)
                             aaa = v(1)
                             vvv = v(4)

                     end if

!
                              
                     if (x+h==x) &
                          pause "stepsize vot signficant, in rkdumb"
                          x = x + h ! time/step advance
! 
                      
!                    period and frequencies finding
                     if (i==3) then
                           ttdum = abs(tt(1)-tt(2))
                           if (itrap.eq.1) then
                              wtheta = 2.0*pi/ttdum ! trapped ions
                           else
                              wtheta = 2.0*pi/ttdum*sign(1.0,ptcha(iname))! passing ions
                           endif
                           !wtheta2 = (thth(2)-thth(1))/ttdum
                           ! this part repeated
                           wtheta=wtheta/rhoh_bar/sqrt(2.0) ! normalized to vh/R0
                           wphi =( (pp(2)-pp(1))/(ttdum))!for v_||<0 
                           wphi = wphi/rhoh_bar/sqrt(2.0) ! normalized to vh/R0
                           ttdum = ttdum * rhoh_bar*sqrt(2.0) ! normalized to R0/vh
!                          initial position for w/o FOW
                           tmpfpol = inrcsval(gridx2,cscoefxpsi,v(1),0)
                           dume = temp1/ekev
                           pol00 = tmpfpol/pw

!                           phase positions  and frequencies
!                           recover data units for spline-interp
                            ttdum = abs(tt(1)-tt(2))
                           !if (itrap.eq.1) then
                           !   wtheta = 2.0*pi/ttdum ! trapped ions
                              !wtheta = 2.0*pi/ttdum*sign(1.0,ptcha(iname)) 
                           !else
                           !   wtheta = 2.0*pi/ttdum*sign(1.0,ptcha(iname)) ! passing ions
                              !wtheta = 2.0*pi/ttdum
                           !endif
                            wtheta = 2.0*pi/ttdum ! shen wei fixed 250321
                            wtheta2 = (thth(2)-thth(1))/ttdum
                            wphi =((pp(2)-pp(1))/(ttdum))
                            dum8 = tt(1)
                            dum9 = tt(2)
                            dtime = ttdum/real(nstep2-1,sp)
                            ! check bounce time
                            !write(1003,*) 'tt1:',tt(1),'tt2:',tt(2),&
                            !              'ptcha:',ptcha(iname),'idx:',iname
                            ! write(1003,*) 'itrap:',itrap,'wtheta:',wtheta,'iname:',iname
                      !write(8000+proc_num,*) 'wtheta:',wtheta,'wtheta2:',wtheta2,&
                      !              'iname:',iname,'itrap:',itrap
!                         coefs of spline-interp for one transit 
                          allocate(cscoefR(4,m-2),cscoefZ(4,m-2),cscoefphi(4,m-2),&
                                  cscoeftheta(4,m-2),cscoefpsir(4,m-2),&
                                  cscoeftime(4,m-2),cscoefvparal(4,m-2))
                          call inrcsnak(tdum(1:m-1),Rdum(1:m-1),cscoefR) 
                          call inrcsnak(tdum(1:m-1),Zdum(1:m-1),cscoefZ)
                          call inrcsnak(tdum(1:m-1),phidum(1:m-1),cscoefphi)
                          call inrcsnak(tdum(1:m-1),thetadum(1:m-1),cscoeftheta)
                          call inrcsnak(thetadum(1:m-1),tdum(1:m-1),cscoeftime)
                          call inrcsnak(tdum(1:m-1),psirdum(1:m-1),cscoefpsir)
                          call inrcsnak(tdum(1:m-1),vparaldum(1:m-1),cscoefvparal)
!
                          allocate(tdum1(nstep2),Rdum1(nstep2),Zdum1(nstep2),&
                                   phidum1(nstep2),phi1dum1(nstep2),&
                                   thetadum1(nstep2),psirdum1(nstep2))
                          allocate(psirdum2(ngyro,nstep2),thetadum2(ngyro,nstep2))
                          
!
                          do p=1,nstep2
                             dum0=(p-1)*dtime+dum8
                             dum5=inrcsval(tdum(1:m-1),cscoefpsir,dum0,0)!psir:r/a
                             fi(p) = dum5
                          enddo
!
!                          r averaging
                           call simp(nstep2,dtime,fi,x00)
                           x00 = x00/ttdum ! averaged r
                           psir_bar(iname)=x00



                
!
!                         poloidal and toroidal freq 
!                         this part repeated for check
                        !if ((itrap==1)) then
                        !    ttdum = abs(tt(1)-tt(2))
                        !    wtheta = 2.0*pi/ttdum
                        !   ! wphi = abs( (pp(2)-pp(1))/(ttdum)) ! trapped
                        !    wphi = ( (pp(2)-pp(1))/(ttdum)) ! trapped
                                                            ! should be positive ?
                        !    dum8 = tt(1)
                        !    dum9 = tt(2)
                        !    dtime = ttdum/real(nstep2-1,sp)
                        !    ! check bounce time 
                        !    !write(1003,*) 'tini:',dum8, 'tfin:',dum9
                        !    !write(1003,*) 'tau_b:',ttdum,'idx:',iname
                        !    !write(1003,*) 'wphi:',wphi,'itrap:',itrap
                        !    !write(1003,*) 'Rmin:',MinVal(Rdum(1:m-1)),&
                        !    !              'Rmax:',maxval(Rdum(1:m-1))
                        !  !if (itrap.eq.1) then
!                       !   bounce angle of trapped ions using secant method
                        !  !   x0 = maxval(thetadum(1:m-1))
                        !  !   dl = 1e-6 ! secant trunction error
                        !  !   b = 0.0
                        !  !   c = PI_D
                        !  !   ddx = (c-b)/100.0
                        !  !   call secant3(dl,x0,ddx,istep,x00,v(6))

                        !  !   dum10 = x0  ! theta_b
                        !  ! else
                        !  !   dum10 = maxval(thetadum(1:m-1)) ! theta_b
                        !  ! endif
                           
!
!
                        ! itrap = 0
                        ! sp-interp for transit time
                       ! else  
                       !   if (abs(vstart(2))==0.0) then
                       !       dum6 = 0.0 ! ini theta 
                       !       if (ptch(l) > 0) then
                       !          dum7 = TWOPI_D ! fin theta passing particles
                       !       else
                       !          dum7 = -TWOPI_D
                       !       endif
                       !       dum8 = inrcsval(thetadum(1:m-1),cscoeftime,dum6,0) 
                       !       dum9 = inrcsval(thetadum(1:m-1),cscoeftime,dum7,0)
                       !       ttdum = abs(dum9-dum8)
                       !       ! check transit time for X > 0
                       !       !write(1003,*) 'dum6:',dum6,'tini:',&
                       !       !               dum8, 'tfin:',dum9,&
                       !       !               'idx:',iname
                       !       !print *, 't0=',dum8,'t2pi=',dum9
                       !       !write(1003,*) 'tau_t_sp:',ttdum,'idx:',iname
                       !       dtime = ttdum/real(nstep2-1,sp)
                       !       wtheta=TWOPI_D/ttdum ! definition
                       !       dum4 = inrcsval(tdum(1:m-1),cscoefphi,dum6,0)
                       !       dum5 = inrcsval(tdum(1:m-1),cscoefphi,dum9,0)
                       !       !wphi = abs(dum5-dum4)/ttdum
                       !       wphi = (dum5-dum4)/ttdum ! take account co-/ctr- passing
                       !       !write(1003,*) 'wphi:',wphi

                       !   else
                       !       dum6 = theta0(l)
                       !       if (ptch(l) > 0) then
                       !          dum7 = TWOPI_D+theta0(l) ! passing particles
                       !       else
                       !          dum7 = -PI_D
                       !       endif
                            
                       !       dum8 = inrcsval(thetadum(1:m-1),cscoeftime,dum6,0)
                       !       dum9 = inrcsval(thetadum(1:m-1),cscoeftime,dum7,0)
                       !       ttdum = abs(dum9-dum8)
                       !       ! check transit time for X < 0
                       !       !write(1003,*) 'dum6:',dum6,'tini:',&
                       !       !               dum8, 'tfin:',dum9,&
                       !       !               'idx:',iname
                       !       !write(1003,*) 'tau_t_sp:',ttdum,'idx:',iname
                       !       dtime = ttdum/real(nstep2-1,sp)
                       !       wtheta=TWOPI_D/ttdum
                       !       dum4 = inrcsval(tdum(1:m-1),cscoefphi,dum6,0)
                       !       dum5 = inrcsval(tdum(1:m-1),cscoefphi,dum9,0)
                       !       !wphi = abs(dum5-dum4)/ttdum
                       !       wphi = (dum5-dum4)/ttdum ! take account co-/ctr- passing
                       !       !write(1003,*) 'wphi:',wphi
                       !   endif
!
                       ! endif


                           ! using sp-interp to obtain, units in  KAEC code
                           ! wtheta1=wtheta/rhoh_bar/sqrt(2.0) ! normalized to vh/R0
                           ! wphi1= wphi/rhoh_bar/sqrt(2.0) ! normalized to vh/R0
                           ! ttdum1 = ttdum * rhoh_bar*sqrt(2.0) ! normalized to R0/vh
                           ! note ttdum, wtheta, wphi in current code orbit_eq unit
                             ttdum1 = ttdum*rhoh_bar*sqrt(2.0) 
                             wtheta1 = wtheta/rhoh_bar/sqrt(2.0)
                             wphi1 = wphi/rhoh_bar/sqrt(2.0)
                             tau_theta1(iname) = ttdum1
                             w_theta1(iname) = wtheta1
                             w_phi1(iname) = wphi1 
                             if (file_write ==1) then
!                       '# x(r/a), E, mube, tau/(R0/vh), omega_theta,omega_phi, kth, itrap'
                                 write(6000+proc_num,443) x00,temp1,vstart(6),&
                                            ttdum1,wtheta1,wphi1,name0(l),itrap
                              endif

                            

!
                         ! sp-interp orb positions 
                         do p=1,nstep2
                            dum0=(p-1)*dtime+dum8
                            tdum1(p)=dum0
                            if (p==1) then
                            phidum1(p)=inrcsval(tdum(1:m-1),cscoefphi,dum0,0) !! phi
                            thetadum1(p)=inrcsval(tdum(1:m-1),cscoeftheta,dum0,0) !! theta
                            phi1dum1(p)=phidum1(p)-wphi*tdum1(p)!phi_tidle = phi- omega_phi*t
                            endif
                           ! reset starts to zeros
                            phidum1(p)=inrcsval(tdum(1:m-1),cscoefphi,dum0,0)-phidum1(1) !! phi
                            thetadum1(p)=inrcsval(tdum(1:m-1),cscoeftheta,dum0,0)-thetadum1(1) !! theta
                            psirdum1(p)=inrcsval(tdum(1:m-1),cscoefpsir,dum0,0) !! psir: r/a
                            phi1dum1(p)=phidum1(p)-wphi*tdum1(p)-phi1dum1(1)!phi_tidle= phi- omega_phi
                            xrdum = psirdum1(p)
                            
                            
                            if (gyro_info == 1) then
                               dume = en(iname)/rhoh_bar/rhoh_bar ! E/ekev
                               dum4 = thetadum1(p) ! theta
                               dum3 = sqrt(2.0 * vstart(6) * dume / b_nor(dum4, xrdum)) * rhoh_bar ! rho_c
                               !dum3 = 0.0
                               dum6 = sqrt(grr(dum4, xrdum))
                               fac = 1.0
                               do  idm = 1, ngyro
                                 temp = TWOPI_D/float(ngyro)*float(idm-1)
                                 !call gyrosub(xrdum, vstart(6), dume, dum4, temp,&
                                 delta_r = dum3/e*cos(temp)*dum6 
                                 dum1 = xrdum + delta_r
                                 ! boundary conditions
                                 if (dum1 > 1.0) then 
                                    dum1 = xrdum - delta_r
                                 endif
                                 if (dum1 < 0.0) then
                                    dum1 = xrdum - delta_r
                                  endif
                                 psirdum2(idm,p) = dum1 ! save positions
                                 delta_theta = dum3/e/dum6*&
                                     (cos(temp)*grt(dum4,xrdum)-&
                                     sin(temp)/rf(dum4,xrdum)/&
                                     jr2(xrdum))
                                 thetadum2(idm,p) = dum4 + delta_theta
      
                               enddo
                            end if
                            ! record sp traj1
                           if (orb_info ==1) then
                             if (iname==iname0) then
                               write(35,441) tdum1(p), rf(thetadum1(p),psirdum1(p)),&
                                         zf(thetadum1(p),psirdum1(p)), &
                                          phidum1(p), psirdum1(p), temp2 ! turn off if nprt particles run suggested
                              endif
                           endif

                         enddo
                        
                         if (gyro_info == 1) then
                            ra = minval(psirdum2)
                            rb = maxval(psirdum2)
                         else
                           ra = minval(psirdum1)
                           rb = maxval(psirdum1)
                         endif
                         !print *, 'ra=',ra,'rb=',rb,'iname:',iname,&
                         !          'delta_r_j:',delta_r_j(idm)
!
                     eps=1.0D-6
                     eps=xini+eps
                              
                     if (kct_info == 1) then
                          xrdum = x00 ! averaged r for w/o FOW
                       if (fow_info == 1) then ! turn on fow eff 20240210
                          do im = 1, nomode ! im-> m
                                do il = 0, nog ! il->l
                                !  range of finite element h_l
                                   if(il.eq.0)then ! fixed 1->0 20230304 yulm
                                       xla=eps
                                       xlb=xg(il+2)
                                   else if(il.eq.(nog-1))then
                                       xla=xg(il-2)
                                       xlb=xg(nog)
                                   else
                                       xla=xg(il-2)
                                       xlb=xg(il+2)
                                   endif
                                   
                                    if(ra <= xlb .and. rb >= xla) then ! nes FOW cond                                     
                                       xrdum = psir_bar(iname)
                                       dume = en(iname)/rhoh_bar/rhoh_bar ! E/ekev
                                       pol0 = inrcsval(gridx2,cscoefxpsi,xrdum,0)/pw !<psi>/pw
                                       dum2 = pfpe(dume,mube(iname),pol0)
                                       dum3 = pfpz1(dume,mube(iname),pol0)
                                       dum4 = dLa(iname)
                                        !  real or complex omg ?
                                        !fdum = (dume)**3*dE0(iname)/ekev*dpz0(iname)/&
                                        !            e/e*dum4*tau_theta1(iname)*&
                                        !            (omg*dum2-float(ntor)*dum3)
                                        fdum = (dume)**3*dE0(iname)/ekev*dpz0(iname)/&
                                               e/e*dum4*tau_theta1(iname)*&
                                               (real(omg1)*dum2-float(ntor)*dum3)
                                        !fdum = (dume)**3*dE0(iname)/ekev*dpz0(iname)/&
                                        !       e/e*dum4*tau_theta1(iname)*&
                                        !       (real(omg1)*dum2) ! test

                                    endif

                                      do ip = 1, nopmode ! ip-> p for matrix element M'_{i,l}^{k,m}

                                        if(ra <= xlb .and. rb >= xla) then ! nes FOW cond

                                         fi=0.0
                                         fi2=0.0
                                         gi=(0.0,0.0)
                                         !$OMP PARALLEL DO PRIVATE(dum0,dum1,dum2,dum3,dum4,dum5,&
                                         !$OMP& dum8,dum11,xrdum) &
                                         !$OMP& REDUCTION(+:local_ifow,sdum,sdum1,&
                                         !$OMP& sdum2,sdum3)
                                         sdum_cplx = (0.0,0.0)
                                         fac = 1.0
                                         do p=1,nstep2
                                             dum0=tdum1(p)
                                             dum4=thetadum1(p) ! theta
                                             dum5=psirdum1(p) ! psir: r/a
                                             ! phi_tidle = phi -  omega_phi*t
                                             dum11 = phi1dum1(p)
                                             xrdum = dum5
                                             !if ((xrdum <= xlb .and.xrdum >=xla)) then ! eno FOW condition 
                                             !!!!! gyro average
                                                if (gyro_info == 1) then
                                                  !call cpu_time(start_time)  ! 开始计时
                                                    dum1_cplx = (0.0,0.0)
                                                    dum2_cplx = (0.0,0.0)
                                                    do  idm = 1, ngyro
                                                       temp = TWOPI_D/float(ngyro)*float(idm-1)
                                                       !call gyrosub(xrdum, vstart(6), dume, dum4, temp,&
                                                       !     delta_r_j(idm), delta_theta_j(idm))
                                                       dum1 = psirdum2(idm,p)

                                                       if ((dum1 <= xlb .and.dum1 >=xla)) then ! eno FOW condition
                                                          dum2 = thetadum2(idm,p)
                                                          dum1_cplx = dum1_cplx + fem(il,mode(im),dum1)*&
                                                               exp(cmplx(0.0,-dum2*float(mode(im))))
                                                          dum2_cplx = dum2_cplx + femp(il,mode(im),dum1)*&
                                                               exp(cmplx(0.0,-dum2*float(mode(im))))  
                                                       endif
                                                     enddo
                                                   !  call cpu_time(end_time)  ! 结束计时
                                                   !  elapsed_time = end_time - start_time
                                                   !  print *, 'Loop & execution&
                                                   !         time: ', elapsed_time, ' seconds'
                                                    dum1_cplx = dum1_cplx/float(ngyro)
                                                    dum2_cplx = dum2_cplx/float(ngyro)
                                                    !if ((dum1_cplx /= cmplx(0.0,0.0)) .or. &
                                                    !   (dum2_cplx /= cmplx(0.0,0.0))) then
                                                       gi(p) = (2.0/b_nor(dum4,xrdum)-vstart(6))*&
                                                            1.0/b_nor(dum4,xrdum)/jf(dum4,xrdum)*&
                                                            (-kappa_t(dum4,xrdum)*dum2_cplx + &
                                                            cmplx(0.0,-float(mode(im)))*&
                                                            kappa_r(dum4,xrdum)*dum1_cplx)*&
                                                            exp(cmplx(0.0,float(ntor)*dum11))*&
                                                            exp(cmplx(0.0,-float(pmode(ip))*wtheta*dum0)) 
                                                       ! another way to integral with summation
                                                       ! sdum_cplx = sdum_cplx + gi(p)
                                                    !endif
                                                else
                                                    if ((xrdum <= xlb .and.xrdum >=xla)) then ! eno FOW condition 
                                                      fi(p)=realfuny(vstart(6),xrdum,dum0,&
                                                          dum4,dum11,wtheta,ntor,il,&
                                                          mode(im),pmode(ip))
                                                      fi2(p)=imagfuny(vstart(6),xrdum,dum0,&
                                                          dum4,dum11,wtheta,ntor,&
                                                          il,mode(im),pmode(ip))
                                                      ! another way to integral with summation
                                                      !sdum_cplx = sdum_cplx + cmplx(fi(p), fi2(p))
                                                    endif
                              
                                                endif
                                             !endif
                                 
                                         enddo
                                         !$OMP END PARALLEL DO 
                                             ik_proc = ik_proc + 1
                                            if (gyro_info == 1) then
                                                 call simp_complex(nstep2, dtime, gi, sdum_cplx)
                                                 ylmp(ip, iarray, im, il) = sdum_cplx / ttdum
                                                 !another way to integral with summation
                                                 !ylmp(ip,iarray,im,il) = sdum_cplx / nstep2
                                                 
                                            else
                                                 call simp(nstep2, dtime, fi, sdum)
                                                 call simp(nstep2, dtime, fi2, sdum1)
                                                 ylmp(ip, iarray, im, il) = cmplx(sdum, sdum1) / ttdum
                                                 !another way to integral with summation
                                                 !ylmp(ip,iarray,im,il) = sdum_cplx / nstep2
                                            endif
                                              fylmp(ip,iarray,im,il) = ylmp(ip,iarray,im,il)/&
                                                         (-omg1+float(pmode(ip))*w_theta1(iname)+&
                                                         float(ntor)*w_phi1(iname))*fdum
                                             ! fylmpp(ip,iarray,im,il)=ylmp(ip,iarray,im,il)/&
                                             ! (-omg1+float(pmode(ip))*w_theta1(iname)+&
                                             !  float(ntor)*w_phi1(iname))**2*fdum/2.0/omg                            
                                            if (file_write ==1) then
                                             ! '# Re(YY*), Im(YY*), m, k, l, i, p '
                                             !write(7000+proc_num,447) sdum,sdum1,mode(im),mode(im),&
                                             !          il,il,pmode(ip)
                                             endif
                                         
                                        endif
                                      enddo

                                enddo
                          enddo
                       else ! turn off fow eff 20240210 add
                          do im = 1, nomode ! im-> m
                             !do im1 = 1, nomode ! im1->k
                                do il = 0, nog ! il->l
                                   !do il1 = 0, nog  ! il1->i

                                      call xrange(il,il1,xa,xb,nx)

                                      if (xrdum < xb .and. xrdum > xa)then !  w/o FOW
                                          xrdum = psir_bar(iname)
                                          dume = en(iname)/rhoh_bar/rhoh_bar ! E/ekev
                                          pol0 = inrcsval(gridx2,cscoefxpsi,xrdum,0)/pw !<psi>/pw
                                          dum2 = pfpe(dume,mube(iname),pol0)
                                          dum3 = pfpz1(dume,mube(iname),pol0)
                                          dum4= dLa(iname)
                                          fdum = (dume)**3*dE0(iname)/ekev*dpz0(iname)/&
                                                  e/e*dum4*tau_theta1(iname)*&
                                                  (omg1*dum2-float(ntor)*dum3)
                                      endif

                                      do ip = 1, nopmode ! ip-> p for matrix element M'_{i,l}^{k,m}

                                        if (xrdum < xb .and. xrdum > xa)then !  w/o FOW

                                         ik_proc = ik_proc + 1
                                         fi=0.0
                                         fi2=0.0
                                         fi3=0.0
                                         fi4=0.0

                                         !$OMP PARALLEL DO PRIVATE(dum0,dum1,dum2,dum3,dum4,dum5,&
                                         !$OMP& dum8,dum11,xrdum) &
                                         !$OMP& REDUCTION(+:local_ifow,sdum,sdum1,&
                                         !$OMP& sdum2,sdum3)
                                          do p=1,nstep2
                                             dum0=tdum1(p)
                                             dum4=thetadum1(p) ! theta
                                             dum5=psirdum1(p) ! psir:r/a
                                             ! phi_tidle = phi - omega_phi*t
                                             !dum11 = phi1dum1(p)
                                             dum11 = phidum1(p)
                                             !xrdum = dum5
                                             fi(p)=realfuny(vstart(6),xrdum,dum0,&
                                                          dum4,dum11,wtheta,ntor,il,&
                                                          mode(im),pmode(ip))
                                             fi2(p)=imagfuny(vstart(6),xrdum,dum0,&
                                                          dum4,dum11,wtheta,ntor,&
                                                          il,mode(im),pmode(ip))
                                             !fi3(p)=realfuny(vstart(6),xrdum,dum0,&! conjecture part
                                             !             dum4,dum11,wtheta,ntor,il1,&
                                             !             mode(im1),pmode(ip))/ttdum
                                             !fi4(p)=-imagfuny(vstart(6),xrdum,dum0,&
                                             !             dum4,dum11,wtheta,ntor,&
                                             !             il1,mode(im1),pmode(ip))/ttdum
                                          enddo
                                         !$OMP END PARALLEL DO
                                          ik_proc = ik_proc + 1
                                          call simp(nstep2,dtime,fi,sdum)
                                          call simp(nstep2,dtime,fi2,sdum1)
                                          !call simp(nstep2,dtime,fi3,sdum2)
                                          !call simp(nstep2,dtime,fi4,sdum3)
                                          sdum = sdum / ttdum
                                          sdum1 = sdum1 / ttdum
                                          ylmp(ip,iarray,im,il)=cmplx(sdum,sdum1)
                                          fylmp(ip,iarray,im,il)=ylmp(ip,iarray,im,il)/&
                                           (-omg1+float(pmode(ip))*w_theta1(iname)+&
                                            float(ntor)*w_phi1(iname))*fdum
                                           !fylmpp(ip,iarray,im,il)=ylmp(ip,iarray,im,il)/&
                                           !(-omg1+float(pmode(ip))*w_theta1(iname)+&
                                           ! float(ntor)*w_phi1(iname))**2*fdum/2.0/omg     
                                          !temp_aint_KCT4D(il,il1,im,im1)=temp_aint_KCT4D(il,il1,im,im1)+&
                                          !           1.0/2.0*cmplx(sdum4,sdum5)/(-omg1+float(pmode(ip))*wtheta1+&
                                          !           float(ntor)*wphi1)*fdum

                                          if (file_write ==1) then
                                            ! '# Re(YY*), Im(YY*), m, k, l, i, p '
                                             !write(7000+proc_num,447)sdum,sdum1,mode(im),mode(im),&
                                             !              il,il,pmode(ip)
                                          endif



                                        else
                                          sdum4 = 0.0
                                          sdum5 = 0.0
                                          temp_aint_KCT4D(il,il1,im,im1)=temp_aint_KCT4D(il,il1,im,im1)+0.0

                                        endif


                                      enddo
                                   !enddo
                                enddo
                             !enddo
                          enddo


                       endif
                     endif
!
                          deallocate(cscoefR)
                          deallocate(cscoefZ)
                          deallocate(cscoefphi)
                          deallocate(cscoeftheta)
                          deallocate(cscoefpsir)
                          deallocate(cscoeftime)
                          deallocate(cscoefvparal)
                          deallocate(tdum1)
                          deallocate(Rdum1)
                          deallocate(Zdum1)
                          deallocate(phidum1)
                          deallocate(phi1dum1)
                          deallocate(thetadum1)
                          deallocate(psirdum1)
                          deallocate(psirdum2)
                          deallocate(thetadum2)

                        exit ! one transit to exit
                     end if
                 enddo
!             
                      
               
                tdum=0.0
                Rdum=0.0
                Zdum=0.0
                phidum=0.0
                thetadum=0.0
                psirdum=0.0
                vparaldum=0.0
!
               enddo
!
               deallocate(tdum)
               deallocate(Rdum)
               deallocate(Zdum)
               deallocate(phidum)
               deallocate(thetadum)
               deallocate(psidum)
               deallocate(psirdum)
               deallocate(vparaldum)
!             
           do l=istart, iend 
                iname = name0(l)
                !print *, 'iname2=',iname, 'pro:', proc_num
                iarray = iname - istart + 1
                !temp = b_nor(theta0(iname),xx0(iname))
                !if (mube(iname)*temp < 1.0) then
                !    Cn_proc = Cn_proc + &
                !              temp/sqrt(1.0-mube(iname)*temp)*&
                !              dLa(iname)
                ! endif
                if (kct_info == 1) then  !
                    if (w_theta1(iname) /= 0.0D0) then ! consider
                                                       !confined particles
                       !xrdum = psir_bar(iname)
                       !dume = en(iname)/rhoh_bar/rhoh_bar ! E/ekev
                       !pol0 = inrcsval(gridx2,cscoefxpsi,xrdum,0)/pw !<psi>/pw
                       !dum2 = pfpe(dume,mube(iname),pol0)
                       !dum3 = pfpz1(dume,mube(iname),pol0)
                       !dum4= dLa(iname)
                       !if (file_write ==1) then
                       !!'#pfpe, pfpz1, dL, dPz,dE, ith'
                       !     write(8000+proc_num,445) &
                       !     dum2,dum3,dum4,dpz0(l),dE0(l),name0(l)
                       !endif

                       !  real or complex omg ?
                       !fdum = (dume)**3*dE0(iname)/ekev*dpz0(iname)/&
                       !       e/e*dum4*tau_theta1(iname)*&
                       !       (omg*dum2-float(ntor)*dum3)
                       !fdum = (dume)**3*dE0(iname)/ekev*dpz0(iname)/&
                       !       e/e*dum4*tau_theta1(iname)*&
                       !       (real(omg)*dum2-float(ntor)*dum3)
!
                          do im = 1, nomode ! im-> m
                             do im1 = 1, nomode ! im1->k
                                do il = 0, nog ! il->l
                                   do il1 = 0, nog  ! il1->i
                                      do ip = 1, nopmode ! ip-> p for
                                                         ! matrix element M'_{i,l}^{k,m}
                                         sdum3 = fylmp(ip,iarray,im,il)
                                         sdum4 = conjg(ylmp(ip,iarray,im1,il1))
                                         !sdum3=(0.001,0.001)
                                         !sdum4=(0.002,0.002)
                                         if ((sdum3 /= cmplx(0.0,0.0)) .and. &
                                             (sdum4 /= cmplx(0.0,0.0))) then

                                            sdum5 = sdum3*sdum4
                                            temp_aint_KCT4D(il,il1,im,im1)=temp_aint_KCT4D(il,il1,im,im1)+&
                                                                           1.0/1.0*sdum5
                                                   
                                            !temp_aint_KCT4Dp(il,il1,im,im1)=temp_aint_KCT4Dp(il,il1,im,im1)+&
                                            !1.0/2.0*sdum4*fylmpp(ip,iarray,im,il)

                                            if (file_write ==1) then
                                             ! '# Re(YY*), Im(YY*), m,
                                             ! k, l, i, p '
                                             !write(7000+proc_num,447) &
                                             !sdum4,sdum3,mode(im),mode(im1),&
                                             !          il,il1,pmode(ip)
                                             endif
                                         !else
                                            !temp_aint_KCT4D(il,il1,im,im1)=temp_aint_KCT4D(il,il1,im,im1)+0.0
                                            !temp_aint_KCT4Dp(il,il1,im,im1)=temp_aint_KCT4Dp(il,il1,im,im1)+0.0

                                       endif  



                                      enddo

                                   enddo
                                enddo
                             enddo
                          enddo
!
                    endif
!
                endif
           enddo                                     


               
          call MPI_REDUCE(temp_aint_KCT4D, aint_KCT4D, size(aint_KCT4D),&
               MPI_DOUBLE_COMPLEX,MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         !call MPI_REDUCE(temp_aint_KCT4Dp, aint_KCT4Dp, size(aint_KCT4Dp),&
         !      MPI_DOUBLE_COMPLEX,MPI_SUM, 0, MPI_COMM_WORLD, ierr)

          call MPI_REDUCE(nprt_proc,nprt1,1,MPI_INTEGER,MPI_SUM,0, &
                        MPI_COMM_WORLD,ierr)
          call MPI_REDUCE(ik_proc,ik,1,MPI_INTEGER,MPI_SUM,0, &
                        MPI_COMM_WORLD,ierr)
         ! call MPI_REDUCE(Cn_proc,Cn1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
         !               MPI_COMM_WORLD,ierr)


          deallocate(temp_aint_KCT4D)
          !deallocate(temp_aint_KCT4Dp)
          deallocate(ylmp)
          deallocate(fylmp)
          !deallocate(fylmpp)
          deallocate(f1)
          deallocate(f2)
          deallocate(f3)
          deallocate(f4)
          deallocate(f5)
          deallocate(f6)
          deallocate(R0)
          deallocate(Z0)
          deallocate(XR0)
          deallocate(xx0)
          deallocate(theta0)
          deallocate(phi0)
          deallocate(en)
          deallocate(name0)
          deallocate(dE0)
          deallocate(ptcha)
          deallocate(ptch)
          deallocate(mube)
          deallocate(pz0)
          deallocate(sgnb)
          deallocate(dxa)
          deallocate(dLa)
          deallocate(dpz0)
          deallocate(w_theta)
          deallocate(w_theta1)
          deallocate(w_phi1)
          deallocate(tau_theta1)
          deallocate(psir_bar)


!



              if (proc_num == 0) then 
                 write(8,*) 'confined ions:', nprt1  ! 20230323
                 write(8,*) 'rows of Y6D_1.dat', ik
                 !write(8,*) 'Cn1:', Cn1
                 write(1001,78)
                 write(1001,*) nprt1,ik, vA_bar
                close(1001)
               endif
             enddo
           else
!
!!!!!!!!!!!! single particle name0 run !!!!!!!!!!!         
             j = 1   
!                                    
             vstart(1) = f1(name0(l),j) ! x
             vstart(2) = f2(name0(l),j) ! theta
             vstart(3) = f3(name0(l),j) ! phi: toroidal angle
             vstart(4) = f4(name0(l),j) ! v_||: v parallel
             vstart(5) = f5(name0(l),j) ! E: energy
             vstart(6) = f6(name0(l),j) ! mube: Lambda
             !print *, 'R',vstart(1)
             
             do i=1,nvar
                  v(i) = vstart(i)
             enddo
             ! initial positions record
             temp = b_nor(v(2),v(1))
             tmpfpol = inrcsval(gridx2,cscoefxpsi,v(1),0)
!                       particle energy with kev unit
             temp1 = (v(4)**2/2.0 + v(5)*v(6)*temp)/&
                         rhoh_bar/rhoh_bar*ekev
!                       particle toroidal canonical momentum  
             temp2 = tmpfpol*e**2 - v(4)/temp   ! P_phi 

             write(1,441) x, rf(v(2),v(1)), zf(v(2),v(1)), v(3), temp1, temp2
             write(2,441) x, v(1), v(2), v(1), qfun(v(1)), v(4)

             XR0(l) = rf(v(2),v(1)) - 1.0
             Z0(l) = zf(v(2),v(1))
             phi0(l) = v(3)
             ptcha(l) = ptch(l)
             x=x1  ! reset initial time for each ion fixed
             m=1
             i=0
             ik=0
             zzz=0.0
             p=0
             xw = 1.0
             !pw = inrcsval(gridx2,cscoefxpsi,xw,0)

             ! initial positions record for sp interp
             tdum(m)=x
             Rdum(m)=rf(v(2),v(1))
             Zdum(m)=zf(v(2),v(1))
             phidum(m)=v(3)
             thetadum(m)=v(2)
             psirdum(m)=v(1)

!                advance of particle with RK4 scheme 
             do k = 1,nstep
                    call derivs1(x,v,dv)
                    call rk4(v,dv,nvar,x,h,v)
                    !pol = cdbbsval2d(knotr,knotz,kr,kz,&
                    !      bscoefpsi,v(1),v(2),0,0)! particle flux psi 
!                   lost condition
                    if (v(1)>xw) then
                        write(3,443)  XR0(l), Z0(l), &
                                     phi0(l), ptcha(l),&
                                     temp1, temp2,name0(l)
                        exit  ! particle loss to last flux surface
                    endif
!                   particle data record
                     if (mod(k,modstep)==0) then
                        m=m+1
                        !temp = cdbbsval2d(knotr,knotz,kr,kz,&
                        !       bscoefscalb,v(1),v(2),0,0) ! B field normalise to B(0) at axis
                        temp = b_nor(v(2),v(1))
                        tmpfpol = inrcsval(gridx2,cscoefxpsi,v(1),0)
                        !tmpfpol = (10.0*log(v(1)**2 + 21.0/11.0))/11.0-&
                        !       (10.0*log(21.0/11.0))/11.0
                        
                        !tmpq = inrcsval(gridpsi1,cscoefq,pol,0) ! q
!                       particle energy with kev unit
                        temp1 = (v(4)**2/2.0 + v(5)*v(6)*temp)/&
                                rhoh_bar/rhoh_bar*ekev  
!                       particle toroidal canonical momentum  
                        !temp2 = tmpfpol*e**2 - v(4)/temp   ! P_phi 
                        !temp2 = v(4)/temp   ! P_phi

                        !v(5) = temp1/ekev*rhoh_bar**2
                        write(1,441) x, rf(v(2),v(1)), zf(v(2),v(1)), v(3), temp1, temp2
                        write(2,441) x, v(1), v(2), v(1), qfun(v(1)), v(4)
                        !write(69,*) b_nor(v(2),v(1)), b_star(v(2),v(1),v(4))
                                          
!
!                      poloidal frequency and toroidal frequency of
!                      particle
!                      as recording duration of cross mid plane 
                       tdum(m)=x
                       Rdum(m)=rf(v(2),v(1))
                       Zdum(m)=zf(v(2),v(1))  
                       phidum(m)=v(3)
                       thetadum(m)=v(2)
                       psirdum(m)=v(1)
                       !write(8,*) 'R=',Rdum(m)

                       if (m>1) then
                              if (zzz>0 .and. Zdum(m)<0) then   ! cross the mid plane by change sign of z
                                 i=i+1
                            ! interporate time linearly  for z=0
                                 tt(i) = ttt*(-Zdum(m))/(zzz-Zdum(m)) + &
                                         x*(-zzz)/(Zdum(m)-zzz)
                            ! interporate phi linearly for z=0
                                 pp(i) = ppp*(-Zdum(m))/(zzz-Zdum(m)) + &
                                         v(3)*(-zzz)/(Zdum(m)-zzz)
                            
                              end if
!
!
                       end if
                       !if (ik == 1) then
                       !   write(410,442) real(m), tdum(m),Rdum(m),Zdum(m),thetadum(m)
                       !end if

                       
                             ttt = x
                             xxx = rf(v(2),v(1))
                             zzz = zf(v(2),v(1))
                             ppp = v(3)
                             aaa = v(1)
                     end if

                     if (x+h==x) &
                          pause "stepsize vot signficant, in rkdumb"
                          x = x + h ! advance time
                     if (i==3) then
                           wtheta=2.0*pi/abs(tt(1)-tt(2)) ! omega_theta
                           wphi =abs( (pp(2)-pp(1))/(tt(2)-tt(1))) ! omega_phi
                           ttdum = abs(tt(1)-tt(2))
                           write(8,*) 'T_tran=',ttdum,'omega_t=',&
                                     wtheta,'omega_p=',wphi,&
                                     'R_min=',MinVal(Rdum(1:m-1))
!               phase position  and frequencies
                            !write(5,443) ttdum,temp1,temp2,vstart(6),wtheta,wphi,name0(l)
!
!               selected positions of particles
                          allocate(cscoefR(4,m-2),cscoefZ(4,m-2),cscoefphi(4,m-2),&
                                  cscoeftheta(4,m-2),cscoefpsir(4,m-2),cscoeftime(4,m-2))
                          call inrcsnak(tdum(1:m-1),Rdum(1:m-1),cscoefR)
                          call inrcsnak(tdum(1:m-1),Zdum(1:m-1),cscoefZ)
                          call inrcsnak(tdum(1:m-1),phidum(1:m-1),cscoefphi)
                          call inrcsnak(tdum(1:m-1),thetadum(1:m-1),cscoeftheta)
                          call inrcsnak(thetadum(1:m-1),tdum(1:m-1),cscoeftime)
                          call inrcsnak(tdum(1:m-1),psirdum(1:m-1),cscoefpsir)

!
!               choosen time to record
                          dum6 = 0.0
                          if (ptch(j) > 0) then
                              dum7 = TWOPI_D ! passing particles
                           else
                              dum7 = -TWOPI_D
                          endif
                          dum8 = inrcsval(thetadum(1:m-1),cscoeftime,dum6,0)
                          dum9 = inrcsval(thetadum(1:m-1),cscoeftime,dum7,0)
                          print *, 't0=',dum8,'t2pi=',dum9
                          ttdum = abs(dum9-dum8)
                          dtime = ttdum/real(nstep2-1,sp)
                          wtheta=TWOPI_D/ttdum
!                       test case
                          do p=1,nstep2
                             tdum1(p)=(p-1)*dtime+dum8
                        !     Rdum1(p)=inrcsval(tdum(1:m-1),cscoefR,tdum1(p),0) ! R
                        !     Zdum1(p)=inrcsval(tdum(1:m-1),cscoefZ,tdum1(p),0) ! Z
                        !     phidum1(p)=inrcsval(tdum(1:m-1),cscoefphi,tdum1(p),0) ! phi
                              thetadum1(p)=inrcsval(tdum(1:m-1),cscoeftheta,tdum1(p),0) ! theta
                        !     psirdum1(p)=inrcsval(tdum(1:m-1),cscoefpsir,tdum1(p),0) ! psir: r/a
                        !     write(4,443) tdum1(p), Rdum1(p), Zdum1(p),&
                        !                  phidum1(p), thetadum1(p), psirdum1(p), p
                        !     ifem = 5
                        !     mfem = 5
                        !     fi(p) = fem(ifem,mfem,psirdum1(p))
                        !     zsquar = cmplx(0.0, -thetadum1(p)+phidum1(p)+(wtheta-wphi)*tdum1(p))
                        !     gi(p) = fi(p)*exp(zsquar)
                              !fi(p) = sin(TWOPI_D/ttdum*tdum1(p))
                              fi(p) = cos(wtheta*tdum1(p))
                        !      write(1002,*) tdum1(p),thetadum1(p),wtheta*tdum1(p)
                          enddo
!
                         ! print *,'dt=',dtime
                         ! call simpc(nstep2,dtime,gi,aint0)
                         call simp(nstep2,dtime,fi,sdum)
                          print *,'sdum=',sdum
                          deallocate(cscoefR)
                          deallocate(cscoefZ)
                          deallocate(cscoefphi)
                          deallocate(cscoeftheta)
                          deallocate(cscoefpsir)
                          deallocate(cscoeftime)

!
                        exit ! one/two transit to exit
                     end if
             enddo
           endif
!        

    
           
       if (proc_num == 0) then 
          call CPU_TIME(ctime1)
          call system_clock(et2,count_rate)
          runt2 = (et2 - st2) * 1.0 / count_rate
          write(8,*) "orbit solve mod cpu time=", ctime1-ctime0,&
                     'for process number',proc_num,&
                     ' of ',numprocs, ' processes'
          write(8,*) "orbit solve mod wall time=", runt2
       endif

              call mpi_barrier(MPI_COMM_WORLD, ierr)
                !  if (ierr /= MPI_SUCCESS) then
                !   call MPI_Abort(MPI_COMM_WORLD, ierr)
                !  end if
             end subroutine main

        end module test2_mod
