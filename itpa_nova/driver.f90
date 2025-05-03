! $Id: driver.f90,v 1.8 2003/12/22 13:29:55 pletzer Exp $
! 2024/03/23 Limin Yu, particle mesh
program driver
  use initial_mod
  implicit none
  integer input

  print *, 'GACO initial driver'
  call ini
  !print *, '0: all test runs'
  !print *, '1: ktae'
  !print *, '2: eigenbessel'
  !print *, '3: eigenwave'
  !print *, '4: singular'
  !print *, '5: Airy'
  !print *, '6: order4'
  !print *, '7: coupled'
  !print *, 'select test case now '
  !read(5,*) input
  
  !select case(input)

  !   case (1)
  !      call ktae_test_run
  !   case (2)
  !      call eigenbessel_test_run
  !   case (3)
  !      call eigenwave_test_run
  !   case (4)
  !      call singular_test_run
  !   case (5)
  !      call Airy_test_run
  !   case (6)
  !      call order4_test_run
  !   case (7)
  !      call coupled_test_run
  !   case default
  !      call ktae_test_run
  !      call eigenbessel_test_run
  !      call eigenwave_test_run
  !      call singular_test_run
  !      call Airy_test_run
  !      call order4_test_run
  !      call coupled_test_run        

  !end select

end program driver

!========================================================================
subroutine ktae_test_run
  use ktae
  implicit none
  type(ktae_obj) :: this
  call ktae_test(this)
end subroutine ktae_test_run

!========================================================================
!========================================================================
subroutine eigenbessel_test_run
  use eigenbessel
  implicit none
  type(eigenbessel_obj) :: this
  call eigenbessel_test(this)
end subroutine eigenbessel_test_run

!========================================================================
!========================================================================
subroutine eigenwave_test_run
  use eigenwave
  implicit none
  type(eigenwave_obj) :: this
  call eigenwave_test(this)
end subroutine eigenwave_test_run

!========================================================================
!========================================================================
subroutine singular_test_run
  use singular
  implicit none
  type(singular_obj) :: this
  call singular_test(this)
end subroutine singular_test_run

!========================================================================
!========================================================================
subroutine order4_test_run
  use order4
  implicit none
  type(order4_obj) :: this
  call order4_test(this)
end subroutine order4_test_run

!========================================================================
subroutine coupled_test_run
  use coupled
  implicit none
  type(coupled_obj) :: this
  call coupled_test(this)
end subroutine coupled_test_run
!========================================================================

subroutine Airy_test_run
  use airy
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  type(airy_obj) :: this

    integer :: Nz, Nj
    real(r8) sigma_dx2, sampling, xmin, xmax, tol, alpha, pi
    complex(r8) :: bc0, bc1
    real(r8), allocatable :: x(:)
    complex(r8), allocatable :: y(:), y_exact(:)
    real(r8) :: error

    integer n, i

    Nz = 20
    Nj = 5
    sigma_dx2 = 5 ! 2
    sampling = 2._r8    
    call airy_init(this, Nz, Nj, sigma_dx2, sampling)

    xmin = 0
    xmax = 1
    call airy_set_domain(this, xmin, xmax)

    tol = 1.e-8_r8
    call airy_set_tol(this, tol)

    pi = acos(-1._r8)
    alpha = 100*pi ! 11*pi/2._r8 ! 0
    call airy_set_equation(this, alpha)

    ! set the boundary conditions
    call gsl_sf_airy_ai((alpha/2._r8)**(2._r8/3._r8)*(2._r8*0._r8-1._r8), bc0)
    call gsl_sf_airy_ai((alpha/2._r8)**(2._r8/3._r8)*(2._r8*1._r8-1._r8), bc1)
    call airy_set_bcs(this, bc0, bc1)

    call airy_assemble(this)
    call airy_solve(this)

    n = this % iter2 % ntot * 5
    allocate(x(n), y(n), y_exact(n))
    do i = 1, n
       x(i) = xmin + (xmax-xmin)*real(i-1,r8)/real(n-1, r8)
       call gsl_sf_airy_ai((alpha/2._r8)**(2._r8/3._r8)*(2._r8*x(i)-1._r8), y_exact(i))
    end do
    call airy_get_solution(this, x, y)
    error = sqrt(sum( abs(y-y_exact)**2 )/n )
    print *,'alpha=', alpha, ' Nz=', Nz, ' Nj=', Nj, 'error = ', error
    call airy_plot(this, x, y)
    
    deallocate(x, y, y_exact)

    call airy_free(this)

  end subroutine Airy_test_run


!========================================================================

subroutine CurlCurl

  ! Example of test problem

  use indexer
  use operators
  use supralu_mod
  use gabor

  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  real(r8), parameter :: twopi = 6.2831853071795862320_r8
  
  type(indexer_obj) :: iter1, iter2
  type(sparseC16_obj) :: amat

  integer ier, all_inda(4), all_indb(7), ia, ib, i1, i2
  integer :: rank, nnz, base

  ! Grid sizes

  ! Gabor positions (range from 0:Nz)
  integer, parameter :: Nz(3) = (/2, 2, 2/)
  ! k-vector grid (Fourier modes range from -Nj:+Nj)
  integer, parameter :: Nj(3) = (/0, 0, 0/)
  ! No of components
  integer, parameter :: Nc = 1
  ! Collocation
  integer, parameter :: Nx(3) = (2*Nj+1)*(Nz+1)-1

  real(r8) :: kwav(3) ! wave-vector
  real(r8) :: zpos(3) ! envelop position
  real(r8) :: xpos(3) ! collocation position
  real(r8) :: delta_k(3), delta_x(3), delta_z(3)
  real(r8) :: sampling(3) ! sampling rate
  real(r8) :: sigma_dz2(3), sigma(3) ! inverse extent of envelopes
  real(r8) :: xmin(3), xmax(3), xmin_eps(3), xmax_eps(3)
  real(r8) :: tol, minus_log_tol ! sparsity tolerance
  complex(r8) :: val, gabor_fct
  complex(r8), allocatable :: x(:)
  logical newcol_flag
  real(r8) :: eps

  integer :: sizes2(7), start2(7)
  integer :: sizes1(4), start1(4)

  ! INPUTS

  ! Box boundaries
  xmin = 0
  xmax = 1

  ! Sampling rate
  sampling = (/1.1_r8, 1.5_r8, 2._r8/)

  ! Gabor inverse extent normalized to delta_z**2. delta_z is the Gabor separation
  sigma_dz2 = 1.5_r8

  delta_z = (xmax-xmin)/Nz
  delta_x = (xmax-xmin)/Nx
  delta_k = twopi/(delta_z * sampling)

  sigma = sigma_dz2 / delta_z**2

  tol = 1.e-10_r8
  minus_log_tol = -log(tol)
  eps = 10*epsilon(1._r8)

  ! ITERATOR SETUP
  
  ! Looping order: space components, then Fourier and finally 
  ! envelope positions
  sizes2 = (/Nc, 2*Nj(1)+1, 2*Nj(2)+1, 2*Nj(3)+1, Nz(1)+1, Nz(2)+1, Nz(3)+1/)
  start2 = (/ 1,  -Nj(1)  ,  -Nj(2)  ,  -Nj(3)  ,       0,       0,       0/)
  call indexer_init(iter2, sizes2, start2, ier)
  ! collocation looping:
  sizes1 = (/Nc, Nx(1)+1, Nx(2)+1, Nx(3)+1/)
  start1 = (/ 1,       0,       0,       0/)
  call indexer_init(iter1, sizes1, start1, ier)

  ! ASSEMBLY

  rank = iter1%ntot
  nnz = rank**2 !max(100000, int(1.0 * rank**2))
  base = 0
  call sparseC16_init(amat, rank, nnz, base, ier)
  call sparseC16_error(amat, 6, ier)

  ! right hand side vector
  allocate(x(iter1%ntot))
  x = 0  

  xmin_eps = xmin + eps
  xmax_eps = xmax - eps

  print *,'iter1%ntot=', iter1%ntot,' iter2%ntot=', iter2%ntot

  do i2 = 1, iter2%ntot
     ! Iterate over Gabors
     newcol_flag = .TRUE.
     call indexer_next(iter2, ier)
     call indexer_get(iter2, all_indb, ier)
     ib = all_indb(1)
     kwav = delta_k * all_indb(2:4)
     zpos = xmin + delta_z * all_indb(5:7)
     print *,'ib=',ib,' kwav=',kwav,' zpos=', zpos
     
     do i1 = 1, iter1%ntot
        ! Iterate over collocation points
        call indexer_next(iter1, ier)
        call indexer_get(iter1, all_inda, ier)
        ia = all_inda(1)
        xpos = xmin + delta_x * all_inda(2:4)
        gabor_fct = GaborN(xpos, zpos, kwav, sigma)

        if( all(xpos <= xmin_eps) ) then
           val = gabor_fct
           call sparseC16_set_next(amat, i1-1, newcol_flag, val, ier)
           call sparseC16_error(amat, 6, ier)
           newcol_flag = .FALSE.
           x(i1) = 0 ! BC at x=xmin
           cycle
        end if
        if( all(xpos >= xmax_eps) ) then
           val = gabor_fct
           call sparseC16_set_next(amat, i1-1, newcol_flag, val, ier)
           call sparseC16_error(amat, 6, ier)
           newcol_flag = .FALSE.
           if(ia==1) x(i1) = 1 ! BC at x=xmax
           cycle
        end if

        if (sum(0.5_r8*sigma*(xpos-zpos)**2) > minus_log_tol) cycle

        call curl_curl(ia, ib, xpos, kwav, zpos, sigma, val, ier)
        val = val * gabor_fct
        print *,'ia=',ia,' ib=',ib,' val=', val
        call sparseC16_set_next(amat, i1-1, newcol_flag, val, ier)
        call sparseC16_error(amat, 6, ier)
        newcol_flag = .FALSE.
     end do
  end do

  ! SOLVE
  
  call sparseC16_save(amat, 'amat.nc', ier)
  call sparseC16_error(amat, 6, ier)
  call sparseC16_solve(amat, x, ier)
  call sparseC16_error(amat, 6, ier)

  call indexer_free(iter1, ier)
  call indexer_free(iter2, ier)
  call sparseC16_free(amat, ier)
  deallocate(x)

end subroutine CurlCurl
