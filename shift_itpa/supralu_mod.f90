module supralu_mod
  use vars_k
  use vars_e
  use vars
  use inrtype
  use test2_mod
  use shared_mod
  ! Column-compressed-storage based sparse matrix container

  implicit none
  integer, parameter, private :: r8 = selected_real_kind(12,100)
  character(*), parameter, private :: version = &
       & '$Id: supralu_mod.f90,v 1.4 2003/12/22 17:24:29 pletzer Exp $'
  
  type sparseC16_obj
     integer, pointer :: irow(:)
     integer, pointer :: jcol_ptr(:)
     complex(r8), pointer :: values(:)
     integer :: size, size_ptr, nnz, col_counter
     integer :: base

     ! SuperLU stuff

     ! column permutation specification
     !   permc_spec = 0: natural ordering 
     !   permc_spec = 1: minimum degree on structure of A'*A
     !   permc_spec = 2: minimum degree on structure of A'+A
     !   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     integer :: permc_spec 
     
  end type sparseC16_obj

  type sparseR8_obj
     integer, pointer :: irow(:)
     integer, pointer :: jcol_ptr(:)
     real(r8), pointer :: values(:)
     integer :: size, size_ptr, nnz, col_counter
     integer :: base

     ! SuperLU stuff

     ! column permutation specification
     !   permc_spec = 0: natural ordering 
     !   permc_spec = 1: minimum degree on structure of A'*A
     !   permc_spec = 2: minimum degree on structure of A'+A
     !   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     integer :: permc_spec 
     
  end type sparseR8_obj

  contains
!===================================================================================

! complex(r8) ::

    subroutine sparseC16_init(this, rank, size, base, ier)

      ! Constructor.
      !
      ! rank: matrix rank (or approximation if not known)
      ! size: >= no of non-zero values to store. If 
      ! size turns out to be too tight, fresh arrays will be 
      ! reallocated and the data will be copied, however this may
      ! can cause performance degradation.
      !
      ! base: indexing base (0 like C or 1 like Fortran)

      implicit none
      type(sparseC16_obj) :: this
      integer, intent(in) :: rank
      integer, intent(in) :: size
      integer, intent(in) :: base
      integer, intent(out) :: ier

      integer jer

      ier = 0

      this % permc_spec = 2 ! seems to work best

      this % nnz = 0
      this % size_ptr = rank+1
      this % size = size
      this % base = base
      this % col_counter = 0
      allocate(this%irow(this%size), stat=jer)
      if(jer /= 0 ) then
         ier = 9
         return
      end if
      allocate(this%values(this%size), stat=jer)
      if(jer /= 0 ) then
         ier = 10
         return
      end if
      allocate(this%jcol_ptr(this%size_ptr), stat=jer)
      if(jer /= 0 ) then
         ier = 11
         return
      end if
      this%irow = 0
      this%jcol_ptr = 0
      this%values = 0
      if(jer /=0 ) ier = 1
      
    end subroutine sparseC16_init
    
    subroutine sparseC16_free(this, ier)

      ! Destructor

      implicit none
      type(sparseC16_obj) :: this
      integer, intent(out) :: ier

      integer jer

      ier = 0
      deallocate(this%irow, stat=jer)
      deallocate(this%values, stat=jer)
      deallocate(this%jcol_ptr, stat=jer)
      if(jer/=0) ier = 2
      
    end subroutine sparseC16_free

    subroutine sparseC16_copy(this, that, ier)

      ! Copy this into that. that should be unitialized.

      implicit none
      type(sparseC16_obj) :: this, that
      integer, intent(out) :: ier

      integer jer

      ier = 0
      that % nnz = this % nnz
      that % size_ptr = this % size_ptr
      that % size = this % size
      that % base = this % base
      that % col_counter = this % col_counter
      allocate(that%irow(that%size), stat=jer)
      if(jer /= 0 ) then
         ier = 9
         return
      end if
      allocate(that%values(that%size), stat=jer)
      if(jer /= 0 ) then
         ier = 10
         return
      end if
      allocate(that%jcol_ptr(that%size_ptr), stat=jer)
      if(jer /= 0 ) then
         ier = 11
         return
      end if
      that % irow = this%irow
      that % jcol_ptr = this%jcol_ptr
      that % values = this%values
      if(jer /=0 ) ier = 12
      
      
    end subroutine sparseC16_copy

    subroutine sparseC16_compare_shapes(this, that, ier)

      ! Compare two sparse matrices

      implicit none
      type(sparseC16_obj) :: this, that
      integer, intent(out) :: ier

      integer i

      ier = 0
      if(this%base /= that%base) then
         ier = 13
         return
      end if
      if(this%nnz /= that%nnz) then
         ier = 14
         return
      end if
      if(this%size_ptr /= that%size_ptr) then
         ier = 15
         return
      end if
      if(this%size /= that%size) then
         ier = 16
         return
      end if
      if(this%col_counter /= that%col_counter) then
         ier = 17
         return
      end if
      do i = 1, this%size_ptr
         if(this%jcol_ptr(i) /= that%jcol_ptr(i)) then
            ier = 18
            return 
         end if
      enddo
      do i = 1, this%size
         if(this%irow(i) /= that%irow(i)) then
            ier = 19
            return      
         end if
      enddo


    end subroutine sparseC16_compare_shapes

    subroutine sparseC16_set_next(this, row, newcol_flag, value, ier)

      ! Set next value. Must iterate down each column 
      ! running from left to right.
      !
    
      implicit none
      type(sparseC16_obj) :: this
      integer, intent(in) :: row 
      logical, intent(in) :: newcol_flag ! .TRUE. if new column, .FALSE. otherwise
      complex(r8), intent(in) :: value 
      integer, intent(out) :: ier

      ier = 0
      this % nnz = this % nnz + 1
      if(this % nnz > this % size) then         
         ier = 3
         return
      endif
      this % irow(this%nnz) = row
      this % values(this%nnz) = value
      if(newcol_flag) then
         this % col_counter = this % col_counter + 1
         this % jcol_ptr(this % col_counter) = this%nnz - 1 + this%base
      endif
      this % jcol_ptr(this % col_counter+1) = this%nnz + this%base
      
    end subroutine sparseC16_set_next

    subroutine sparseC16_get_size(this, nrows, ncols, ier)

      ! Return no of rows and columns

      implicit none
      type(sparseC16_obj) :: this
      integer, intent(out) :: nrows, ncols
      integer, intent(out) :: ier

      ier = 0
      nrows = maxval(this % irow)
      ncols = this % col_counter
    end subroutine sparseC16_get_size

    subroutine sparseC16_solve(this, rhs, ier)

      ! Solve matrix system
      
      implicit none
      type(sparseC16_obj) :: this
      complex(r8), intent(inout) :: rhs(:) ! right-hand side vector, solution vector upon return
      integer, intent(out) :: ier

      integer n, jer
      integer :: handle(8) ! opaque handle

      ier = 0
      n = size(rhs)

      ! note: we pass here address (1) of vectors to force passing the 
      ! arguments by reference. 
      call zsupralu_new(handle, &
           & this%values(1), this%irow(1), this%jcol_ptr(1), this%nnz, n, &
           & jer)
      call zsupralu_colperm(handle, this % permc_spec, jer)
      call zsupralu_lu(handle, jer)
      call zsupralu_solve(handle, rhs(1), jer)
      call zsupralu_del(handle, jer)

      if(jer/=0) ier = 8

    end subroutine sparseC16_solve

    subroutine sparseC16_det(this, res_mantissa, res_exponent, ier)

      ! Return the determinant of A in the form of 
      !
      ! res_mantissa * 2**res_exponent

      implicit none
      type(sparseC16_obj) :: this
      complex(r8), intent(out) :: res_mantissa
      integer, intent(out) :: res_exponent
      integer, intent(out) :: ier
      
      integer :: handle(8) ! opaque handle
      integer :: n, jer

      ier = 0
      n = this % col_counter
      call zsupralu_new(handle, &
           & this%values(1), this%irow(1), this%jcol_ptr(1), this%nnz, n, &
           & jer)
      call zsupralu_colperm(handle, this % permc_spec, jer)
      call zsupralu_lu(handle, jer)
      call zsupralu_determinant(handle, res_mantissa, res_exponent, jer)
      call zsupralu_del(handle, jer)
      if(jer/=0) ier = 20

    end subroutine sparseC16_det

   subroutine sparseC16_eigen(Amat, Bmat, x, lambda, tol, nmax, ier)

      ! Solve A * x = lambda B * x by inverse iterations

      implicit none
      type(sparseC16_obj) :: Amat, Bmat ! Must have the same sparsity pattern!
      complex(r8), intent(inout) :: lambda ! initial guess/computed eigenvalue
      complex(r8), intent(inout) :: x(:) ! initial egenvector guess/computed eigenvector
      real(r8), intent(inout) :: tol ! tolerance
      integer, intent(inout) :: nmax ! max number of iterations
      integer, intent(out) :: ier

      integer n, i, i1, i2, i3
      complex(r8), allocatable :: rhs(:)
      real(r8) :: residue, norm
      complex(r8) :: newv_b_v, newv_b_newv
      type(sparseC16_obj) :: Cmat

      ier = 0
      n = size(x)
      allocate(rhs(n))
      call sparseC16_copy(Amat, Cmat, ier)
      !nmax = 1
      do i = 1, nmax
         Cmat % values = Amat % values - lambda * Bmat % values
         call sparseC16_A_dot_x(Cmat, x, rhs, ier)
         residue = sqrt( dot_product(rhs, rhs) )
         if(residue < tol) exit

         call sparseC16_A_dot_x(Bmat, x, rhs, ier)
         x = rhs
         call sparseC16_save(Cmat, 'cmat.nc', ier)
         call sparseC16_solve(Cmat, x, ier)
         newv_b_v = dot_product(x, rhs)
         call sparseC16_A_dot_x(Bmat, x, rhs, ier)
         newv_b_newv = dot_product(x, rhs)
         lambda = lambda + newv_b_v/newv_b_newv
         print *,'iteration=', i,' lambda=', lambda, ' residue=', residue
         norm = sqrt(dot_product(x, x))
         x = x/norm
      end do
      write(3100,*) lambda
      close(3100)
!
      i1 = 101
      i2 = 1
      i3 = 0
      open(unit=1005,file='fort.1005',status='replace')
      write(1005,*) nog,i2
      write(1005,*) lambda
      write(1005,*) i3,real(rho),real(rho),imag(rho)
      close(1005)
!
      res = residue
      nmax = i
      tol = residue

      call sparseC16_free(Cmat, ier)
      deallocate(rhs)

    end subroutine sparseC16_eigen

    subroutine sparseC16_eigen_k(Amat, B0mat, x, lambda,matsize,rank, base,tol, nmax, ier)

      ! Solve A(lambda) * x = lambda B * x by inverse iterations

      implicit none
      type(sparseC16_obj) :: Amat, Bmat ! Must have the same sparsity pattern! 
      complex(r8), intent(inout) :: lambda ! initial guess/computed eigenvalue
      complex(r8), intent(inout) :: x(:) ! initial egenvector guess/computed eigenvector
      real(r8), intent(inout) :: tol ! tolerance
      integer, intent(inout) :: nmax ! max number of iterations      
      integer, intent(out) :: ier
      integer ngrid,i1,i2,ii,jj
      integer n, i, j, k, ii1, ii2, ii3
      integer matsize,rank,base
      complex(r8), allocatable :: rhs(:)
      real(r8) :: residue, norm
      complex(r8) :: newv_b_v, newv_b_newv
      !complex(r8) :: rho
      type(sparseC16_obj) :: Cmat,B0mat,Dmat

      ier = 0
      n = size(x)
      if (allocated(rhs)) then
         deallocate(rhs)
      endif
      allocate(rhs(n))
      call sparseC16_copy(Amat, Cmat, ier)
      call sparseC16_copy(B0mat, Bmat, ier)      
     do k = 1, nmax
         if (k>1) then
           ! call main
            call matrix
            ngrid=nomode*(nog-1)-1
            call sparseC16_init(Amat, rank, matsize, base, ier)
            call sparseC16_init(Dmat, rank, matsize, base, ier)
            do j=0,ngrid
              i1 = 0
              i2 = ngrid
              jj=j+1
              do i=i1,i2
                 ii=i+1
                 if(i.eq.i1)then
                    call sparseC16_set_next(Amat, i, .TRUE., aa(ii,jj), ier)
                    call sparseC16_set_next(Dmat, i, .TRUE., dd(ii,jj), ier)
                 else
                    call sparseC16_set_next(Amat, i, .FALSE., aa(ii,jj), ier)
                    call sparseC16_set_next(Dmat, i, .FALSE., dd(ii,jj), ier)
                 endif
              enddo
            enddo      
           Bmat % values = B0mat % values - Dmat % values
         endif 
         Cmat % values = Amat % values - lambda * B0mat % values
         call sparseC16_A_dot_x(Cmat, x, rhs, ier)
         residue = sqrt( dot_product(rhs, rhs) )
         if(residue < tol) exit
         
         call sparseC16_A_dot_x(Bmat, x, rhs, ier)
         x = rhs
         call sparseC16_save(Cmat, 'cmat.nc', ier)
         call sparseC16_solve(Cmat, x, ier)
         newv_b_v = dot_product(x, rhs)
         call sparseC16_A_dot_x(Bmat, x, rhs, ier)
         newv_b_newv = dot_product(x, rhs)
         lambda = lambda + newv_b_v/newv_b_newv
         print *,'iteration=', k,' lambda=', lambda, ' residue=', residue
         
         norm = sqrt(dot_product(x, x))
         x = x/norm
     end do
     write(3100,*) lambda
     close(3100)
     !
      ii1 = 101
      ii2 = 1
      ii3 = 0
      open(unit=1005,file='fort.1005',status='replace')
      write(1005,*) nog,ii2
      write(1005,*) lambda
      write(1005,*) ii3,ii3,ii3,ii3
      close(1005)
!
      nmax = k
      tol = residue

      call sparseC16_free(Cmat, ier)
      deallocate(rhs)

    end subroutine sparseC16_eigen_k

    subroutine sparseC16_A_dot_x(this, x, res, ier)

      ! matrix . x => res

      implicit none
      type(sparseC16_obj) :: this
      complex(r8), intent(in) :: x(:) ! vector
      complex(r8), intent(out) :: res(:) ! result
      integer, intent(out) :: ier

      integer i, k, rank, j

      ier = 0
      res = 0
      rank = maxval(this % irow) + 1 - this%base

      
      do j = this%base+1, rank
         do k = this%jcol_ptr(j)+1-this%base, this%jcol_ptr(j+1)-this%base
            i = this%irow(k) + 1 - this%base
            res(i) = res(i) + this%values(k)*x(j)
         enddo
      enddo
      
    end subroutine sparseC16_A_dot_x

    subroutine sparseC16_save(this, filename, ier)

    ! Save state in netCDF file 'filename'

    use ezcdf
    implicit none
    type(sparseC16_obj) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer ncid, jer
    ier = 0

    call cdf_open(ncid, filename, 'w', jer)

    call cdf_define(ncid, 'version', (/len_trim(version),0,0/), 'CHAR', jer)
    call cdf_define(ncid, 'size', this%size, jer)
    call cdf_define(ncid, 'size_ptr', this%size_ptr, jer)
    call cdf_define(ncid, 'nnz', this%nnz, jer)
    call cdf_define(ncid, 'col_counter', this%col_counter, jer)
    call cdf_define(ncid, 'base', this%base, jer)
    call cdf_define(ncid, 'irow', this%irow(1:this%nnz), jer)
    call cdf_define(ncid, 'jcol_ptr', this%jcol_ptr, jer)
    call cdf_define(ncid, 'values', this%values(1:this%nnz), jer)

    call cdf_write(ncid, 'version', version, jer)
    call cdf_write(ncid, 'size', this%size, jer)
    call cdf_write(ncid, 'size_ptr', this%size_ptr, jer)
    call cdf_write(ncid, 'nnz', this%nnz, jer)
    call cdf_write(ncid, 'col_counter', this%col_counter, jer)
    call cdf_write(ncid, 'base', this%base, jer)
    call cdf_write(ncid, 'irow', this%irow(1:this%nnz), jer)
    call cdf_write(ncid, 'jcol_ptr', this%jcol_ptr, jer)
    call cdf_write(ncid, 'values', this%values(1:this%nnz), jer)

    call cdf_close(ncid, jer) 
    if(jer /= 0) ier = 4

    end subroutine sparseC16_save

    subroutine sparseC16_load(this, filename, ier)

      ! Load state from netCDF file 'filename'

    use ezcdf
    implicit none
    type(sparseC16_obj) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer ncid, jer
    ier = 0

    call cdf_open(ncid, filename, 'r', jer)

    call cdf_read(ncid, 'size', this%size, jer)
    call cdf_read(ncid, 'size_ptr', this%size_ptr, jer)
    call cdf_read(ncid, 'nnz', this%nnz, jer)
    call cdf_read(ncid, 'col_counter', this%col_counter, jer)
    call cdf_read(ncid, 'base', this%base, jer)
    allocate(this % irow(this%size), stat=jer)
    allocate(this % jcol_ptr(this%size_ptr), stat=jer)
    allocate(this % values(this%size), stat=jer)
    if(jer /= 0) then 
       ier = 1
       return
    endif
    call cdf_read(ncid, 'irow', this%irow, jer)
    call cdf_read(ncid, 'jcol_ptr', this%jcol_ptr, jer)
    call cdf_read(ncid, 'values', this%values, jer)

    call cdf_close(ncid, jer) 
    if(jer /= 0) ier = 5
      
    end subroutine sparseC16_load


    subroutine sparseC16_test(this, ier)

      ! Test unit

      implicit none
      type(sparseC16_obj) :: this, that, amat, bmat, dmat
      integer, intent(out) :: ier

      integer rank, size, base, i, nmax,ngrid,i1,i2,j
      !integer nomode,ii,jj
      integer ii,jj
      integer n,n2,noiter,iter,icase
      parameter(n=800,n2=6*n)
      real(r8) :: tol,abig,bsm,yy,rhor,rhor1,rhor2,drho,rhoi
      real(r8), allocatable :: yr(:)
      complex(r8), allocatable :: rhs(:), x(:), res(:), y(:)
      complex(r8) :: lambda
!      complex(r8) aa(0:1000)
      !complex(r8) :: rho
      real(r8) :: bbm,bbp,eps1,eps,xx,am,h,alamr,xxm,xxp,yya
      real(r8) :: lambr,lambr1,lambr2,dlamb,lambi
      real(r8) :: xau,fau1,fau2,fau3,eta1,eta2,deta
      real(r8) :: q01,q02,dq,qp1,qp2,dqp,fau1a
      complex(r8) :: fau4,zin,zout
      real(r8) :: x1,x2,dxx,xy,xend1
      integer :: nxx,idon
      complex(r8) :: gg1,gg2,hh1,hh2

      ikct = 0 ! 0: simple iteration for fixed lambda in Amat 
               ! 1: newton iteration for Amat as function of lambda      
      if (ikct.eq.1) then
        allocate(aa(n2,n2),bb(n2,n2),dd(n2,n2))
      else
        allocate(aa(n2,n2),bb(n2,n2)) 
      endif
!
      !call profiledat
      idon=0
!test
      if(idon.eq.1)then
      write(6,*)"input nxx,x1,x2"
      read(5,*)nxx,x1,x2
      dxx=(x2-x1)/float(nxx)
      do i=1,nxx
      xy=x1+dxx*float(i-1)
      gg1=gfun(xy,1)
      gg2=gfun(xy,2)
      hh1=gfun1(xy,1)
      hh2=gfun1(xy,2)
!      write(6,600)xy,gg1,hh1,gg2,hh2
 600  format(9(2x,1e8.2))
      enddo

      write(6,*)"zout=",zout
      stop
      endif

      nog=50
      write(0,*)"input nog,icase"
      !read(5,*)nog,icase
      read(1005,*) nog,icase
!test

!      do i=1,201
      do i=1,401
      xau=0.0025*float(i-1)
      fau1=den(xau)
      fau1a=den1(xau)
      fau2=tfun(xau)
      fau3=qfun(xau)
      fau4=gfun(xau,1)
!      write(0,113)xau,fau1,fau2,fau3,real(fau4),imag(fau4)
      write(51,114)xau,fau1,fau1a
 113  format(6(2x,1e12.6))
 114  format(3(2x,1e12.6))
      enddo
!      stop

! define grid points

       !call setgrid
       


      ier = 0
      print *, 'testing sparseC16...', version

      base = 0 ! C-like
!      ngrid=100
      !nomode=10
      ngrid=nomode*(nog-1)-1
      !print *,'nomode=',nomode,'nog=',nog,'ngrid=',ngrid
!      read(23,*)ngrid
      abig=1.0e+12_r8
!      h=2.0D0*xend/float(ngrid+2)
      xend1=xend-xini
      h=2.0D0*xend1/float(ngrid+2)
      bsm=h**2
      rank = ngrid+1
      !size = nomode*7*rank-12*nomode
      size = rank **2
      ! matrix to store
      ! [ 0 2 x 7 x ]
      ! [ 1 3 x x x ]
      ! [ x x 4 x x ]
      ! [ x x 5 8 x ]
      ! [ x x 6 x 9 ]
      ! the result should be 
      ! irow = [0,1,0,1,2,3,4,0,3,4]
      ! jcol_ptr = [0,2,4,7,9]
      ! values = [0,1,2,3,4,5,6,7,8,9]
      print *, '...assembling'

      am=1.0
      eps=0.2
      if (allocated(y)) then
         deallocate(y)
      endif
      if (allocated(yr)) then
         deallocate(yr)
      endif
      allocate(y(ngrid+1))
      allocate(yr(ngrid+1))

      noiter=2


      lambda=cmplx(0.138,0.0D0)
      rhor1=0.003D0
      rhor2=0.003D0
      rhor=0.003D0
      rhoi=0.0D0
!      etai=0.003D0
!      icase=5

      if(icase.eq.1)then
      write(0,*)"input inital guess for lambda"
      !read(5,*)lambda
      read(1005,*) lambda
      write(0,*)"input noiter,rhor1,rhor2,rhoi"
      !read(5,*)noiter,rhor1,rhor2,rhoi
      read(1005,*) noiter,rhor1,rhor2,rhoi
      close(1005)
      print*,'nog',nog,'lambda=',lambda
      drho=(rhor2-rhor1)/float(noiter+1)
      else if(icase.eq.2)then
      write(0,*)"input rhor,rhoi"
      !read(5,*)rhor,rhoi
      read(1005,*) rhor,rhoi
      write(0,*)"input noiter,lambr1,lambr2,lambi"
      !read(5,*)noiter,lambr1,lambr2,lambi
      read(1005,*)noiter,lambr1,lambr2,lambi
      close(1005)
      dlamb=(lambr2-lambr1)/float(noiter)
      else if(icase.eq.3)then
      read(5,*)lambda
      write(0,*)"input noiter,eta1,eta2"
      read(5,*)noiter,eta1,eta2
      deta=(eta2-eta1)/float(noiter)
      else if(icase.eq.4)then
      read(5,*)lambda
      write(0,*)"input noiter,q01,q02"
      read(5,*)noiter,q01,q02
      dq=(q02-q01)/float(noiter)
      else if(icase.eq.5)then
      read(5,*)lambda
      write(0,*)"input noiter,qp1,qp2"
      read(5,*)noiter,qp1,qp2
      dqp=(qp2-qp1)/float(noiter)
      endif


     !call set_couple_coefficients
      write(20,*)noiter+1
!      write(19,*)nomode,noiter,nog
      !rho= (0.0,0.0)
      !call matrix(rho)

      do iter=1,noiter+1
         if(icase.eq.1)then
         rhor=rhor1+drho*float(iter-1)
      rho=rhor*cmplx(1.0,rhoi)
      else if(icase.eq.2)then
        rho=rhor*cmplx(1.0,rhoi)
         lambr=lambr1+dlamb*float(iter-1)
      lambda=cmplx(lambr,lambi)
      else if(icase.eq.3)then
        rho=rhor*cmplx(1.0,rhoi)
      etai=eta1+deta*float(iter-1)
      else if(icase.eq.4)then
        rho=rhor*cmplx(1.0,rhoi)
      q0=q01+dq*float(iter-1)
      else if(icase.eq.5)then
        rho=rhor*cmplx(1.0,rhoi)
      cq1=qp1+dqp*float(iter-1)
      endif
      
     ! call main ! test
      call matrix
     
      call sparseC16_init(that, rank, size, base, ier)

      do j=0,ngrid
      !i1=(j/nomode)*nomode-3*nomode
      !i2=(j/nomode)*nomode+4*nomode-1
      i1 = 0
      i2 = ngrid
      !if(i1.lt.0)i1=0
      !if(i2.gt.ngrid)i2=ngrid

!      read(23,*)(aa(i),i=i1,i2)
!      write(0,*)(aa(i),i=i1,i2)
      jj=j+1
      do i=i1,i2
      ii=i+1
      if(i.eq.i1)then
      call sparseC16_set_next(that, i, .TRUE., bb(ii,jj), ier)
      else
      call sparseC16_set_next(that, i, .FALSE.,bb(ii,jj), ier)
      endif
      enddo
!      write(0,*)"j,bbb",j
      enddo

      call sparseC16_init(this, rank, size, base, ier)

      do j=0,ngrid

      !i1=(j/nomode)*nomode+1-3*nomode-1
      !i2=(j/nomode)*nomode+4*nomode-1
      !if(i1.lt.0)i1=0
      !if(i2.gt.ngrid)i2=ngrid
      i1 = 0
      i2 = ngrid

!c      read(23,*)(aa(i),i=i1,i2)
      jj=j+1
      do i=i1,i2
      ii=i+1
      if(i.eq.i1)then
      call sparseC16_set_next(this, i, .TRUE., aa(ii,jj), ier)
      else
      call sparseC16_set_next(this, i, .FALSE.,aa(ii,jj), ier)
      endif
      enddo
!      write(0,*)"j,aa",j
      enddo


      tol = 1.e-12_r8
      nmax = 10
!
!      write(0,*)"input lambda"
!      read(5,*)lambda
!      lambda = cmplx(alamr,0.0D0,r8)
      call random_number(yr)
      do j=1,nomode
      do i=0,ngrid,nomode
      alamr=sin(3.14159D0*h*float(i))
      y(i+j)=cmplx(alamr,0.0D0,r8)
      enddo
      enddo
      write(0,*)"before eigen"
      if (ikct.eq.1) then
         call sparseC16_eigen_k(this, that, y, lambda, size, rank, base, tol, nmax, ier)
      else
         call sparseC16_eigen(this, that, y, lambda, tol, nmax, ier)
      endif
      print *,'lambda = ', lambda
      print *,'tol=', tol
      print *,'nmax=', nmax

!      write(20,201)"q0=",q0,"qp=",cq1,"lambda=",real(lambda),imag(lambda)
      write(20,201)"rhoi=",imag(rho),"rho=",real(rho),"lambda=",real(lambda),imag(lambda)
!      write(20,*)
 201  format(1a5,1e12.6,2x,1a5,1e12.6,2x,1a10,1e12.6,2x,1e12.6)
      
      do j=1,nomode
      do i=j,ngrid+1,nomode
      ii=1+(i-j)/nomode
      yy=float(ii)/(float(ngrid+1)/float(nomode)+1.0D0)
      yy=yy*(xend-xini)+xini
      yya=real(y(i))
      write(19,100)yy,yya
      enddo
      enddo
      enddo
 100  format(2(3x,1e12.6))

      print *,'-------------------------------------------------'


!c      deallocate(rhs, x, res)
      deallocate(y)
      deallocate(yr)
      deallocate(aa)
      deallocate(bb)
      deallocate(C1tt) ! inertia_tt
      deallocate(C1tr) ! inertia_tr
      deallocate(C1rr) ! inertia_rr
      deallocate(C2rr) ! curvature_rr
      deallocate(C2tr) ! curvature_tr
      deallocate(C2tt) ! curvature_tt
      deallocate(C3rt) ! current_rt
      deallocate(C3tt) ! current_tt
      deallocate(C4rt) ! pressure_rt
      deallocate(C4tt) ! pressure_tt
      deallocate(C4rr) ! pressure_rr
      deallocate(C4tr) ! pressure_tr
      if (ikct.eq.1) then
        deallocate(dd)
      endif

      call sparseC16_free(this, ier)
      call sparseC16_free(that, ier)


    end subroutine sparseC16_test


      subroutine matrix
      implicit none
!      integer n,m,i1,i2,j,i,nx,im,jm,ir,jr,m1,m2,ng,nomode
      integer n,m,i1,i2,j,i,nx,im,jm,ir,jr,m1,m2,ng

      integer n2
      parameter(n=800,n2=6*n)
     
      real(r8) :: g,beta,f
      real(r8) :: xa,xb
  
      !complex(r8) :: aint,rho
      complex(r8) :: aint
      complex(r8) :: aint_c,aint_m,aint_k,aint_kp
      !integer mode(100)
      real(r8) ::ctime0,ctime1,ctime2,ctime3
      !mode(1) = 1
      !mode(2) = 2
      !mode(3) = 3
      !mode(4) = 4
      !mode(5) = 5
      !mode(6) = 6
      !mode(7) = 7
      !mode(8) = 8
      !mode(9) = 9
      !mode(10) = 10

!      write(0,*)"input,rho"
      g=0.0D0
!      read(5,*)rho

      !nomode=10
!      nomode=2
      
      ng=(nog-1)*nomode

!      write(24,*)ng-1
      
 !    
 !     call CPU_TIME(ctime0)
      call set_couple_coefficients
 !     call CPU_TIME(ctime1)
!      write(*,*)'sub set_couple_coefficients costs cpu time=',ctime1-ctime0
!      write(*,*)'begin checking:'
!      call check1  !nog=20 for checking
!      call check2(0.25D0)
!      write(*,*)'end checking'
      call check4
!
      do j=1,ng
      jr=(j-1)/nomode+1
      jm=j-((j-1)/nomode)*nomode

      call irange(ng,nomode,j,i1,i2)

      do i=i1,i2
      ir=(i-1)/nomode+1
      im=i-((i-1)/nomode)*nomode

      m1=mode(jm)
      m2=mode(im)
!
      call xrange(jr,ir,xa,xb,nx)
!  
      call matbyy(jr,ir,xa,xb,nx,m1,m2,aint,aint_c)
       !bb(j,i)=aint
       write(25,*) real(aint_c)
       write(27,*) j,i,real(aint)
    enddo
   enddo
!
   close(25)
   close(27)       

!!      call irange(n-1,j,i1,i2)
!!   write(*,*) 'before loop1'
   do j=1,ng
      jr=(j-1)/nomode+1
      jm=j-((j-1)/nomode)*nomode

      call irange(ng,nomode,j,i1,i2)
!    write(*,*) 'j=',j
      do i=i1,i2
!  
      ir=(i-1)/nomode+1
      im=i-((i-1)/nomode)*nomode
      m1=mode(jm)
      m2=mode(im)

      call xrange(jr,ir,xa,xb,nx)
      if (ikct.eq.1) then
        call mataayy(jr,ir,xa,xb,nx,m1,m2,aint,aint_c,aint_m,aint_k,aint_kp)
      else
        call matayy(jr,ir,xa,xb,nx,m1,m2,aint,aint_c,aint_m,aint_k)
      endif

      write(24,*) real(aint_c)
      write(26,*) j,i,real(aint),aimag(aint)
!  

       enddo
   enddo
!
      close(24)
      close(26)
!
   do j=1,ng
      jr=(j-1)/nomode+1
      jm=j-((j-1)/nomode)*nomode

      do i=1,ng
!
      ir=(i-1)/nomode+1
      im=i-((i-1)/nomode)*nomode
      m1=mode(jm)
      m2=mode(im)

      call xrange(jr,ir,xa,xb,nx)
      call matbyy(jr,ir,xa,xb,nx,m1,m2,aint,aint_c)
       bb(j,i)=aint
      write(270,*) j,i,real(aint)
      if (ikct.eq.1) then
        call mataayy(jr,ir,xa,xb,nx,m1,m2,aint,aint_c,aint_m,aint_k,aint_kp)
        aa(j,i)=aint
        dd(j,i)=aint_kp
      else
        call matayy(jr,ir,xa,xb,nx,m1,m2,aint,aint_c,aint_m,aint_k)
        aa(j,i)=aint  
      endif
      write(260,*) j,i,real(aint_m)
      write(261,*) j,i,real(aint_k),aimag(aint_k)
!
       enddo
   enddo
      close(260)
      close(261)
      close(270)

!  
      if(1.eq.-1)then
      do i=1,ng
!      do j=1,ng
      write(0,101)(real(bb(j,i)),j=1,ng)
      write(0,*)
      write(0,*)
!      enddo
      enddo
 101  format(16(2x,1e10.4))
      write(0,*)
      write(0,*)
      do i=1,ng
!      do j=1,ng
      write(0,101)(real(aa(j,i)),j=1,ng)
      write(0,*)
      write(0,*)
!      enddo
      enddo
      endif

!  
!      return
      end subroutine matrix
!  

      subroutine matbyy(j,i,xa,xb,nx,m1,m2,aint,aint_c)
      implicit none
      integer i,j,nx,m1,m2,nxo2,jj,j1
      real(r8) :: xa,xb,h,am1,am2,g,eps,eps2
      real(r8) :: x0,x1,x2
!      real(r8) :: fem,femp,fempp,den,gamma1,delta1,theta1,lammbda1
      !complex(r8) :: aint,f0,f1,f2,rho
      complex(r8) :: aint,f0,f1,f2
      complex(r8) :: aint_c,aint1,aint2
!  
!      eps=0.20D0
!      eps=0.375

     
!      eps2=-eps/4.0
      nxo2=nx/2
      h=(xb-xa)/float(nx)
      aint_c=(0.0D0,0.0D0)
      aint=(0.0D0,0.0D0)
      f0=(0.0D0,0.0D0)
      f1=(0.0D0,0.0D0)
      f2=(0.0D0,0.0D0)
      
!  
      am1=float(m1)
      am2=float(m2)
       j1=m1-m2
        do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!        write(30,*)'x0=',x0,'x1=',x1,'x2=',x2
!  
        f0=den(x0)*(fem(j,m1,x0)*fem(i,m2,x0)*gamma1(x0,m1,m2,j1)+ &
          fem(j,m1,x0)*femp(i,m2,x0)*delta1(x0,m1,m2,j1)- &
          theta1(x0,m1,m2,j1)*femp(j,m1,x0)*fem(i,m2,x0)+ &
          lammbda1(x0,m1,m2,j1)*femp(j,m1,x0)*femp(i,m2,x0))
        f1=den(x1)*(fem(j,m1,x1)*fem(i,m2,x1)*gamma1(x1,m1,m2,j1)+ &
          fem(j,m1,x1)*femp(i,m2,x1)*delta1(x1,m1,m2,j1)- &
          theta1(x1,m1,m2,j1)*femp(j,m1,x1)*fem(i,m2,x1)+ &
          lammbda1(x1,m1,m2,j1)*femp(j,m1,x1)*femp(i,m2,x1))
        f2=den(x2)*(fem(j,m1,x2)*fem(i,m2,x2)*gamma1(x2,m1,m2,j1)+ &
          fem(j,m1,x2)*femp(i,m2,x2)*delta1(x2,m1,m2,j1)- &
          theta1(x2,m1,m2,j1)*femp(j,m1,x2)*fem(i,m2,x2)+ &
          lammbda1(x2,m1,m2,j1)*femp(j,m1,x2)*femp(i,m2,x2))
!  
        aint=aint+h/3.0D0*(f0+4.0D0*f1+f2)
        enddo

       do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!        write(30,*)'x0=',x0,'x1=',x1,'x2=',x2
!
        f0=den(x0)*(lammbda1(x0,m1,m2,j1)*femp(j,m1,x0)*femp(i,m2,x0))
        f1=den(x1)*(lammbda1(x1,m1,m2,j1)*femp(j,m1,x1)*femp(i,m2,x1))
        f2=den(x2)*(lammbda1(x2,m1,m2,j1)*femp(j,m1,x2)*femp(i,m2,x2))
!
        aint_c=aint_c+h/3.0D0*(f0+4.0D0*f1+f2)
       enddo

!        return
        end subroutine matbyy
!
      subroutine matayy(j,i,xa,xb,nx,m1,m2,aint,aint_c,aint_m,aint_k)
        implicit none
        real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
        integer i,j,nx,m1,m2,nxo2,j1,ii,jj,mm1,mm2
        real(r8) :: xa,xb,h,am1,am2,g,eps,eps2
        real(r8) :: x0,x1,x2,dum
        !complex(r8) :: aint,f0,f1,f2,aint1,aint2,aint3,aint4,rho,rho2
        complex(r8) :: aint,f0,f1,f2,aint1,aint2,aint3,aint4,aint5,rho2
        complex(r8) :: aint_c,aint_c1,aint_c2,aint_c3
        complex(r8) :: aint_m,aint_k
        complex(r8) :: g0,g1,g2
!
        nxo2=nx/2
        h=(xb-xa)/float(nx)
        rho2=rho**2
        aint=(0.0D0,0.0D0)
        aint_c=(0.0D0,0.0D0)
        aint1=(0.0D0,0.0D0)
        aint2=(0.0D0,0.0D0)
        aint3=(0.0D0,0.0D0)
        aint4=(0.0D0,0.0D0)
        aint5=(0.0D0,0.0D0)
        aint_c1=(0.0D0,0.0D0)
        aint_c2=(0.0D0,0.0D0)
        aint_c3=(0.0D0,0.0D0)
        f0=(0.0D0,0.0D0)
        f1=(0.0D0,0.0D0)
        f2=(0.0D0,0.0D0)

!
        am1=float(m1)
        am2=float(m2)
        j1=m1-m2
        do jj=1,nxo2
          x0=xa+h*float(jj-1)*2.0D0
          x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
          x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
          f0=F(x0)*(gamma2(x0,m1,m2,j1)*akfem(j,m1,x0)*akfem(i,m2,x0)- &
             theta2(x0,m1,m2,j1)*akfem(i,m2,x0)*akfemp(j,m1,x0)+ &
             delta2(x0,m1,m2,j1)*akfemp(i,m2,x0)*akfem(j,m1,x0)+ &
             lammbda2(x0,m1,m2,j1)*akfemp(i,m2,x0)*akfemp(j,m1,x0))
          f1=F(x1)*(gamma2(x1,m1,m2,j1)*akfem(j,m1,x1)*akfem(i,m2,x1)- &
             theta2(x1,m1,m2,j1)*akfem(i,m2,x1)*akfemp(j,m1,x1)+ &
             delta2(x1,m1,m2,j1)*akfemp(i,m2,x1)*akfem(j,m1,x1)+ &
             lammbda2(x1,m1,m2,j1)*akfemp(i,m2,x1)*akfemp(j,m1,x1))
          f2=F(x2)*(gamma2(x2,m1,m2,j1)*akfem(j,m1,x2)*akfem(i,m2,x2)- &
             theta2(x2,m1,m2,j1)*akfem(i,m2,x2)*akfemp(j,m1,x2)+ &
             delta2(x2,m1,m2,j1)*akfemp(i,m2,x2)*akfem(j,m1,x2)+ &
             lammbda2(x2,m1,m2,j1)*akfemp(i,m2,x2)*akfemp(j,m1,x2))
!
          aint1=aint1+h/3.0D0*(f0+4.0D0*f1+f2)

        enddo

        do jj=1,nxo2
           x0=xa+h*float(jj-1)*2.0D0
           x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
           x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
           f0=F(x0)*(lammbda2(x0,m1,m2,j1)*akfemp(i,m2,x0)*akfemp(j,m1,x0))
           f1=F(x1)*(lammbda2(x1,m1,m2,j1)*akfemp(i,m2,x1)*akfemp(j,m1,x1))
           f2=F(x2)*(lammbda2(x2,m1,m2,j1)*akfemp(i,m2,x2)*akfemp(j,m1,x2))
!
           aint_c1=aint_c1+h/3.0D0*(f0+4.0D0*f1+f2)

        enddo
!
        do jj=1,nxo2
           x0=xa+h*float(jj-1)*2.0D0
           x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
           x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
        !f0=ak1(x0,m2)*(theta3(x0,m1,m2,j1)*fem(i,m2,x0)* &
        ! femp(j,m1,x0)+gamma3(x0,m1,m2,j1)*fem(i,m2,x0)*fem(j,m1,x0))
        !f1=ak1(x1,m2)*(theta3(x1,m1,m2,j1)*fem(i,m2,x1)* &
        !  femp(j,m1,x1)+gamma3(x1,m1,m2,j1)*fem(i,m2,x1)*fem(j,m1,x1))
        !f2=ak1(x2,m2)*(theta3(x2,m1,m2,j1)*fem(i,m2,x2)* &
        !  femp(j,m1,x2)+gamma3(x2,m1,m2,j1)*fem(i,m2,x2)*fem(j,m1,x2))
           f0=theta3(x0,m1,m2,j1)*(am1*fem(j,m1,x0)* &
              akfemp(i,m2,x0)+am2*akfem(i,m2,x0)*femp(j,m1,x0))
           f1=theta3(x1,m1,m2,j1)*(am1*fem(j,m1,x1)* &
              akfemp(i,m2,x1)+am2*akfem(i,m2,x1)*femp(j,m1,x1))
           f2=theta3(x2,m1,m2,j1)*(am1*fem(j,m1,x2)* &
              akfemp(i,m2,x2)+am2*akfem(i,m2,x2)*femp(j,m1,x2))
!
           aint2=aint2+h/3.0D0*(f0+4.0D0*f1+f2)
        enddo
!
       do jj=1,nxo2
          x0=xa+h*float(jj-1)*2.0D0
          x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
          x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
          f0=(lammbda4(x0,m1,m2,j1)*femp(j,m1,x0)*femp(i,m2,x0))/jr2(x0)
          f1=(lammbda4(x1,m1,m2,j1)*femp(j,m1,x1)*femp(i,m2,x1))/jr2(x1)
          f2=(lammbda4(x2,m1,m2,j1)*femp(j,m1,x2)*femp(i,m2,x2))/jr2(x2)

!
          aint_c2=aint_c2+h/3.0D0*(f0+4.0D0*f1+f2)
      enddo
!
       do jj=1,nxo2
          x0=xa+h*float(jj-1)*2.0D0
          x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
          x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
          f0=(theta4(x0,m1,m2,j1)*femp(j,m1,x0)*fem(i,m2,x0)+ &
             gamma4(x0,m1,m2,j1)*fem(j,m1,x0)*fem(i,m2,x0)+ &
             lammbda4(x0,m1,m2,j1)*femp(j,m1,x0)*femp(i,m2,x0)- &
             delta4(x0,m1,m2,j1)*femp(i,m2,x0)*fem(j,m1,x0))/jr2(x0)
          f1=(theta4(x1,m1,m2,j1)*femp(j,m1,x1)*fem(i,m2,x1)+ &
             gamma4(x1,m1,m2,j1)*fem(j,m1,x1)*fem(i,m2,x1)+ &
             lammbda4(x1,m1,m2,j1)*femp(j,m1,x1)*femp(i,m2,x1)- &
             delta4(x1,m1,m2,j1)*femp(i,m2,x1)*fem(j,m1,x1))/jr2(x1)
          f2=(theta4(x2,m1,m2,j1)*femp(j,m1,x2)*fem(i,m2,x2)+ &
             gamma4(x2,m1,m2,j1)*fem(j,m1,x2)*fem(i,m2,x2)+ &
             lammbda4(x2,m1,m2,j1)*femp(j,m1,x2)*femp(i,m2,x2)- &
             delta4(x2,m1,m2,j1)*femp(i,m2,x2)*fem(j,m1,x2))/jr2(x2)


          aint3=aint3+h/3.0D0*(f0+4.0D0*f1+f2)
       enddo
!
       if (j1==0) then
         do jj=1,nxo2
            x0=xa+h*float(jj-1)*2.0D0
            x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
            x2=xa+h*(float(jj-1)*2.0D0+2.0D0)       
!
            g0=gfun(x0,m1)
            g1=gfun(x1,m1)
            g2=gfun(x2,m1)
!
            f0=g0*rho2*(fempp(i,m1,x0)+femp(i,m1,x0)*(1.0/x0+denp(x0)/den(x0)) &
               -am1**2/x0**2*fem(i,m1,x0))* &
                  x0*(fempp(j,m1,x0)+femp(j,m1,x0)/x0-am1**2/x0**2*fem(j,m1,x0))

            f1=g1*rho2*(fempp(i,m1,x1)+femp(i,m1,x1)*(1.0/x1+denp(x1)/den(x1)) &
               -am1**2/x1**2*fem(i,m1,x1))* &
                  x1*(fempp(j,m1,x1)+femp(j,m1,x1)/x1-am1**2/x1**2*fem(j,m1,x1))

            f2=g2*rho2*(fempp(i,m1,x2)+femp(i,m1,x2)*(1.0/x2+denp(x2)/den(x2)) &
               -am1**2/x2**2*fem(i,m1,x2))* &
                  x2*(fempp(j,m1,x2)+femp(j,m1,x2)/x2-am1**2/x2**2*fem(j,m1,x2))
!
            aint5=aint5+h/3.0D0*(f0+4.0D0*f1+f2)
        
!
         enddo
           
       endif
!
!       fix radial grid index from 0 to nog
!
       mm1 = m1 - midx
       mm2 = m2 - midx
       aint4 = 1.0/4.0*b0_h*aint_KCT4D(i,j,mm2,mm1) ! factor (2*pi)^2 reduced
       !aint4 = 1.0/4.0*b0_h*aimag(aint_KCT4D(i,j,mm2,mm1))
       !  aint=aint1-aint2+2*aint3/e**2
       aint_k=aint4/e**2
        ! aint_m=aint1-aint2+2*aint3/e**2
       aint_m=aint1+aint2+2*aint3/e**2+aint5
       !aint_m=aint1+2.0*aint3/e**2+aint5
       aint=aint_m+aint_k
       !   aint = aint1 - aint2 + 2.0*aint3/e**2 + aint4/e**2

       aint_c=aint_c1+2.0*aint_c2/e**2

       end subroutine matayy
!
      subroutine mataayy(j,i,xa,xb,nx,m1,m2,aint,aint_c,aint_m,aint_k,aint_kp)
      implicit none
      real,parameter ::pi=3.14159265358979, two_pi=6.28318530717959d0
      integer i,j,nx,m1,m2,nxo2,j1,ii,jj,mm1,mm2
      real(r8) :: xa,xb,h,am1,am2,g,eps,eps2
      real(r8) :: x0,x1,x2,dum
!      real(r8) :: ak1,akfemp,akfem,gamma2,delta2,theta2,lammbda2,F, &
!	  theta3,gamma3,delta4,gamma4,theta4,lammbda4,jr2,fem,femp
      !complex(r8) :: aint,f0,f1,f2,aint1,aint2,aint3,aint4,rho 
      complex(r8) :: aint,f0,f1,f2,aint1,aint2,aint3,aint4
      complex(r8) :: aint_c,aint_c1,aint_c2,aint_c3
      complex(r8) :: aint_m,aint_k,aint_kp
!        e=0.1
!  
!      eps=0.20D0
!      eps=0.375
       
     
!      eps2=-eps/4.0
      nxo2=nx/2
      h=(xb-xa)/float(nx)
      aint=(0.0D0,0.0D0)
      aint_c=(0.0D0,0.0D0)
      aint1=(0.0D0,0.0D0)
      aint2=(0.0D0,0.0D0)
      aint3=(0.0D0,0.0D0)
      aint4=(0.0D0,0.0D0)
      aint_c1=(0.0D0,0.0D0)
      aint_c2=(0.0D0,0.0D0)
      aint_c3=(0.0D0,0.0D0)

      f0=(0.0D0,0.0D0)
      f1=(0.0D0,0.0D0)
      f2=(0.0D0,0.0D0)
      
!  
      am1=float(m1)
      am2=float(m2)
       j1=m1-m2
        do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
        f0=F(x0)*(gamma2(x0,m1,m2,j1)*akfem(j,m1,x0)*akfem(i,m2,x0)- &
          theta2(x0,m1,m2,j1)*akfem(i,m2,x0)*akfemp(j,m1,x0)+ &
          delta2(x0,m1,m2,j1)*akfemp(i,m2,x0)*akfem(j,m1,x0)+ &
          lammbda2(x0,m1,m2,j1)*akfemp(i,m2,x0)*akfemp(j,m1,x0))
        f1=F(x1)*(gamma2(x1,m1,m2,j1)*akfem(j,m1,x1)*akfem(i,m2,x1)- &
          theta2(x1,m1,m2,j1)*akfem(i,m2,x1)*akfemp(j,m1,x1)+ &
          delta2(x1,m1,m2,j1)*akfemp(i,m2,x1)*akfem(j,m1,x1)+ &
          lammbda2(x1,m1,m2,j1)*akfemp(i,m2,x1)*akfemp(j,m1,x1))
        f2=F(x2)*(gamma2(x2,m1,m2,j1)*akfem(j,m1,x2)*akfem(i,m2,x2)- &
          theta2(x2,m1,m2,j1)*akfem(i,m2,x2)*akfemp(j,m1,x2)+ &
          delta2(x2,m1,m2,j1)*akfemp(i,m2,x2)*akfem(j,m1,x2)+ &
          lammbda2(x2,m1,m2,j1)*akfemp(i,m2,x2)*akfemp(j,m1,x2))
!       
        aint1=aint1+h/3.0D0*(f0+4.0D0*f1+f2)

        enddo

      do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
        f0=F(x0)*(lammbda2(x0,m1,m2,j1)*akfemp(i,m2,x0)*akfemp(j,m1,x0))
        f1=F(x1)*(lammbda2(x1,m1,m2,j1)*akfemp(i,m2,x1)*akfemp(j,m1,x1))
        f2=F(x2)*(lammbda2(x2,m1,m2,j1)*akfemp(i,m2,x2)*akfemp(j,m1,x2))
!
        aint_c1=aint_c1+h/3.0D0*(f0+4.0D0*f1+f2)

        enddo

!        return
!
        do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
        !f0=ak1(x0,m2)*(theta3(x0,m1,m2,j1)*fem(i,m2,x0)* &
        ! femp(j,m1,x0)+gamma3(x0,m1,m2,j1)*fem(i,m2,x0)*fem(j,m1,x0))
        !f1=ak1(x1,m2)*(theta3(x1,m1,m2,j1)*fem(i,m2,x1)* &
        !  femp(j,m1,x1)+gamma3(x1,m1,m2,j1)*fem(i,m2,x1)*fem(j,m1,x1))
        !f2=ak1(x2,m2)*(theta3(x2,m1,m2,j1)*fem(i,m2,x2)* &
        !  femp(j,m1,x2)+gamma3(x2,m1,m2,j1)*fem(i,m2,x2)*fem(j,m1,x2))
        f0=theta3(x0,m1,m2,j1)*(am1*fem(j,m1,x0)* &
              akfemp(i,m2,x0)+am2*akfem(i,m2,x0)*femp(j,m1,x0))
        f1=theta3(x1,m1,m2,j1)*(am1*fem(j,m1,x1)* &
              akfemp(i,m2,x1)+am2*akfem(i,m2,x1)*femp(j,m1,x1))
        f2=theta3(x2,m1,m2,j1)*(am1*fem(j,m1,x2)* &
              akfemp(i,m2,x2)+am2*akfem(i,m2,x2)*femp(j,m1,x2))
!
        aint2=aint2+h/3.0D0*(f0+4.0D0*f1+f2)
		enddo
!
       do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
        f0=(lammbda4(x0,m1,m2,j1)*femp(j,m1,x0)*femp(i,m2,x0))/jr2(x0)
        f1=(lammbda4(x1,m1,m2,j1)*femp(j,m1,x1)*femp(i,m2,x1))/jr2(x1)
        f2=(lammbda4(x2,m1,m2,j1)*femp(j,m1,x2)*femp(i,m2,x2))/jr2(x2)

!
        aint_c2=aint_c2+h/3.0D0*(f0+4.0D0*f1+f2)
         enddo

!
        do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!
        f0=(theta4(x0,m1,m2,j1)*femp(j,m1,x0)*fem(i,m2,x0)+ &
            gamma4(x0,m1,m2,j1)*fem(j,m1,x0)*fem(i,m2,x0)+ &
            lammbda4(x0,m1,m2,j1)*femp(j,m1,x0)*femp(i,m2,x0)- &
            delta4(x0,m1,m2,j1)*femp(i,m2,x0)*fem(j,m1,x0))/jr2(x0)
        f1=(theta4(x1,m1,m2,j1)*femp(j,m1,x1)*fem(i,m2,x1)+ &
            gamma4(x1,m1,m2,j1)*fem(j,m1,x1)*fem(i,m2,x1)+ &
            lammbda4(x1,m1,m2,j1)*femp(j,m1,x1)*femp(i,m2,x1)- &
            delta4(x1,m1,m2,j1)*femp(i,m2,x1)*fem(j,m1,x1))/jr2(x1)
        f2=(theta4(x2,m1,m2,j1)*femp(j,m1,x2)*fem(i,m2,x2)+ &
            gamma4(x2,m1,m2,j1)*fem(j,m1,x2)*fem(i,m2,x2)+ &
            lammbda4(x2,m1,m2,j1)*femp(j,m1,x2)*femp(i,m2,x2)- &
            delta4(x2,m1,m2,j1)*femp(i,m2,x2)*fem(j,m1,x2))/jr2(x2)
       
        
        aint3=aint3+h/3.0D0*(f0+4.0D0*f1+f2)
		enddo
!
!       fix radial grid index from 0 to nog
!       
        mm1 = m1 - midx
        mm2 = m2 - midx
        aint4 = 1.0/4.0*b0_h*aint_KCT4D(i,j,mm2,mm1) ! factor (2*pi)^2 reduced
        !dum = aimag(aint4)                        ! betah(0)
        !write(271,*) dum
       !  aint=aint1-aint2+2*aint3/e**2
         aint_k=aint4/e**2
         aint_kp = 1.0/4.0*b0_h*aint_KCT4Dp(i,j,mm2,mm1)/e**2
        ! aint_m=aint1-aint2+2*aint3/e**2   
         aint_m=aint1+aint2+2*aint3/e**2
         aint=aint_m+aint_k       
       !  aint=aint1-aint2+2*aint3/e**2+cmplx(0.0,dum/e**2) ! 1e-2 : betah(0)
       !   aint = aint1 - aint2 + 2.0*aint3/e**2 + aint4/e**2
 
        aint_c=aint_c1+2.0*aint_c2/e**2
         
        end subroutine mataayy
!
!
     subroutine matbxx(j,i,xa,xb,nx,m,aint)
      implicit none
      integer i,j,nx,m,nxo2,jj
      real(r8) :: xa,xb,h,am,g
      real(r8) :: x0,x1,x2
!      real(r8) :: fem,femp,fempp,den
      complex(r8) :: aint,f0,f1,f2
      nxo2=nx/2
      h=(xb-xa)/float(nx)
      aint=(0.0D0,0.0D0)
      f0=(0.0D0,0.0D0)
      f1=(0.0D0,0.0D0)
      f2=(0.0D0,0.0D0)
      
!  
      am=float(m)
        do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!  
        f0=x0*den(x0)*femp(i,m,x0)*femp(j,m,x0)+ &
           den(x0)*am**2/x0* fem(i,m,x0)*fem(j,m,x0)
        f1=x1*den(x1)*femp(i,m,x1)*femp(j,m,x1)+ &
           den(x1)*am**2/x1* fem(i,m,x1)*fem(j,m,x1)
        f2=x2*den(x2)*femp(i,m,x2)*femp(j,m,x2)+ &
           den(x2)*am**2/x2* fem(i,m,x2)*fem(j,m,x2)
!  
        aint=aint+h/3.0D0*(f0+4.0D0*f1+f2)
        enddo
!        return
        end subroutine matbxx


      subroutine matbepsxx(j,i,xa,xb,nx,m1,m2,aint)
      implicit none
      integer i,j,nx,m1,m2,nxo2,jj
      real(r8) :: xa,xb,h,am1,am2,g,eps,eps2
      real(r8) :: x0,x1,x2
!      real(r8) :: fem,femp,fempp,den
      complex(r8) :: aint,f0,f1,f2
!  
!      eps=0.20D0
      eps=0.375

     
      eps2=-eps/4.0
      nxo2=nx/2
      h=(xb-xa)/float(nx)
      aint=(0.0D0,0.0D0)
      f0=(0.0D0,0.0D0)
      f1=(0.0D0,0.0D0)
      f2=(0.0D0,0.0D0)
      
!  
      am1=float(m1)
      am2=float(m2)
        do jj=1,nxo2
        x0=xa+h*float(jj-1)*2.0D0
        x1=xa+h*(float(jj-1)*2.0D0+1.0D0)
        x2=xa+h*(float(jj-1)*2.0D0+2.0D0)
!  
        f0=x0**2*den(x0)*femp(i,m1,x0)*femp(j,m2,x0)+ &
          eps2*den(x0)*am1*am2*fem(i,m1,x0)*fem(j,m2,x0)
        f1=x1**2*den(x1)*femp(i,m1,x1)*femp(j,m2,x1)+ &
          eps2*den(x1)*am1*am2*fem(i,m1,x1)*fem(j,m2,x1)
        f2=x2**2*den(x2)*femp(i,m1,x2)*femp(j,m2,x2)+ &
          eps2*den(x2)*am1*am2* fem(i,m1,x2)*fem(j,m2,x2)
!  
        aint=aint+h/3.0D0*(f0+4.0D0*f1+f2)*eps
        enddo
!        return
        end subroutine matbepsxx
!

    subroutine sparseC16_error(this, nunit, ier)
      
      ! Error handler
      
      implicit none
      type(sparseC16_obj) :: this
      integer, intent(in) :: nunit ! output channel
      integer, intent(in) :: ier
      
      character(*), parameter :: header = 'SparseC16::'

    select case(ier)
       case(1)
          write(nunit, *) header//'ERROR while initializing'
       case(2)
          write(nunit, *) header//'ERROR while deallocating'
       case(3)
          write(nunit, *) &
               & header//'ERROR counter exceeded declared size ', this%size,' rank=', this%size_ptr-1
       case(4)
          write(nunit, *) header//'ERROR while saving'
       case(5)
          write(nunit, *) header//'ERROR while loading'
       case(6)
          write(nunit, *) header//'ERROR occurred in test unit'
       case(7)
          write(nunit, *) header//'ERROR method not implemented'
       case(8)
          write(nunit, *) header//'ERROR occurred in solve'
       case(9)
          write(nunit, *) header//'ERROR occurred when allocating "irow"'
       case(10)
          write(nunit, *) header//'ERROR occurred when allocating "values"'
       case(11)
          write(nunit, *) header//'ERROR occurred when allocating "jcol_ptr"'
       case(12)
          write(nunit, *) header//'ERROR occurred when copying'
       case(13)
          write(nunit, *) header//'WARNING "base" mismatch'
       case(14)
          write(nunit, *) header//'WARNING "nnz" mismatch'
       case(15)
          write(nunit, *) header//'WARNING "size_ptr" mismatch'
       case(16)
          write(nunit, *) header//'WARNING "size" mismatch'
       case(17)
          write(nunit, *) header//'WARNING "col_counter" mismatch'
       case(18)
          write(nunit, *) header//'WARNING "jcol_ptr" array differ'
       case(19)
          write(nunit, *) header//'WARNING "irow" array differ'
       case(20)
          write(nunit, *) header//'ERROR occurred  in det'
       
       case default
    end select    
      
      
    end subroutine sparseC16_error
    
!===================================================================================

! real*8

    subroutine sparseR8_init(this, rank, size, base, ier)

      ! Constructor.
      !
      ! rank: matrix rank (or approximation if not known)
      ! size: >= no of non-zero values to store. If 
      ! size turns out to be too tight, fresh arrays will be 
      ! reallocated and the data will be copied, however this may
      ! can cause performance degradation.
      !
      ! base: indexing base (0 like C or 1 like Fortran)

      implicit none
      type(sparseR8_obj) :: this
      integer, intent(in) :: rank
      integer, intent(in) :: size
      integer, intent(in) :: base
      integer, intent(out) :: ier

      integer jer

      ier = 0

      this % permc_spec = 2 ! seems to work best

      this % nnz = 0
      this % size_ptr = rank+1
      this % size = size
      this % base = base
      this % col_counter = 0
      allocate(this%irow(this%size), stat=jer)
      if(jer /= 0 ) then
         ier = 9
         return
      end if
      allocate(this%values(this%size), stat=jer)
      if(jer /= 0 ) then
         ier = 10
         return
      end if
      allocate(this%jcol_ptr(this%size_ptr), stat=jer)
      if(jer /= 0 ) then
         ier = 11
         return
      end if
      this%irow = 0
      this%jcol_ptr = 0
      this%values = 0
      if(jer /=0 ) ier = 1
      
    end subroutine sparseR8_init
    
    subroutine sparseR8_free(this, ier)

      ! Destructor

      implicit none
      type(sparseR8_obj) :: this
      integer, intent(out) :: ier

      integer jer

      ier = 0
      deallocate(this%irow, stat=jer)
      deallocate(this%values, stat=jer)
      deallocate(this%jcol_ptr, stat=jer)
      if(jer/=0) ier = 2
      
    end subroutine sparseR8_free

    subroutine sparseR8_copy(this, that, ier)

      ! Copy this into that. that should be unitialized.

      implicit none
      type(sparseR8_obj) :: this, that
      integer, intent(out) :: ier

      integer jer

      ier = 0
      that % nnz = this % nnz
      that % size_ptr = this % size_ptr
      that % size = this % size
      that % base = this % base
      that % col_counter = this % col_counter
      allocate(that%irow(that%size), stat=jer)
      if(jer /= 0 ) then
         ier = 9
         return
      end if
      allocate(that%values(that%size), stat=jer)
      if(jer /= 0 ) then
         ier = 10
         return
      end if
      allocate(that%jcol_ptr(that%size_ptr), stat=jer)
      if(jer /= 0 ) then
         ier = 11
         return
      end if
      that % irow = this%irow
      that % jcol_ptr = this%jcol_ptr
      that % values = this%values
      if(jer /=0 ) ier = 12
      
      
    end subroutine sparseR8_copy

    subroutine sparseR8_compare_shapes(this, that, ier)

      ! Compare two sparse matrices

      implicit none
      type(sparseR8_obj) :: this, that
      integer, intent(out) :: ier

      integer i

      ier = 0
      if(this%base /= that%base) then
         ier = 13
         return
      end if
      if(this%nnz /= that%nnz) then
         ier = 14
         return
      end if
      if(this%size_ptr /= that%size_ptr) then
         ier = 15
         return
      end if
      if(this%size /= that%size) then
         ier = 16
         return
      end if
      if(this%col_counter /= that%col_counter) then
         ier = 17
         return
      end if
      do i = 1, this%size_ptr
         if(this%jcol_ptr(i) /= that%jcol_ptr(i)) then
            ier = 18
            return 
         end if
      enddo
      do i = 1, this%size
         if(this%irow(i) /= that%irow(i)) then
            ier = 19
            return      
         end if
      enddo


    end subroutine sparseR8_compare_shapes

    subroutine sparseR8_set_next(this, row, newcol_flag, value, ier)

      ! Set next value. Must iterate down each column 
      ! running from left to right.
      !
    
      implicit none
      type(sparseR8_obj) :: this
      integer, intent(in) :: row 
      logical, intent(in) :: newcol_flag ! .TRUE. if new column, .FALSE. otherwise
      real(r8), intent(in) :: value 
      integer, intent(out) :: ier

      ier = 0
      this % nnz = this % nnz + 1
      if(this % nnz > this % size) then         
         ier = 3
         return
      endif
      this % irow(this%nnz) = row
      this % values(this%nnz) = value
      if(newcol_flag) then
         this % col_counter = this % col_counter + 1
         this % jcol_ptr(this % col_counter) = this%nnz - 1 + this%base
      endif
      this % jcol_ptr(this % col_counter+1) = this%nnz + this%base
      
    end subroutine sparseR8_set_next

    subroutine sparseR8_get_size(this, nrows, ncols, ier)

      ! Return no of rows and columns

      implicit none
      type(sparseR8_obj) :: this
      integer, intent(out) :: nrows, ncols
      integer, intent(out) :: ier

      ier = 0
      nrows = maxval(this % irow)
      ncols = this % col_counter
    end subroutine sparseR8_get_size

    subroutine sparseR8_solve(this, rhs, ier)

      ! Solve matrix system
      
      implicit none
      type(sparseR8_obj) :: this
      real(r8), intent(inout) :: rhs(:) ! right-hand side vector, solution vector upon return
      integer, intent(out) :: ier

      integer n, jer
      integer :: handle(8) ! opaque handle

      ier = 0
      n = size(rhs)

      ! note: we pass here address (1) of vectors to force passing the 
      ! arguments by reference. 
      call dsupralu_new(handle, &
           & this%values(1), this%irow(1), this%jcol_ptr(1), this%nnz, n, &
           & jer)
      call dsupralu_colperm(handle, this % permc_spec, jer)
      call dsupralu_lu(handle, jer)
      call dsupralu_solve(handle, rhs(1), jer)
      call dsupralu_del(handle, jer)

      if(jer/=0) ier = 8

    end subroutine sparseR8_solve

    subroutine sparseR8_det(this, res_mantissa, res_exponent, ier)

      ! Return the determinant of A in the form of 
      !
      ! res_mantissa * 2**res_exponent

      implicit none
      type(sparseR8_obj) :: this
      real(r8), intent(out) :: res_mantissa
      integer, intent(out) :: res_exponent
      integer, intent(out) :: ier
      
      integer :: handle(8) ! opaque handle
      integer :: n, jer

      ier = 0
      n = this % col_counter
      call dsupralu_new(handle, &
           & this%values(1), this%irow(1), this%jcol_ptr(1), this%nnz, n, &
           & jer)
      call dsupralu_colperm(handle, this % permc_spec, jer)
      call dsupralu_lu(handle, jer)
      call dsupralu_determinant(handle, res_mantissa, res_exponent, jer)
      call dsupralu_del(handle, jer)
      if(jer/=0) ier = 20

    end subroutine sparseR8_det

    subroutine sparseR8_A_dot_x(this, x, res, ier)

      ! matrix . x => res

      implicit none
      type(sparseR8_obj) :: this
      real(r8), intent(in) :: x(:) ! vector
      real(r8), intent(out) :: res(:) ! result
      integer, intent(out) :: ier

      integer i, k, rank, j

      ier = 0
      res = 0
      rank = maxval(this % irow) + 1 - this%base

      
      do j = this%base+1, rank
         do k = this%jcol_ptr(j)+1-this%base, this%jcol_ptr(j+1)-this%base
            i = this%irow(k) + 1 - this%base
            res(i) = res(i) + this%values(k)*x(j)
         enddo
      enddo
      
    end subroutine sparseR8_A_dot_x

    subroutine sparseR8_save(this, filename, ier)

    ! Save state in netCDF file 'filename'

    use ezcdf
    implicit none
    type(sparseR8_obj) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer ncid, jer
    ier = 0

    call cdf_open(ncid, filename, 'w', jer)

    call cdf_define(ncid, 'version', (/len_trim(version),0,0/), 'CHAR', jer)
    call cdf_define(ncid, 'size', this%size, jer)
    call cdf_define(ncid, 'size_ptr', this%size_ptr, jer)
    call cdf_define(ncid, 'nnz', this%nnz, jer)
    call cdf_define(ncid, 'col_counter', this%col_counter, jer)
    call cdf_define(ncid, 'base', this%base, jer)
    call cdf_define(ncid, 'irow', this%irow(1:this%nnz), jer)
    call cdf_define(ncid, 'jcol_ptr', this%jcol_ptr, jer)
    call cdf_define(ncid, 'values', this%values(1:this%nnz), jer)

    call cdf_write(ncid, 'version', version, jer)
    call cdf_write(ncid, 'size', this%size, jer)
    call cdf_write(ncid, 'size_ptr', this%size_ptr, jer)
    call cdf_write(ncid, 'nnz', this%nnz, jer)
    call cdf_write(ncid, 'col_counter', this%col_counter, jer)
    call cdf_write(ncid, 'base', this%base, jer)
    call cdf_write(ncid, 'irow', this%irow(1:this%nnz), jer)
    call cdf_write(ncid, 'jcol_ptr', this%jcol_ptr, jer)
    call cdf_write(ncid, 'values', this%values(1:this%nnz), jer)

    call cdf_close(ncid, jer) 
    if(jer /= 0) ier = 4

    end subroutine sparseR8_save

    subroutine sparseR8_load(this, filename, ier)

      ! Load state from netCDF file 'filename'

    use ezcdf
    implicit none
    type(sparseR8_obj) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer ncid, jer
    ier = 0

    call cdf_open(ncid, filename, 'r', jer)

    call cdf_read(ncid, 'size', this%size, jer)
    call cdf_read(ncid, 'size_ptr', this%size_ptr, jer)
    call cdf_read(ncid, 'nnz', this%nnz, jer)
    call cdf_read(ncid, 'col_counter', this%col_counter, jer)
    call cdf_read(ncid, 'base', this%base, jer)
    allocate(this % irow(this%size), stat=jer)
    allocate(this % jcol_ptr(this%size_ptr), stat=jer)
    allocate(this % values(this%size), stat=jer)
    if(jer /= 0) then 
       ier = 1
       return
    endif
    call cdf_read(ncid, 'irow', this%irow, jer)
    call cdf_read(ncid, 'jcol_ptr', this%jcol_ptr, jer)
    call cdf_read(ncid, 'values', this%values, jer)

    call cdf_close(ncid, jer) 
    if(jer /= 0) ier = 5
      
    end subroutine sparseR8_load

    subroutine sparseR8_test(this, ier)

      ! Test unit

      implicit none
      type(sparseR8_obj) :: this, that, amat, bmat
      integer, intent(out) :: ier

      integer rank, size, base, i, nmax
      ! **not needed: eigen vector C16 only (dmc): real(r8) :: tol
      real(r8), allocatable :: rhs(:), x(:), res(:), y(:)
      real(r8) :: lambda, adet, bdet, adet_mantissa, bdet_mantissa
      integer :: adet_exponent, bdet_exponent
      real(r8) :: cum_error

      ier = 0
      cum_error = 0
      print *, 'testing sparseR8...', version

      base = 0 ! C-like
      rank = 5
      size = 10
      call sparseR8_init(this, rank, size, base, ier)
      ! matrix to store
      ! [ 0 2 x 7 x ]
      ! [ 1 3 x x x ]
      ! [ x x 4 x x ]
      ! [ x x 5 8 x ]
      ! [ x x 6 x 9 ]
      ! the result should be 
      ! irow = [0,1,0,1,2,3,4,0,3,4]
      ! jcol_ptr = [0,2,4,7,9]
      ! values = [0,1,2,3,4,5,6,7,8,9]
      print *, '...assembling'
      call sparseR8_set_next(this, 0, .TRUE., 0._r8, ier)
      call sparseR8_set_next(this, 1, .FALSE., 1._r8, ier)
      call sparseR8_set_next(this, 0, .TRUE., 2._r8, ier)
      call sparseR8_set_next(this, 1, .FALSE., 3._r8, ier)
      call sparseR8_set_next(this, 2, .TRUE., 4._r8, ier)
      call sparseR8_set_next(this, 3, .FALSE., 5._r8, ier)
      call sparseR8_set_next(this, 4, .FALSE., 6._r8, ier)
      call sparseR8_set_next(this, 0, .TRUE., 7._r8, ier)
      call sparseR8_set_next(this, 3, .FALSE., 8._r8, ier)
      call sparseR8_set_next(this, 4, .TRUE., 9._r8, ier)

      print *, '...saving'
      call sparseR8_save(this, 'sparseR8.nc', ier)
      call sparseR8_error(this, 6, ier)
      cum_error = cum_error + ier

      print *,'...freeing'
      call sparseR8_free(this, ier)
      cum_error = cum_error + ier

      print *,'...loading'
      call sparseR8_load(this, 'sparseR8.nc', ier)
      call sparseR8_error(this, 6, ier)
      cum_error = cum_error + ier

      print *, '...saving'
      call sparseR8_save(this, 'sparseR8_2.nc', ier)
      call sparseR8_error(this, 6, ier)
      cum_error = cum_error + ier

      allocate(rhs(rank), x(rank), res(rank))
      rhs = (/ (i, i=1,rank) /)
      x = rhs
      print *, '...solving'
      call sparseR8_solve(this, x, ier)      
      cum_error = cum_error + ier
      call sparseR8_A_dot_x(this, x, res, ier)
      cum_error = cum_error + ier
      print *,'error = ', sqrt(abs(sum((res - rhs)**2)/rank))
      cum_error = cum_error + sqrt(abs(sum((res - rhs)**2)/rank))

      call sparseR8_copy(this, that, ier)
      call sparseR8_error(this, 6, ier)      
      cum_error = cum_error + ier

      that % values = 2 * this % values

      call sparseR8_free(this, ier)
      cum_error = cum_error + ier
      call sparseR8_free(that, ier)
      cum_error = cum_error + ier

!!$>> amat=[1. -1 0; -1 2 3; 0 1 1.];
!!$>> bmat=[1 0.5 0; 0.5 1 0.5; 0 0.5 1];
!!$>> det(amat)
!!$
!!$ans =
!!$
!!$    -2
!!$
!!$>> det(bmat)
!!$
!!$ans =
!!$
!!$    0.5000


      call sparseR8_init(amat, 3, 7, 0, ier)
      call sparseR8_set_next(amat, 0, .TRUE.,  1._r8, ier)
      call sparseR8_set_next(amat, 1, .FALSE.,-1._r8, ier)
      call sparseR8_set_next(amat, 0, .TRUE., -1._r8, ier)
      call sparseR8_set_next(amat, 1, .FALSE., 2._r8, ier)
      call sparseR8_set_next(amat, 2, .FALSE., 1._r8, ier)
      call sparseR8_set_next(amat, 1, .TRUE.,  3._r8, ier)
      call sparseR8_set_next(amat, 2, .FALSE., 1._r8, ier)
      cum_error = cum_error + ier

      call sparseR8_init(bmat, 3, 7, 0, ier)
      call sparseR8_set_next(bmat, 0, .TRUE.,  1._r8, ier)
      call sparseR8_set_next(bmat, 1, .FALSE., .5_r8, ier)
      call sparseR8_set_next(bmat, 0, .TRUE.,  .5_r8, ier)
      call sparseR8_set_next(bmat, 1, .FALSE., 1._r8, ier)
      call sparseR8_set_next(bmat, 2, .FALSE., .5_r8, ier)
      call sparseR8_set_next(bmat, 1, .TRUE.,  .5_r8, ier)
      call sparseR8_set_next(bmat, 2, .FALSE., 1._r8, ier)
      cum_error = cum_error + ier

      amat % permc_spec = 2
      call sparseR8_det(amat, adet_mantissa, adet_exponent, ier)
      call sparseR8_det(bmat, bdet_mantissa, bdet_exponent, ier)
      adet = adet_mantissa * 2._r8**adet_exponent
      bdet = bdet_mantissa * 2._r8**bdet_exponent
      print *,'det A = ',adet,' det B = ', bdet
      cum_error = cum_error + abs(adet - real(-2,r8))
      cum_error = cum_error + abs(bdet - 0.5_r8)

      deallocate(rhs, x, res)
      call sparseR8_free(amat, ier)
      call sparseR8_free(bmat, ier)

      print *,'***********************************************************'
      print *,'CUMULATIVE ERROR = ', cum_error
      if(cum_error > 1.e-10_r8) then
         print *,' **TEST FAILED**'
      else
         print *,' TEST SUCCEEDED '
      endif
      print *,'***********************************************************'

    end subroutine sparseR8_test

    subroutine sparseR8_error(this, nunit, ier)
      
      ! Error handler
      
      implicit none
      type(sparseR8_obj) :: this
      integer, intent(in) :: nunit ! output channel
      integer, intent(in) :: ier
      
      character(*), parameter :: header = 'SparseR8::'

    select case(ier)
       case(1)
          write(nunit, *) header//'ERROR while initializing'
       case(2)
          write(nunit, *) header//'ERROR while deallocating'
       case(3)
          write(nunit, *) &
               & header//'ERROR counter exceeded declared size ', this%size,' rank=', this%size_ptr-1
       case(4)
          write(nunit, *) header//'ERROR while saving'
       case(5)
          write(nunit, *) header//'ERROR while loading'
       case(6)
          write(nunit, *) header//'ERROR occurred in test unit'
       case(7)
          write(nunit, *) header//'ERROR method not implemented'
       case(8)
          write(nunit, *) header//'ERROR occurred in solve'
       case(9)
          write(nunit, *) header//'ERROR occurred when allocating "irow"'
       case(10)
          write(nunit, *) header//'ERROR occurred when allocating "values"'
       case(11)
          write(nunit, *) header//'ERROR occurred when allocating "jcol_ptr"'
       case(12)
          write(nunit, *) header//'ERROR occurred when copying'
       case(13)
          write(nunit, *) header//'WARNING "base" mismatch'
       case(14)
          write(nunit, *) header//'WARNING "nnz" mismatch'
       case(15)
          write(nunit, *) header//'WARNING "size_ptr" mismatch'
       case(16)
          write(nunit, *) header//'WARNING "size" mismatch'
       case(17)
          write(nunit, *) header//'WARNING "col_counter" mismatch'
       case(18)
          write(nunit, *) header//'WARNING "jcol_ptr" array differ'
       case(19)
          write(nunit, *) header//'WARNING "irow" array differ'
       case(20)
          write(nunit, *) header//'ERROR occurred  in det'
       
       case default
    end select    
      
      
    end subroutine sparseR8_error

end module supralu_mod
