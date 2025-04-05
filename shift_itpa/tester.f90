program tester
  use mpi
  use test2_mod
  use vars
  implicit none
  integer input, i
  integer :: ierr, numprocs, proc_num
  real :: ctime0,ctime1,tol
  ! read data and compute aint_KCT4D parallelly
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,proc_num,ierr)
  if (proc_num == 0) then
    write(*,*) 'kaec ',numprocs,' processes computing...'
  endif

  call CPU_TIME(ctime0)
 tol = 1.e-12
 do  i = 1,3
  ! calculate matrix A(omg)
  call main

  !call CPU_TIME(ctime0)
 
  ! solve matrix serially
  if (proc_num==0) then
  print *, 'GACO module tester'
  print *, '0: all test runs'
  print *, '1: plotmtv'
  print *, '2: indexer'
  print *, '2: gabor'
  print *, '3: kzgrid'
  print *, '4: xgrid'
  print *, '5: sparseC16'
    input=6
  select case(input)

     case(1)
        call plotmtv_test_run
     case(2)
        call indexer_test_run
     case(3)
        call gabor_test_run
     case(4)
        call kzgrid_test_run
     case(5)
        call xgrid_test_run
     case(6)
        call sparseC16_test_run
     case default
        call plotmtv_test_run
        call indexer_test_run
        call gabor_test_run
        call kzgrid_test_run
        call xgrid_test_run
        call sparseC16_test_run
  end select
  !if(res < tol) exit
  endif
enddo
  call CPU_TIME(ctime1)
 if (proc_num == 0) then
  write(*,*)'eigen solve mod cpu time=',&
             ctime1-ctime0
  print *,'residue', res ,'nmax=',i
 endif
  
  call mpi_finalize(ierr)
end program tester

subroutine plotmtv_test_run
  use plotmtv
  implicit none
  type(plotmtv_obj) :: this_plotmtv

  call plotmtv_test(this_plotmtv)
  
end subroutine plotmtv_test_run

subroutine indexer_test_run
  use indexer
  implicit none
  integer ier
  type(indexer_obj) :: this_indexer

  call indexer_test(this_indexer, ier)
  if(ier /=0 ) then
     print *,'indexer_test **FAILED**'
  else
     print *,'indexer_test SUCCEDED'
  endif

end subroutine indexer_test_run

subroutine gabor_test_run
  use gabor
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)
  real(r8), parameter :: twopi = 6.2831853071795865_r8
  real(r8) :: x, z, delta_x, delta_k, sigma, ak
  real(r8), dimension(3) :: xs, zs, delta_xs, delta_ks, sigmas, aks

  x=0.212346_r8; delta_x = 0.1_r8; z=0.1_r8;
  sigma=1._r8/delta_x**2
  delta_k = twopi/(delta_x * 1.3_r8)
  ak = 2*delta_k
  print *,'k=',ak, ' x=',x,' z=', z
  print *,'Gabor0(x, z, ak, sigma)=', Gabor0(x, z, ak, sigma)
  xs = (/0._r8, x, 1._r8/)
  zs =  (/0._r8, z, 1._r8/)
  aks = (/0._r8, ak, 2*ak/)
  sigmas = sigma
  print *,'ks=',aks
  print *,'xs=',xs
  print *,'zs=', zs
  print *,'sigmas=', sigmas
  print *,'GaborN(xs, zs, aks, sigmas)=', GaborN(xs, zs, aks, sigmas)
  
end subroutine gabor_test_run

subroutine kzgrid_test_run
  use kzgrid
  implicit none
  integer ier
  type(kzgrid_obj) :: this_kzgrid  
  
  call kzgrid_test(this_kzgrid, ier)
  if(ier /=0 ) then
     print *,'kzgrid_test **FAILED**'
  else
     print *,'kzgrid_test SUCCEDED'
  endif

end subroutine kzgrid_test_run

subroutine xgrid_test_run
  use xgrid
  implicit none
  integer ier
  type(xgrid_obj) :: this_xgrid  
  
  call xgrid_test(this_xgrid, ier)
  if(ier /=0 ) then
     print *,'xgrid_test **FAILED**'
  else
     print *,'xgrid_test SUCCEDED'
  endif

end subroutine xgrid_test_run
  
subroutine sparseC16_test_run
  use supralu_mod
  implicit none
  integer ier
  type(sparseC16_obj) :: this_sparse

  call sparseC16_test(this_sparse, ier)
  if(ier /=0 ) then
     print *,'sparseC16_test **FAILED**'
  else
     print *,'sparseC16_test SUCCEDED'
  endif
  
end subroutine sparseC16_test_run
