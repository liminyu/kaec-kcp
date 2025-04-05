          module initial_mod
           contains

            subroutine ini
             use vars
             use vars1
             use global
             use vars_k
             use vars_e
             use inrtype
             use splines
             use supralu_mod
             use orbit_eq_mod
             ! generate mesh_xal_arc.dat of particle mesh
             ! with calling check_mesh3
             !  
!            equilibrium for device configuration
             call check_eq ! using gfile, bkg, rmaj input
!            parameters for device and particles
             call profiledat
!             call profcoefdat
            ! call setgrid
!            flux surface making for gfile
!             call contourpoints
!             call contourpointsr
!             coefficients of Fouier components of R and Z with straight
!            field line coordinates of flux
!             call tensordat
!            flux surface recording
             !call check3
             !call check_mesh
             !call check_mesh1
             call check_flux1
             !call check_flux2
            ! call check_mesh2
            ! call check_mesh3
             call check_mesh5
           end subroutine ini
          end module initial_mod
