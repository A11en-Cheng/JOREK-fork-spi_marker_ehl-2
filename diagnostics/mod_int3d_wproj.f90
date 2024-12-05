!> Calculates 3D integrals and boundary fluxes. including the aux_node_list 
!! only for post_proc
module mod_int3d_wproj

#define NOMPIVERSION=1

  use constants
  use mod_parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use tr_module
  use phys_module
  use mod_interp
  use convert_character
  use mpi_mod
  use mod_expression
  use mod_plasma_functions
  use corr_neg

  implicit none
  
  private
  
  public :: int3d_wproj 
  
  contains


subroutine int3d_wproj(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, aux_node_list, expr_list, res, units)

!$ use omp_lib
 
implicit none

type (type_node_list),        intent(in)    :: node_list
type (type_element_list),     intent(in)    :: element_list   
type (type_bnd_node_list),    intent(in)    :: bnd_node_list
type (type_bnd_element_list), intent(in)    :: bnd_elm_list   
type (type_node_list),        intent(in)    :: aux_node_list
type (t_expr_list),           intent(in)    :: expr_list
real*8,                    intent(inout)    :: res(:)
integer,                      intent(in)    :: units


! --- Local variables
type (type_element)      :: element, elm_k
type (type_node)         :: nodes(n_vertex_max), aux_nodes(n_vertex_max), node_k, aux_node_k
type (type_bnd_element)  :: bndelem
type (type_surface_list) :: surface_list

real*8  :: psi_axis, psi_bnd
real*8  :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss),        x_ss(n_gauss,n_gauss),        x_tt(n_gauss,n_gauss),        x_st(n_gauss,n_gauss)
real*8  :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss),        y_ss(n_gauss,n_gauss),        y_tt(n_gauss,n_gauss),        y_st(n_gauss,n_gauss)
real*8  :: eq_g(n_plane,0:n_var,n_gauss,n_gauss), eq_s(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_t(n_plane,0:n_var,n_gauss,n_gauss), eq_p(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_ss(n_plane,0:n_var,n_gauss,n_gauss), eq_tt(n_plane,0:n_var,n_gauss,n_gauss), eq_st(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_sp(n_plane,0:n_var,n_gauss,n_gauss), eq_tp(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_aux_g(n_plane,0:n_var,n_gauss,n_gauss), eq_aux_s(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_aux_t(n_plane,0:n_var,n_gauss,n_gauss), eq_aux_p(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_spp(n_plane,0:n_var,n_gauss,n_gauss), eq_tpp(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: wgauss_copy(n_gauss)
real*8  :: psi_axisym(n_gauss,n_gauss)

real*8  :: xt, t_norm, rho_norm, t_norm2
integer :: i, j, k, in, ms, mt, mp, iv, inode, ife, n_elements, i_elm_axis, i_elm_xpoint(2), ifail, i_phi
integer :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
integer :: k_vertex, k_dof, k_node, k_dir, k_dir_perp, m_bndelem, dir_perp(2), mv1, m_elm

real*8  :: xjac, xjac_R, xjac_Z, BigR, wst, delta_phi, phi
real*8  :: rimp0, rimp0_corr
real*8  :: local_radiation, total_radiation
real*8  :: local_radiation_phi(n_plane), total_radiation_phi(n_plane), rad_TPF_phi(n_plane)
real*8  :: local_nimp_phi(n_plane), total_nimp_phi(n_plane), nimp_TPF_phi(n_plane)
real*8  :: m_i_over_m_imp, m_imp

#ifndef NOMPIVERSION
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI procs
n_cpu = max(n_cpu,1)
#else
n_cpu = 1
#endif

local_radiation       = 0.d0
local_radiation_phi   = 0.d0
local_nimp_phi        = 0.d0
delta_phi     = 2.d0 * PI / float(n_plane) / float(n_period)
wgauss_copy = wgauss

ife_delta = ceiling(float(element_list%n_elements) / n_cpu)
ife_min   =      my_id     * ife_delta + 1
ife_max   = min((my_id +1) * ife_delta, element_list%n_elements)

!$omp parallel default(none)                                                                   &
!$omp   shared(element_list, node_list, H, H_s, H_t, HZ, HZ_p, ife_min, ife_max, xpoint,       &
!$omp          H_ss, H_tt, H_st, HZ_pp, aux_node_list, use_marker, my_id, delta_phi, GAMMA,    &
!$omp          local_radiation, local_radiation_phi, imp_type, local_nimp_phi, index_main_imp, &
!$omp          wgauss_copy, central_density, central_mass)                                                                    &
!$omp   private(ife,iv,inode,element,nodes,aux_nodes,i,j, k,in, mp, ms, mt, phi,          &
!$omp           x_g, y_g, x_s, y_s, x_t, y_t, xjac, xjac_R, xjac_Z, eq_g, eq_s, eq_t, eq_p,    &
!$omp           x_ss, x_tt, x_st, y_ss, y_tt, y_st, eq_ss, eq_tt, eq_st, eq_sp, eq_tp,         &
!$omp           eq_spp, eq_tpp, psi_axisym, eq_aux_g, eq_aux_s, eq_aux_t, eq_aux_p,            &
!$omp           wst, BigR, rimp0, rimp0_corr, m_i_over_m_imp, m_imp, omp_nthreads,omp_tid)


#ifdef OPENMP
omp_nthreads = omp_get_num_threads()
omp_tid      = omp_get_thread_num()
#else
omp_nthreads = 1
omp_tid      = 0
#endif

!$omp do reduction(+:local_radiation, local_radiation_phi, local_nimp_phi)

do ife = ife_min, ife_max

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode         = element%vertex(iv)
    nodes(iv)     = node_list%node(inode)
    if (use_marker) aux_nodes(iv) = aux_node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0; x_ss(:,:)    = 0.d0; x_tt(:,:)    = 0.d0; x_st(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0; y_ss(:,:)    = 0.d0; y_tt(:,:)    = 0.d0; y_st(:,:)    = 0.d0;
  psi_axisym(:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_degrees

      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)
          x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)
          y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
          y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
          y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)


        enddo
      enddo
    enddo
  enddo

  eq_g(:,:,:,:) = 0.d0; eq_s(:,:,:,:) = 0.d0; eq_t(:,:,:,:) = 0.d0; eq_p(:,:,:,:) = 0.d0; eq_ss(:,:,:,:) = 0.d0; eq_tt(:,:,:,:) = 0.d0; eq_st(:,:,:,:) = 0.d0; 
  eq_sp(:,:,:,:) = 0.d0; eq_tp(:,:,:,:) = 0.d0; eq_spp(:,:,:,:) = 0.d0; eq_tpp(:,:,:,:) = 0.d0;
  eq_aux_g(:,:,:,:) = 0.d0; eq_aux_s(:,:,:,:) = 0.d0; eq_aux_t(:,:,:,:) = 0.d0; eq_aux_p(:,:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_degrees

      do mp=1,n_plane
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              do in=1,n_tor
                eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
                eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
                eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
                eq_sp(mp,k,ms,mt) = eq_sp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ_p(in,mp)
                eq_tp(mp,k,ms,mt) = eq_tp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ_p(in,mp)
                eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
                eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)
                eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
                if (use_marker) then
                  eq_aux_g(mp,k,ms,mt) = eq_aux_g(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                  eq_aux_s(mp,k,ms,mt) = eq_aux_s(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
                  eq_aux_t(mp,k,ms,mt) = eq_aux_t(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
                  eq_aux_p(mp,k,ms,mt) = eq_aux_p(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
                endif

                eq_spp(mp,k,ms,mt) = eq_spp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)*HZ_pp(in,mp)
                eq_tpp(mp,k,ms,mt) = eq_tpp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)*HZ_pp(in,mp)
              enddo
            enddo

          enddo
        enddo
      enddo

    enddo
  enddo
  

  !--------------------------------------------------- sum over the Gaussian integration points
  do mp=1,n_plane

    phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)

    do ms=1, n_gauss
      do mt=1, n_gauss

        wst  = wgauss_copy(ms)*wgauss_copy(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        xjac_R = ( x_ss(ms,mt) * y_t(ms,mt)**2 - 2*x_st(ms,mt) * y_s(ms,mt)*y_t(ms,mt) + x_tt(ms,mt) * y_s(ms,mt)**2 + x_s(ms,mt) * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )   &
                 + x_t(ms,mt) * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac
        xjac_Z = ( y_ss(ms,mt) * x_t(ms,mt)**2 - 2*y_st(ms,mt) * x_s(ms,mt)*x_t(ms,mt) + y_tt(ms,mt) * x_s(ms,mt)**2 + y_s(ms,mt) * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )   &
                 + y_t(ms,mt) * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac
        BigR = x_g(ms,mt)

        rimp0  = eq_g(mp,var_rhoimp,ms,mt)
        rimp0_corr = corr_neg_dens(rimp0, (/ 0.d-5, 1.d-5 /)) ! Correction for negative rimp0
      

!-------------------------------------------
! --- Radiation and ionization power
! ------------------------------------------
        !-------------------------------------------
        ! Atomic physics parameters for Impurities
        !-------------------------------------------

        select case ( trim(imp_type(index_main_imp)) )
          case('D2')
            m_i_over_m_imp = central_mass/2.  ! Deuterium mass = 2 u
            m_imp          = 2.
          case('Ar')
            m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u
            m_imp          = 40.
          case('Ne')
            m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u
            m_imp          = 20.
          case('Fe')
            m_i_over_m_imp = central_mass/56. ! Iron mass = 56 u
          case('W')
            m_i_over_m_imp = central_mass/184. ! Tungsten mass = 184 u
          case default
            write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown (in mod_injection_source.f90) !!'
            write(*,*) '=> We assume the gas is D2.'
            m_i_over_m_imp = central_mass/2.
            m_imp          = 2.
        end select

        if (use_marker) then
            local_nimp_phi(mp) = local_nimp_phi(mp) +  (rimp0_corr*m_i_over_m_imp*central_density) * bigr * xjac * wst * delta_phi  
            local_radiation_phi(mp) = local_radiation_phi(mp) +  (eq_aux_g(mp,2,ms,mt)/(MU_ZERO * (GAMMA-1.) * sqrt(MU_ZERO * CENTRAL_MASS * MASS_PROTON * CENTRAL_DENSITY * 1.d20))) * bigr * xjac * wst * delta_phi  
            local_radiation         = local_radiation         +  (eq_aux_g(mp,2,ms,mt)/(MU_ZERO * (GAMMA-1.) * sqrt(MU_ZERO * CENTRAL_MASS * MASS_PROTON * CENTRAL_DENSITY * 1.d20))) * bigr * xjac * wst * delta_phi 
        endif

      enddo
    enddo
  enddo
enddo
!$omp end do
!$omp end parallel


! --- gather contribution from all MPI processes
#ifndef NOMPIVERSION
call MPI_AllReduce(local_radiation, total_radiation,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(local_radiation_phi, total_radiation_phi,n_plane,&
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(local_nimp_phi, total_nimp_phi,n_plane,&
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#else /* NOMPIVERSION */
total_radiation      = local_radiation
total_radiation_phi  = local_radiation_phi
total_nimp_phi       = local_nimp_phi
#endif /* NOMPIVERSION */
total_radiation     = n_period * total_radiation
total_radiation_phi = n_period * total_radiation_phi
total_nimp_phi      = n_period * total_nimp_phi



if (my_id .eq. 0) then 

  ! --- Get current time
  if (index_now > 0) then
    xt = xtime(index_now)
  else
    xt = 0.d0
  endif

  if (use_marker) then          !output radiation at each plane  TODO: SI units
    do i_phi = 1, n_plane
      rad_TPF_phi(i_phi) = total_radiation_phi(i_phi)/(sum(total_radiation_phi)/n_plane)
      nimp_TPF_phi(i_phi) = total_nimp_phi(i_phi)/(sum(total_nimp_phi)/n_plane)
    end do
    open(20,file="total_radiation_phi.dat",action="write",position="append")
    write(20,*) '# time', t_now
    do i_phi = 1, n_plane
      write(20,'(4e14.6)') total_radiation_phi(i_phi)/1.d6, rad_TPF_phi(i_phi), total_nimp_phi(i_phi), nimp_TPF_phi(i_phi)
    end do
    close (20)

    open(20,file="total_radiation_vs_rad_TPF.dat",action="write",position="append")
    write(20,'(i14,4e14.6)') index_now, xt, maxval(rad_TPF_phi), total_radiation, maxval(nimp_TPF_phi) 
    close (20)
  end if

endif !--- my_id
end subroutine int3d_wproj 

#undef NOMPIVERSION

end module mod_int3d_wproj
