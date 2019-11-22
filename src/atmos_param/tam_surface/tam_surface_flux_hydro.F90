 module tam_surface_flux_mod

!-----------------------------------------------------------------------
!  9/13 - surf_flux for tam. JML
!-----------------------------------------------------------------------

use fms_mod,          only: error_mesg, FATAL, close_file, mpp_pe, mpp_root_pe, &
                            write_version_number, file_exist, check_nml_error,  &
                            open_namelist_file, stdlog, read_data, write_data, nullify_domain
use time_manager_mod, only: time_type

use monin_obukhov_mod, only: mo_drag, mo_profile
use constants_mod,     only: cp_air, rdgas, rvgas, grav, kappa, stefan, hlv
use sat_vapor_pres_mod, only: lookup_es

implicit none
private
  
public :: tam_surf_flux_2d

!-----------------------------------------------------------------------
! Namelist:
!-----------------------------------------------------------------------

  logical :: do_mo_drag       = .true.
  
  !Fluxes options
  real    :: gust_min         = 0.1
  real    :: cd_drag_cnst     = 0.003
  real    :: sfc_heat_flx_amp = 1.0
  
  namelist /tam_surface_nml/ do_mo_drag, gust_min, &
                               cd_drag_cnst, sfc_heat_flx_amp
  
  !-----------------------------------------------------------------------

  character(len=13) :: mod_name = 'tam_surface_flux'
  real    :: missing_value = -1.e10

!-----------------------------------------------------------------------

 character(len=*), parameter :: version = '$Id: surface_flux.F90,v 2013/09 JML Exp $'
 character(len=*), parameter :: tagname = '$Name: tam $'


!-----------------------------------------------------------------------

    contains

!-----------------------------------------------------------------------
   subroutine surf_flux_1d (t_atm, q_atm, u_atm, v_atm, p_atm, z_atm,   &
     p_surf, t_surf, q_surf_in, rough_mom, rough_heat, rough_moist, gust,  &
     flux_t, flux_q, flux_u, flux_v, cd_m, cd_t, cd_q,                  &
     w_atm, u_star, b_star, q_star,                                     &
     dhdt_surf, dedt_surf, dedq_surf,                                   &
     dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm,                        &
     is, js, dt, Time, do_gwe, evap_thresh)

  integer, intent(in)                      :: is, js
  real, intent(in)                         :: dt, evap_thresh
  type(time_type), intent(in)              :: Time

  real, intent(in), dimension(:) :: t_atm, q_atm, u_atm, v_atm,  &
     p_atm, z_atm, p_surf, t_surf, rough_mom, rough_heat, rough_moist, gust

  real, intent(out), dimension(:) :: flux_t, flux_q, flux_u, flux_v,    &
       cd_m, cd_t, cd_q,                                                &
       dhdt_surf, dedt_surf, dedq_surf,                                 &
       dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm,                      &
       w_atm, u_star, b_star, q_star

  real, intent(in), dimension(:) :: q_surf_in
  logical, intent(in) :: do_gwe

!local:

  real, dimension(size(t_atm(:))) :: thv_atm, th_atm, tv_atm,  thv_surf,  &
       q_sat, q_sat1, p_sat, p_sat1, p_ratio, u_diff, v_diff,  &
       rho_drag, drag_t, drag_m, drag_q, rho, dw_atmdu, dw_atmdv, w_gust, &
       u_surf,   v_surf

  real, dimension(size(q_surf_in)) :: q_surf
  
  real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp

  integer :: i

   
!-----------------------------------------------------------------------
! Saturation specific humidity
!-----------------------------------------------------------------------

  !Sat vapor pressure (in Pa), perturbed sat vapor pressure
  call lookup_es(t_surf,p_sat)
  call lookup_es(t_surf+del_temp,p_sat1)
  
  q_sat  = Rdgas/Rvgas * (p_sat / (p_surf - (1.0-(Rdgas/Rvgas)*p_sat)))
  q_sat1 = Rdgas/Rvgas * (p_sat1 / (p_surf - (1.0-(Rdgas/Rvgas)*p_sat1)))
  
  where (q_surf_in .gt. evap_thresh) 
    q_surf = q_sat
  elsewhere
    q_surf = q_atm !q_surf_in/100. * q_sat 
  end where 


!-----------------------------------------------------------------------
! Generate data needed by Monin-Obukhov
!-----------------------------------------------------------------------

  p_ratio = (p_surf/p_atm)**kappa
  
  tv_atm = t_atm * (1.+(1.-Rdgas/Rvgas)/(Rdgas/Rvgas)*q_atm) !virtual temperature
  th_atm = t_atm * p_ratio                                   !potential temperature
  thv_atm = tv_atm * p_ratio                                 !virtual potential temperature
  thv_surf= t_surf*(1.+(1.-Rdgas/Rvgas)/(Rdgas/Rvgas)*q_surf)!surface virtual temperature

!-----------------------------------------------------------------------
! Velocity components relative to surface
!-----------------------------------------------------------------------

  u_surf = 0.0
  v_surf = 0.0
  u_diff = u_surf - u_atm
  v_diff = v_surf - v_atm

  w_gust = max(gust,gust_min)
  w_atm = sqrt(u_diff*u_diff + v_diff*v_diff + w_gust*w_gust)

  dw_atmdu = u_diff/w_atm
  dw_atmdv = v_diff/w_atm

!-----------------------------------------------------------------------
! Monin-Obukhov similarity theory
!-----------------------------------------------------------------------

  if (do_mo_drag) then
    call mo_drag (thv_atm, thv_surf, z_atm, rough_mom, rough_heat,  &
                  rough_moist, w_atm, cd_m, cd_t, cd_q, u_star, b_star)
  else
    cd_m = cd_drag_cnst
    cd_t = cd_drag_cnst
    cd_q = cd_drag_cnst
    u_star = 0.0
    b_star = 0.0
  end if

!-----------------------------------------------------------------------
! Surface layer drag coefficients
!-----------------------------------------------------------------------

  drag_t = cd_t * w_atm
  drag_q = cd_q * w_atm
  drag_m = cd_m * w_atm

!-----------------------------------------------------------------------
! Fluxes
!-----------------------------------------------------------------------

  rho = p_atm/(rdgas * tv_atm)

  rho_drag = cp_air * drag_t * rho * sfc_heat_flx_amp
  flux_t = rho_drag * (t_surf - th_atm)      !Sensible heat (W/m**2)
  dhdt_surf = rho_drag                       !d(sensible_heat_flux)/d(surf_temp)
  dhdt_atm = -rho_drag * p_ratio             !d(sensible_heat_flux)/d(atm_temp)

  rho_drag = drag_q * rho
  flux_q = rho_drag * (q_surf - q_atm)       !Methane vapor flux (kg/(m**2 s))

  where ((flux_q*dt/2. .gt. (q_surf_in)) .and. (q_surf_in .gt. evap_thresh))
    flux_q = (q_surf_in - evap_thresh)*2./dt !SPF
  endwhere

  if (do_gwe) then
    where (q_surf_in .le. evap_thresh)
      flux_q = rho_drag * (q_sat - q_atm)
    endwhere
  end if


  !where (q_surf_in .le. 100.)
    !dedq_surf = 0.0                          !d(latent_heat_flux)/d(surf_mixing_ratio)
    !dedt_surf = 0.0                          !d(latent_heat_flux)/d(surf_temp)
    !dedq_atm  = 0.0                          !d(latent_heat_flux)/d(atm_mixing_ratio)
  !elsewhere
    dedq_surf = 0.0
    dedt_surf = rho_drag * (q_sat1 - q_sat) * del_temp_inv
    dedq_atm  = -rho_drag
  !end where
 
  q_star = flux_q / (u_star * rho)           !moisture scale
  
!-----------------------------------------------------------------------
! Stresses
!-----------------------------------------------------------------------

  rho_drag = drag_m * rho
  flux_u = rho_drag * u_diff                 !zonal stress component (Nt/m**2)
  flux_v = rho_drag * v_diff                 !meridional stress component

  dtaudu_atm = -cd_m * rho * (dw_atmdu * u_diff + w_atm)
  dtaudv_atm = -cd_m * rho * (dw_atmdv * v_diff + w_atm)

!-----------------------------------------------------------------------
  end subroutine surf_flux_1d
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine tam_surf_flux_2d (t_atm, q_atm, u_atm, v_atm, p_atm, z_atm,     &
     p_surf, t_surf, q_surf, rough_mom, rough_heat, rough_moist, gust,  &
     flux_t, flux_q, flux_u, flux_v, cd_m, cd_t, cd_q,                  &
     w_atm, u_star, b_star, q_star,                                     &
     dhdt_surf, dedt_surf, dedq_surf,                                   &
     dhdt_atm, dedq_atm, dtaudu_atm, dtaudv_atm,                        &
     is, js, dt, Time, do_gwe, evap_thresh)

  ! ---- arguments -----------------------------------------------------------

  integer, intent(in)                     :: is, js
  real, intent(in)                        :: dt, evap_thresh
  type(time_type), intent(in)             :: Time

  real, intent(in),  dimension(:,:) ::                        &
       t_atm,     q_atm,   u_atm,     v_atm,                  &
       p_atm,     z_atm,                                      &
       p_surf,    t_surf,                                     &
       rough_mom, rough_heat, rough_moist, gust

  real, intent(out), dimension(:,:) ::                       &
       flux_t,    flux_q,     flux_u,  flux_v,               &
       cd_m,      cd_t,       cd_q,                          &
       dhdt_surf, dedt_surf,  dedq_surf,                     &
       dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,        &
       w_atm,     u_star,     b_star,    q_star

  real, intent(inout), dimension(:,:) :: q_surf
  logical, intent(in)                :: do_gwe


  ! ---- local vars -----------------------------------------------------------
  integer :: j

  do j = 1, size(t_atm,2)
     call surf_flux_1d (t_atm(:,j), q_atm(:,j), u_atm(:,j), v_atm(:,j),   &
          p_atm(:,j), z_atm(:,j), p_surf(:,j), t_surf(:,j), q_surf(:,j),  &
          rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), gust(:,j),   &
          flux_t(:,j), flux_q(:,j), flux_u(:,j), flux_v(:,j),             &
          cd_m(:,j), cd_t(:,j), cd_q(:,j), w_atm(:,j),                    &
          u_star(:,j), b_star(:,j), q_star(:,j),                          &
          dhdt_surf(:,j), dedt_surf(:,j), dedq_surf(:,j),                 &
          dhdt_atm(:,j), dedq_atm(:,j), dtaudu_atm(:,j), dtaudv_atm(:,j), &
          is, js, dt, Time, do_gwe, evap_thresh)
  end do

  end subroutine tam_surf_flux_2d
  
end module tam_surface_flux_mod  
