module tam_surface_mod

!-----------------------------------------------------------------------
!  9/13 - tam surface module, written by JML
!    Contains:
!       - tam_surf_temp: updates surface temperature
!       - module variable tam_tgrnd, the ground temperature 
!-----------------------------------------------------------------------
!  10/16 - Added optionally heterogenous surface properties depending on
!          liquid content of grid cell (called from tam_physics)
!          Grid cells with enough liquid are considered lakes with
!          different thermal properties, and also adjust convectively
!-----------------------------------------------------------------------

use fms_mod,          only: error_mesg, FATAL, close_file, mpp_pe, mpp_root_pe, &
                            write_version_number, file_exist, check_nml_error,  &
                            open_namelist_file, stdlog, read_data, write_data, nullify_domain
use transforms_mod,   only: grid_domain 
use diag_manager_mod, only: diag_axis_init, register_diag_field, send_data
use field_manager_mod,only: MODEL_ATMOS
use tracer_manager_mod,only:get_tracer_index
use time_manager_mod, only: time_type
use vert_diff_mod,    only: surf_diff_type
use constants_mod,    only: stefan, cp_air, HLv

implicit none
private

public :: tam_surface_init, tam_surface_end, tam_surf_temp, tam_tgrnd

!-----------------------------------------------------------------------
! Namelist:
!-----------------------------------------------------------------------

  integer :: nlayers       = 10        !Number of ground layers
  real :: cnfac            = .5        !C-N parameter "alpha" (.5 or 1)
  real :: capr             = 1.0 !.34  !Factor for surface depth, 0<.34<1
  real :: thermal_cond_reg = 0.1       !Thermal conductivity (W/m/K), porous icy regolith = .1, rock-ice = 3.5
  real :: thermal_cond_liq = 0.3
  real :: c_v_reg          = 1.12e6    !Volumetric heat capacity [J/m3/K],c_v=c_p*rho=1400*800
  real :: c_v_liq          = 1.8e6     !Volumetric heat capacity [J/m3/J],c_v=c_p*rho=4000*450
  
 
  !Initial surface temperature (cold start)
  real    :: tsurf_init       = 93.0
  
  namelist /tam_surface_nml/ nlayers, cnfac, capr,                &
                               thermal_cond_reg, thermal_cond_liq,  & 
                               c_v_reg, c_v_liq, tsurf_init
                               
  integer  :: id_zgrnd, id_hf, id_tgrnd, sphum                           
  
!-----------------------------------------------------------------------
  character(len=13) :: mod_name = 'tam_surface'
  character(len=*), parameter :: version = '$Id: tam_surface.F90,v 2012/09 JML Exp $'
  character(len=*), parameter :: tagname = '$Name: tam $'
  logical :: module_is_initialized = .false.
  real    :: missing_value = -1.e10
  real, dimension(:,:,:),   allocatable, save  ::  tam_tgrnd
!-----------------------------------------------------------------------

    contains

!-----------------------------------------------------------------------

  subroutine tam_surface_init(lon, lat, axes, Time)
  
   integer, intent(in)              :: axes(4)
   type(time_type), intent(in)      :: Time
   real, dimension(:,:), intent(in) :: lon, lat
  
!-----------------------------------------------------------------------
  
   integer                 :: unit, io, ierr, id, jd, i, j, k, is, js, ie, je
   real,dimension(nlayers) :: z
   logical                 :: used
   character(len=128)      :: filename

!-----------------------------------------------------------------------
! This subroutine is called from tam_physics_init
!-----------------------------------------------------------------------

      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=tam_surface_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'tam_surface_nml')
         enddo
  10     call close_file (unit)
      endif

!-----------------------------------------------------------------------
! Write version info and namelist to log file
!-----------------------------------------------------------------------

  call write_version_number (version,tagname)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=tam_surface_nml)

!-----------------------------------------------------------------------
! Allocate storage for surface properties
!-----------------------------------------------------------------------

  id = size(lon,1)
  jd = size(lat,2)

  allocate (tam_tgrnd(id,jd,nlayers))
  
!-----------------------------------------------------------------------
! Get tracer index for specific humidity
!-----------------------------------------------------------------------
  
  sphum = get_tracer_index( MODEL_ATMOS, 'sphum')

!-----------------------------------------------------------------------
! Read restart file
!-----------------------------------------------------------------------

  filename= 'INPUT/surface.res.nc' 
  if (file_exist( trim( filename ) ) ) then
    call nullify_domain()
    if((mpp_pe() == mpp_root_pe())) print *,'tam_surface module: Reading surface restart file: ',  trim( filename )
    call read_data( trim( filename ), 'tam_tgrnd', tam_tgrnd, grid_domain)
  else
    tam_tgrnd(:,:,:)  = tsurf_init
  end if

!-----------------------------------------------------------------------
! Create subsurface axis ID and register fields
!-----------------------------------------------------------------------

  do k = 1,nlayers
    z(k) = (exp(.5*(real(k)-.5))-1.)
  end do
  id_zgrnd = diag_axis_init('z_ground',z,'m','z','Subsurface level',direction=-1)
  
  id_tgrnd = register_diag_field ( mod_name, 'tam_tgrnd', (/ axes(1), axes(2), id_zgrnd /), &
                   Time, 'Ground temperature (K)', missing_value=missing_value )
  !Note that the surface temperature (i.e. top-level ground temperature) is also
  !output by tam_physics
  
  id_hf = register_diag_field ( mod_name, 'hf', axes(1:2), Time, &
                  'Net flux into the ground', 'W/m2',  &
                   missing_value=missing_value     )

  module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine tam_surface_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  subroutine tam_surf_temp (is, js, dt, Time, liquid, albedoi, &
                        radflx, flux_t, flux_q, dhdt_surf, &
                        dedt_surf, dedq_surf, dhdt_atm,    &
                        dedq_atm, Surf_diff,               &
                        var_surf_prop, threshold_liq, liqch4_dens, &
                        convective_lakes, do_methane_table, &
                        table_height, topo, porosity)
                            
   integer, intent(in)                :: is, js
   real, intent(in)                   :: dt
   type(surf_diff_type), intent(inout):: Surf_diff
   type(time_type),      intent(in)   :: Time
   real, dimension(:,:), intent(in)   :: radflx, flux_t, flux_q
   real, dimension(:,:), intent(in)   :: dhdt_surf, dedt_surf, dedq_surf
   real, dimension(:,:), intent(in)   :: dhdt_atm, dedq_atm
   real, dimension(:,:), intent(in)   :: liquid, albedoi 
   real, dimension(:,:), intent(in)   :: table_height, topo
   logical,              intent(in)   :: var_surf_prop,convective_lakes,do_methane_table
   real,                 intent(in)   :: threshold_liq, liqch4_dens, porosity
   
! Local: --------------------------------------------------------------- 
   integer :: i, j, k, used
   real    :: dzp, dzm
   real, dimension(nlayers) :: at, bt, ct, rt, fact, fn, t_ground
   
   real, dimension(size(radflx,1),size(radflx,2)) :: emiss
   real, dimension(size(radflx,1),size(radflx,2)) :: gamma_t, gamma_q, fn_t, fn_q, &
         en_t, en_q, alpha_t, alpha_q, alpha_lw, beta_t, beta_q, beta_lw, &
         hf, dhfdT
   
   real, dimension(nlayers) :: z  !layer depth (m)
   real, dimension(nlayers) :: zi !layer interface depth
   real, dimension(nlayers) :: dz !layer thickness
   real, dimension(nlayers) :: cv !heat capacity (J/m2/K)
   real, dimension(nlayers) :: tk !thermal conductivity

!-----------------------------------------------------------------------
! Net heat flux into the surface; derivatives of flux wrt temperature  
!-----------------------------------------------------------------------

  emiss(:,:) = 1. - albedoi(:,:)

  gamma_t = 1./(1. - Surf_diff%dtmass * (Surf_diff%dflux_t + dhdt_atm/cp_air))
  gamma_q = 1./(1. - Surf_diff%dtmass * (Surf_diff%dflux_tr(:,:,sphum) + dedq_atm))
  
  fn_t = gamma_t * (Surf_diff%delta_t + Surf_diff%dtmass * flux_t/cp_air)   
  fn_q = gamma_q * (Surf_diff%delta_tr(:,:,sphum) + Surf_diff%dtmass * flux_q)
  
  en_t = gamma_t * Surf_diff%dtmass * dhdt_surf/cp_air
  en_q = gamma_q * Surf_diff%dtmass * dedt_surf
  
  alpha_t = flux_t/cp_air + dhdt_atm * fn_t/cp_air
  alpha_q = flux_q + dedq_atm * fn_q
  alpha_lw= emiss * stefan * tam_tgrnd(:,:,1)**4.
  
  beta_t = dhdt_surf/cp_air + dhdt_atm * en_t/cp_air
  beta_q = dedt_surf + dedq_atm * en_q
  beta_lw= 4. * emiss * stefan * tam_tgrnd(:,:,1)**3.
  
  hf    = radflx - alpha_t*cp_air - alpha_q*HLv - alpha_lw
  dhfdT = -(beta_t*cp_air + beta_q*HLv + beta_lw)
  
!-----------------------------------------------------------------------
! Set up layer depths  
!----------------------------------------------------------------------- 
  
  do k = 1,nlayers
    z(k) = (exp(.5*(real(k)-.5))-1.)
  end do
 
  k = 1
  dz(k) = .5*(z(k)+z(k+1))
  zi(k) = .5*(z(k)+z(k+1))

  do k = 2,nlayers-1 
    dz(k) = .5*(z(k+1)-z(k-1))
    zi(k) = .5*(z(k)+z(k+1))
  end do
  
  k = nlayers
  dz(k) = z(k)-z(k-1)
  zi(k) = z(k)+.5*dz(k)
  
!-----------------------------------------------------------------------
! Enter loops over lat/lon; set up layer properties
!-----------------------------------------------------------------------

  do i = 1,size(tam_tgrnd,1)
    do j = 1,size(tam_tgrnd,2)
    
      do k = 1, nlayers
        cv(k) = dz(k) * c_v_reg
        tk(k) = thermal_cond_reg
        if (var_surf_prop) then
          !if (liquid(i,j) .gt. threshold_liq) cv(k) = dz(k) * c_v_liq
          !if (liquid(i,j) .gt. threshold_liq) tk(k) = thermal_cond_liq 
          
          if (z(k) .lt. liquid(i,j)/liqch4_dens) then 
              cv(k) = dz(k) * c_v_liq
              tk(k) = thermal_cond_liq
          else if (do_methane_table .and. (z(k) .gt. (liquid(i,j)/liqch4_dens + &
                     max(0.0,topo(i,j)-table_height(i,j))) ) ) then 
              cv(k) = dz(k)* (c_v_reg*(1.0-porosity) + c_v_liq*porosity)
              tk(k) = thermal_cond_reg*(1.0-porosity) + thermal_cond_liq*porosity
          end if

        endif
      enddo

      t_ground(:) = tam_tgrnd(i,j,:)
      
!-----------------------------------------------------------------------
! If grid is a "lake," perform convective adjustment  
!-----------------------------------------------------------------------
   
   
      if (convective_lakes) then
        do k = 1, nlayers-1
          if (z(k+1) .lt. liquid(i,j)/liqch4_dens) then
            if (t_ground(k) .lt. t_ground(k+1)) then          
              t_ground(k)   = (t_ground(k)*dz(k) + t_ground(k+1)*dz(k+1))/(dz(k)+dz(k+1))
              t_ground(k+1) = t_ground(k)
            endif
          endif
        enddo
      endif
  
!-----------------------------------------------------------------------
! Set up the matrix elements  
!-----------------------------------------------------------------------
  
      !Factor dt/cv
      k = 1
      zi(0) = 0
      fact(k) = dt / cv(k) * dz(k) / (0.5*(z(k)-zi(k-1)+capr*(z(k+1)-zi(k-1))))
      do k = 2, nlayers
        fact(k) = dt/cv(k)
      enddo
  
      !Conductive Flux
      do k = 1, nlayers - 1
        fn(k) = tk(k)*(t_ground(k+1)-t_ground(k))/(z(k+1)-z(k))
      enddo
      fn(nlayers) = 0.
  
      k = 1
      dzp = z(k+1)-z(k)
      at(k) = 0.
      bt(k) = 1. + fact(k)*((1.-cnfac)*tk(k)/dzp - dhfdT(i,j))
      ct(k) =  -(1.-cnfac)*fact(k)*tk(k)/dzp
      rt(k) = t_ground(k) +  fact(k)*( hf(i,j) - dhfdT(i,j)*t_ground(k) + cnfac*fn(k) )
  
      do k = 2, nlayers - 1
        dzm = (z(k)-z(k-1))
        dzp = (z(k+1)-z(k))

        at(k) =   - (1.-cnfac)*fact(k)* tk(k-1)/dzm
        bt(k) = 1.+ (1.-cnfac)*fact(k)*(tk(k)/dzp + tk(k-1)/dzm)
        ct(k) =   - (1.-cnfac)*fact(k)* tk(k)/dzp

        rt(k) = t_ground(k) + cnfac*fact(k)*( fn(k) - fn(k-1) )
      enddo
  
      k = nlayers
      dzm = (z(k)-z(k-1))
      at(k) =   - (1.-cnfac)*fact(k)*tk(k-1)/dzm
      bt(k) = 1.+ (1.-cnfac)*fact(k)*tk(k-1)/dzm
      ct(k) = 0.
      rt(k) = t_ground(k) - cnfac*fact(k)*fn(k-1)
  
!-----------------------------------------------------------------------
! Solve matrix, update ground temperature  
!-----------------------------------------------------------------------  

      call Tridiagonal (size(at), at, bt, ct, rt, &
                        t_ground(1:nlayers))

      Surf_diff%delta_t(i,j) = fn_t(i,j) + en_t(i,j) * (t_ground(1) - tam_tgrnd(i,j,1))
      Surf_diff%delta_tr(i,j,sphum) = fn_q(i,j) + en_q(i,j) * (t_ground(1) - tam_tgrnd(i,j,1))
      
      tam_tgrnd(i,j,:) = t_ground(:)
      
    end do
  end do
  
  if (id_hf > 0) used = send_data (id_hf, hf, Time, is, js)
  if (id_tgrnd > 0) used = send_data (id_tgrnd, tam_tgrnd, Time, is, js)
 
  end subroutine tam_surf_temp 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine tam_surface_end
  
   character(len=128)    :: filename

    filename = 'RESTART/surface.res.nc'
    call nullify_domain()
    if (mpp_pe() == mpp_root_pe()) print *, 'tam_surface module: Writing to surface restart file: ', &
                                             trim (filename)
    call write_data(trim(filename), 'tam_tgrnd', tam_tgrnd, grid_domain)

    deallocate(tam_tgrnd)
    
    module_is_initialized = .false.

  end subroutine tam_surface_end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  subroutine Tridiagonal (n, a, b, c, r, u )

    integer, intent(in) :: n
    real, intent(in)    :: a(1:n), b(1:n), c(1:n), r(1:n)
    real, intent(out)   :: u(1:n)

    integer :: j
    real    :: gam(1:n)
    real    :: bet

  bet = b(1)
  u(1) = r(1) / bet
  do j = 2, n
     gam(j) = c(j-1) / bet
     bet = b(j) - a(j) * gam(j)
     u(j) = (r(j) - a(j)*u(j-1)) / bet
  enddo

  do j = n-1, 1, -1
     u(j) = u(j) - gam(j+1) * u(j+1)
  enddo

  end subroutine Tridiagonal  
  
end module tam_surface_mod
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
