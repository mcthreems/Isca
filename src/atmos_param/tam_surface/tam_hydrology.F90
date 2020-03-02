module tam_hydrology_mod
!-----------------------------------------------------------------------
! tam hydrology module - SPF 2016/2017
!	Called from tam_physics_mod
!	Computes surface runoff, infiltration, and methane table adjustments
!-----------------------------------------------------------------------

use constants_mod,              only: pi, Radius, grav, liq_dens
use fms_mod,                    only: mpp_npes, error_mesg, FATAL, NOTE, file_exist,   &
                                      open_namelist_file, check_nml_error, &
                                      mpp_pe, mpp_root_pe, close_file,     &
                                      write_version_number, stdlog,        &
                                      read_data, write_data, field_size, nullify_domain
use transforms_mod,             only: grid_domain, area_weighted_global_mean, &
                                      trans_grid_to_spherical, trans_spherical_to_grid, &
                                      get_spec_domain, get_eigen_laplacian,  &
                                      get_num_fourier, get_num_spherical, get_grid_domain, &
                                      get_grid_boundaries
use mpp_domains_mod,            only: mpp_get_global_domain, mpp_define_layout, mpp_define_domains, &
                                      mpp_update_domains, mpp_get_compute_domain, domain2D, mpp_get_data_domain, &
                                      mpp_global_field,mpp_global_min
use topog_regularization_mod,   only: compute_lambda,regularize
use topography_mod,             only: get_topog_mean
use spherical_fourier_mod,      only: get_wts_lat,compute_gradient_cos
use horiz_interp_mod,           only: horiz_interp
                                      
implicit none
private

public :: tam_hydrology_init, tam_hydrology_end, tam_hydrology_driver

type(domain2D), save, public :: topo_domain

logical :: module_initialized = .false.
logical  :: first_hydro_call    = .true.

integer :: is,ie,isd,ied,js,je,jsd,jed,pe,npes, nx, ny
integer :: ms, me, ns, ne 

real, dimension(:,:), allocatable, save  ::  dtd,lat_halo

!-----------------------------------------------------------------------
! Namelist:
!-----------------------------------------------------------------------

  !Methane density:
  !Note that, since this is a parameter, it is not in the actual namelist
  !real,parameter :: liq_dens   = 450.    !Density of liquid methane (km/m3)
 
! Topography
  logical :: do_flat_topo = .true.      !Use a default flat topography
  character(len=64) :: topog_file    = 'topo.nc'   !Otherwise, use input file

! Surface runoff
  logical :: do_topo_runoff  = .false.  !Do runoff at all? 
  logical :: input_speed     = .false.  !Use 'run_speed' parameter below?
  logical :: linear_catch    = .true.  !Calculate runoff speed based on linear reservoir approximation, using concentration time
  real :: basin_length       = 150000.     !Median basin length from Horvath et al. 2016
  real :: run_speed          = 1e-8       !Runoff speed (m/s) = translates to 4.5e-6 kg/m2/s if input_speed is called
 
! Infiltration 
  logical :: do_global_infil = .false.  !Do "infiltration," removing surface liquid everywhere?
  logical :: infil_first     = .true.  !Compute infiltration before runoff
  real    :: hyd_cond        = 1.0e-8  !Horizontal hydraulic conductivity (m/s)
  real    :: infil_rate      = 1.0e-6  !Vertical hydraulic conductivity (m/s)
  real    :: infil_thresh    = 10.0     !Threshold under which there is no infiltration (kg/m2)
 
! Subsurface flow
  logical :: do_darcy        = .true.  !Do fractional descent flow of methane table
  logical :: do_diff         = .false. !Do real diffusion of methane table
  
! Other options
  logical :: init_lakes       = .false.  !Initialize surface liquid to observed maria?
  logical :: init_gauss       = .false.   !Initialize gaussian methane table distribution?
  logical :: do_leapfrog      = .false.    !THIS IS UNTESTED, would recommend to keep false 
  real     :: robert_surf_liq     = 0.04
  real     :: RAW_surf_liq        = 0.53


  namelist /hydrology_nml/ input_speed,linear_catch,run_speed,infil_first,   &
                           basin_length,hyd_cond,infil_rate,  &
                           do_diff,do_darcy,do_flat_topo, &
                           topog_file,do_topo_runoff,do_global_infil,infil_thresh, &
                           init_lakes,init_gauss,do_leapfrog, &
                           robert_surf_liq,RAW_surf_liq               
                               
  
!-----------------------------------------------------------------------
  character(len=13) :: mod_name = 'tam_hydrology'
  character(len=*), parameter :: version = '$Id: hydrology.F90,v 2017/01 SPF Exp $'
  character(len=*), parameter :: tagname = '$Name: tam $'
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine tam_hydrology_init(zsurf,topo_out,liq_begin,table_begin,height_begin,init_surf,init_table,porosity,lat,lon)
!-----------------------------------------------------------------------
! Initialize surface liquids based on low topography (called in tam_physics)
!-----------------------------------------------------------------------

  real, intent(in),    dimension(:,:)   :: zsurf,lat,lon
  real, intent(inout),  dimension(:,:)  :: topo_out,liq_begin,table_begin,height_begin
  real, intent(in)   :: init_surf,init_table,porosity

  real, dimension(size(zsurf,1),size(zsurf,2)) :: surf_height,surf_geopotential,init_depth 
  integer :: i, j, global_num_lon, global_num_lat
  integer, dimension(4) :: siz
  complex, allocatable, dimension(:,:) :: spec_tmp

  real, allocatable, dimension(:) :: blon,blat
  logical :: topo_file_exists

  !For gaussian option
  real    :: olon,olat,wlon,wlat,rlon,rlat
  real    :: tpi,dtr,dx,dy,xx,yy

  integer :: layout(2)
 
  integer                 :: unit, io, ierr
!-----------------------------------------------------------------------
! Write namelist variables
!-----------------------------------------------------------------------


  if (file_exist('input.nml')) then
    unit = open_namelist_file ( )
    ierr = 1; do while (ierr /= 0)
    read (unit, nml = hydrology_nml, iostat = io, end = 10)
    ierr = check_nml_error (io, 'hydrology_nml')
    end do
    10 call close_file (unit)
  end if
  call write_version_number (version,tagname)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=hydrology_nml)


!-----------------------------------------------------------------------
! Open topo file to set up data domain 
!-----------------------------------------------------------------------

 if (.not. do_flat_topo) then
	 if(file_exist('INPUT/'//trim(topog_file))) then
	   call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat) 
	   call field_size('INPUT/'//trim(topog_file), 'zsurf', siz)
	   if ( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) then
		 call read_data('INPUT/'//trim(topog_file), 'zsurf', surf_height, grid_domain)
	   else
		 call error_mesg ('tam_hydrology_init', &
						  'Topography file contains data on wrong grid', FATAL)
	   endif

	   !Ensure all topo heights are positive to make methane table comparisons easier 
	   if (mpp_global_min(grid_domain,surf_height) .le. 0.0) surf_height = surf_height + &
			 abs(mpp_global_min(grid_domain,surf_height)) 
   
	   !Spectrally truncate the topography
	   call get_spec_domain(ms, me, ns, ne)
	   allocate(spec_tmp(ms:me, ns:ne))
	   call trans_grid_to_spherical(surf_height,spec_tmp)
	   call trans_spherical_to_grid(spec_tmp,surf_height)
	   deallocate(spec_tmp)

	 else
	   call error_mesg('tam_hydrology_init','INPUT/'//trim(topog_file)//' does not exist', FATAL)
	 end if
else
	 ! Get flat topography
	   call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat) 
	   call get_grid_domain(is, ie, js, je)
	   allocate(blon(is:ie+1), blat(js:je+1))
	   call get_grid_boundaries(blon, blat)
	   ! topo_file_exists = get_topog_mean(blon, blat, surf_height)
! 	   if (.not.topo_file_exists) then 
! 		 call error_mesg('tam_hydrology_init','topography data file does not exist', FATAL)
! 	   endif
	   surf_height = 1.0					!mmm just setting a flat topography with all points at height 1
	   surf_geopotential = grav*surf_height

	!    Spectrally truncate the flat topography 
		call get_spec_domain(ms, me, ns, ne)
		allocate(spec_tmp(ms:me, ns:ne))
		call trans_grid_to_spherical(surf_geopotential,spec_tmp)
		call trans_spherical_to_grid(spec_tmp,surf_geopotential)
		deallocate(spec_tmp)
		deallocate(blon, blat)
		surf_height = surf_geopotential/grav
		if (mpp_global_min(grid_domain,surf_height) .lt. 0.0) surf_height = &
			 surf_height + abs(mpp_global_min(grid_domain,surf_height))
end if

!-----------------------------------------------------------------------
!  Initialize surface and methane table heights
!-----------------------------------------------------------------------

 liq_begin = 0.0; table_begin = 0.0
 topo_out = surf_height  !output topography that is used for calculating runoff, etc.

 if (init_lakes) then  !Impose very rough approximation of observed seas
      
      dtr=2.0*pi/360.
      init_depth = 0.0

      do j=1,size(lat,2)
        do i=1,size(lon,1)
          if ((lon(i,1) .le. 130.*dtr) .and. (lat(1,j) .ge. 50.*dtr)) then
            init_depth(i,j) = 750. - surf_height(i,j) 
          else if ((lat(1,j) .le. -70.*dtr) .and. (lat(1,j) .ge. -75.*dtr) .and. (lon(i,1) .ge. 165.*dtr) .and. (lon(i,1) .le. 185.*dtr)) then
            init_depth(i,j) = 800. - surf_height(i,j)
          end if
        end do
      end do
      where (init_depth .le. 0.0) init_depth = init_surf
      liq_begin = init_depth*liq_dens 
 
      table_begin = init_table*liq_dens*porosity
      height_begin = table_begin/liq_dens/porosity
      if (height_begin(i,j) .ge. surf_height(i,j)) &
           table_begin(i,j) = surf_height(i,j)*liq_dens*porosity


 else if (init_gauss) then  !Note these gaussian parameters are hard-coded
  tpi = 2.0*pi
  dtr = tpi/360.
  olon = 180.*dtr
  olat = 0.*dtr
  wlon = 30.*dtr
  wlat = 30.*dtr
  rlon = 0.
  rlat = 0.
  !compute gaussian
    do j=1,size(lat,2)
      dy = abs(lat(1,j) - olat)   ! dist from y origin
      yy = max(0., dy-rlat)/wlat
      do i=1,size(lon,1)
        dx = abs(lon(i,1) - olon) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0., dx-rlon)/wlon
        table_begin(i,j) = init_table*liq_dens*porosity*exp(-xx**2 - yy**2)
        !table_begin(i,j) = init_table*liq_dens*porosity
        height_begin(i,j) = table_begin(i,j)/liq_dens/porosity
        if (height_begin(i,j) .ge. surf_height(i,j)) then
           table_begin(i,j) = surf_height(i,j)*liq_dens*porosity
           liq_begin(i,j) = (height_begin(i,j) - surf_height(i,j))*liq_dens + init_surf*liq_dens
        else
           liq_begin(i,j) = init_surf*liq_dens
        end if
      enddo
    enddo
   else   !DEFAULT INITIALIZATION: constant equipotential = init_table; 
          ! and surface liquid is everywhere = init_surf
    call error_mesg('tam_hydrology','Default surface and subsurface initialization', NOTE)
     do j=1,size(lat,2)
      do i=1,size(lon,1)
        !if((j == 16 .or. j==17) .and. (i == 32)) then 
         !   print *,'making two cells of different heights' 
            table_begin(i,j) = (init_table)*liq_dens*porosity
            !table_begin(i,j) = (surf_height(i,j))*liq_dens*porosity
            height_begin(i,j) = table_begin(i,j)/liq_dens/porosity
            if (height_begin(i,j) .gt. surf_height(i,j)) then
             table_begin(i,j) = surf_height(i,j)*liq_dens*porosity
             liq_begin(i,j) = (height_begin(i,j) - surf_height(i,j))*liq_dens   
            else
             liq_begin(i,j) = init_surf*liq_dens
            end if 
            liq_begin(i,j) = init_surf*liq_dens   
         !end if
      enddo
    enddo

   end if

!-----------------------------------------------------------------------
! Set up data domain using 1 halo point in latitude 
!-----------------------------------------------------------------------
 
  nx = global_num_lon
  ny = global_num_lat

  pe = mpp_pe()
  npes = mpp_npes()

 
  !2D decomposition along X and Y 
 ! call mpp_define_layout( (/1,nx,1,ny/),npes,layout)
 ! call mpp_define_domains( (/1,nx,1,ny/), layout, topo_domain,xhalo=1,yhalo=1 )
  
  !1D decomposition along Y
  layout = (/1,npes/)
  call mpp_define_domains( (/1,nx,1,ny/), layout, topo_domain, yhalo=1)

  call mpp_get_data_domain( topo_domain, isd, ied, jsd, jed) !isd = is = 0, ied = ie = 64; jsd = js-1, jed = je+1
  call mpp_get_compute_domain( topo_domain, is, ie, js, je )
  allocate( dtd(is:ie,jsd:jed) )  !dtd = Data Topography Domain 
  allocate( lat_halo(is:ie,jsd:jed)) 
  dtd = 0.0; lat_halo = 0.0
  do i = is,ie
    do j = js,je
      dtd(i,j) = surf_height(i-is+1,j-js+1)
      lat_halo(i,j) = lat(i-is+1,j-js+1)
    end do
  end do

  !Fill in halo points
  call mpp_update_domains( dtd, topo_domain)       !now dtd(:jsd) and dtd(:,jed) are filled in with appropriate values
  call mpp_update_domains( lat_halo, topo_domain)

  !if (js .eq. 1) then 
    print *,'global min (dtd): ',mpp_global_min(grid_domain,dtd(:,js:je))
  !endif 

  if ((linear_catch) .and. (input_speed)) &
      call error_mesg('hydrology', &
                      'linear_catch and input_speed cannot both be .true.', FATAL)
 
  if ((do_diff) .and. (do_darcy)) &
      call error_mesg('hydrology', &
                      'do_diff and do_darcy cannot both be .true.', FATAL)
  
  if ((init_lakes) .and. (init_gauss)) &
      call error_mesg('hydrology', &
                      'init_lakes and init_gauss cannot both be .true.', FATAL)
 

  !hyd_diff = hyd_cond/(liq_dens*grav*beta_comp*porosity)

  module_initialized = .true.
  if (mpp_pe()==mpp_root_pe()) print *, 'Hydrology initialized.'
  
end subroutine tam_hydrology_init

subroutine tam_hydrology_end

   deallocate(dtd)
   deallocate(lat_halo)


   module_initialized = .false.

end subroutine tam_hydrology_end


subroutine tam_hydrology_driver(surf_liq,run_res,meth_table,height,runoff,infiltration,gwe,run_out,dis_out,recharge,&
                            subflow,tot_rain,evap,porosity,lat,lon,current,previous,future,delta_t,dt,do_gwe,   &
                            do_methane_table,evap_thresh)

  real, intent(inout),  dimension(:,:,:)  :: surf_liq,run_res,meth_table,height
  real, intent(inout),  dimension(:,:  )  :: runoff,infiltration,gwe,run_out,dis_out,recharge,subflow
  real, intent(in),     dimension(:,:  )  :: tot_rain,evap,lat,lon
  real, intent(in)                        :: delta_t,dt,porosity,evap_thresh
  logical, intent(in)                     :: do_gwe,do_methane_table

  integer, intent(in)                     :: current,previous
  integer, intent(inout)                  :: future
  
  real, dimension(ny) :: wts_lat   !latitude weights
  real, dimension(size(surf_liq,1),size(surf_liq,2))     :: filt,surf_liq_dt,meth_table_dt,height_dt,diffused, &
                                                            diag_height,sat_flow,resflow,discharge,removal
  integer                                 :: i,j
  

  if(previous .eq. current) then
    future = size(surf_liq,3) + 1 - current
  else
    future = previous
  end if

  call get_wts_lat(wts_lat)
  gwe = 0.0
  
!----------------------------------------------------------------------------
!  Begin by adding precip and evap, adjusting for groundwater evap
!----------------------------------------------------------------------------

  if (do_gwe .and. do_methane_table) then

     where (evap*delta_t .ge. surf_liq(:,:,current))
       meth_table(:,:,current) = meth_table(:,:,current) - (evap*delta_t - surf_liq(:,:,current))
       gwe = (evap*delta_t - surf_liq(:,:,current))/delta_t
       surf_liq(:,:,current) = 0.0
     elsewhere
       surf_liq(:,:,current) = surf_liq(:,:,current) - evap*delta_t
     endwhere

     surf_liq(:,:,current) = surf_liq(:,:,current) + tot_rain*delta_t/dt

  else
	
    if (do_leapfrog) then
       call perform_leap(surf_liq,tot_rain*delta_t/dt - evap*delta_t,previous,current,future)
    else
       surf_liq(:,:,current) = surf_liq(:,:,current) + tot_rain*delta_t/dt - evap*delta_t
    end if
	
  end if


!----------------------------------------------------------------------------
!  Diagnostically determine height based off of whether table is saturated
!----------------------------------------------------------------------------

  where (meth_table(:,:,current)/liq_dens/porosity .ge. dtd(:,js:je) - 1.e-10) !1e-10 to avoid numerical issues with small numbers
    diag_height = dtd(:,js:je) + surf_liq(:,:,current)/liq_dens
  elsewhere
    diag_height = meth_table(:,:,current)/liq_dens/porosity
  endwhere

!----------------------------------------------------------------------------
! Determine topo-driven runoff and infiltration
!----------------------------------------------------------------------------

  infiltration = 0.0; run_out = 0.0; dis_out = 0.0; recharge = 0.0; removal = 0.0; discharge = 0.0
  
  if (do_topo_runoff) then   
    call topo_runoff(runoff,removal,resflow,infiltration,discharge,surf_liq(:,:,current),run_res(:,:,current),&
          diag_height,tot_rain*delta_t/dt,evap*delta_t,lon,porosity,dt,do_methane_table)
    
  else if (do_global_infil) then
    where (surf_liq(:,:,current) .gt. infil_thresh)
      infiltration = min(surf_liq(:,:,current)-infil_thresh,infil_rate*dt*liq_dens)
    endwhere
    runoff = 0.0; resflow = 0.0; removal = 0.0; discharge = 0.0
  else 
    runoff = 0.0; resflow = 0.0; removal = 0.0; infiltration = 0.0; discharge = 0.0
  end if

   run_res(:,:,current) = run_res(:,:,current) + resflow*delta_t/dt
  !run_res and removal are always positive
 
!----------------------------------------------------------------------------
! Adjust table height with added infiltration and runoff
!----------------------------------------------------------------------------
  if (do_methane_table) then
   do i=1,size(lon,1)
    do j=1,size(lat,2)
      if (diag_height(i,j) .ge. dtd(is+i-1,js+j-1)) then
           !  where saturated, adjust height with discharge (no infil possible) 
          
         height_dt(i,j) = (discharge(i,j) - removal(i,j))/liq_dens  !+ runoff(i,j) if that way 
         
         if (diag_height(i,j) + height_dt(i,j) .ge. dtd(is+i-1,js+j-1)) then
            surf_liq_dt(i,j) = height_dt(i,j)*liq_dens
            meth_table_dt(i,j) = 0.0
         else  !this should never happen
            print *,'saturated table went below after runoff, height(i,j),height_dt(i,j),seep(i,j),surf(i,j),topo(i,j): ', diag_height(i,j),height_dt(i,j),removal(i,j)/liq_dens,surf_liq(i,j,current)/liq_dens,dtd(is+i-1,js+j-1),i,js+j-1
            height_dt(i,j) = height_dt(i,j)/porosity + (diag_height(i,j)-dtd(is+i-1,js+j-1))*(1.0/porosity - 1.0)
            surf_liq_dt(i,j) = -(diag_height(i,j)-dtd(is+i-1,js+j-1))*liq_dens
            meth_table_dt(i,j) = (height_dt(i,j) + diag_height(i,j) - dtd(is+i-1,js+j-1))*liq_dens*porosity
         end if
         
         infiltration(i,j) = 0.0  !to ensure no infil occurs when saturated
     
      else !unsaturated
           !  infiltration(i,j) = constant infiltration rate (hyd_cond)
 
         height_dt(i,j) = infiltration(i,j)/liq_dens/porosity
         if (diag_height(i,j) + height_dt(i,j) .ge. dtd(is+i-1,js+j-1)) then
           height_dt(i,j) = (dtd(is+i-1,js+j-1) - diag_height(i,j)) + (height_dt(i,j) + diag_height(i,j) - dtd(is+i-1,js+j-1))*porosity &
                + (surf_liq(i,j,current)-infiltration(i,j)-removal(i,j)+discharge(i,j))/liq_dens
           meth_table_dt(i,j) = (dtd(is+i-1,js+j-1)-diag_height(i,j))*porosity*liq_dens
           surf_liq_dt(i,j) = - meth_table_dt(i,j) - removal(i,j) + discharge(i,j)
         else
           height_dt(i,j) = height_dt(i,j)
           surf_liq_dt(i,j) = - removal(i,j) + discharge(i,j) - infiltration(i,j) 
           meth_table_dt(i,j) = infiltration(i,j)
         end if
       

      end if


    end do
   end do
  
   diag_height = diag_height + height_dt*delta_t/dt

  else
   height_dt = 0.0
   surf_liq_dt = - removal + discharge - infiltration
  end if

  dis_out = discharge/dt 
  run_out = removal/dt
  infiltration = infiltration/dt

!-----------------------------------------------------------------------
! Update subsurface
!-----------------------------------------------------------------------
  !Update surface to know what height changes to make when methane table reaches topo
   if (do_leapfrog) then
     call perform_leap(surf_liq,surf_liq_dt,previous,current,future)
   else
     surf_liq(:,:,current) = surf_liq(:,:,current) + surf_liq_dt*delta_t/dt
   end if


  !DO METHANE TABLE
  if (do_methane_table) then
   !Perform methane table flow and outcrop to surface
   diffused = 0.0; sat_flow = 0.0; subflow = 0.0
   call do_subsurface_flow(diffused,sat_flow,subflow,meth_table(:,:,current),diag_height,surf_liq(:,:,current), &
        lon,porosity,dt)
   meth_table_dt = meth_table_dt + subflow
   if (do_leapfrog) then
   ! Perform Leapfrog time-stepping with RAW filter
    call perform_leap(meth_table,meth_table_dt,previous,current,future)
   else    
     diag_height = diag_height + diffused*delta_t/dt
   end if

  else 
	
    if (do_leapfrog) then
      subflow = 0.0; sat_flow = 0.0; diffused = 0.0
      meth_table_dt = infiltration*dt
      call perform_leap(meth_table,meth_table_dt,previous,current,future)
    else
      subflow = 0.0; sat_flow = 0.0; diffused = 0.0
      meth_table_dt = infiltration*dt
      !!Add infiltration just to maintain mass conservation in printed inventories
      !meth_table(:,:,current) = meth_table(:,:,current) + infiltration*delta_t/dt 
    end if

  end if
 
!-----------------------------------------------------------------------
! Update surface liquid
!-----------------------------------------------------------------------
 
   surf_liq_dt = sat_flow

   if (do_leapfrog) then
     call perform_leap(meth_table,meth_table_dt,previous,current,future)
   else
     meth_table(:,:,current) = meth_table(:,:,current) + meth_table_dt*delta_t/dt
   end if

   if (do_leapfrog) then
     call perform_leap(surf_liq,surf_liq_dt,previous,current,future)
   else
     surf_liq(:,:,current) = surf_liq(:,:,current) + surf_liq_dt*delta_t/dt
   end if

  recharge = sat_flow/dt
  subflow = subflow/dt

  height(:,:,current) = diag_height

   do i=1,size(lon,1)
    do j=1,size(lat,2)
       if (surf_liq(i,j,current) .lt. 0.0) print *,'negative surf: ',surf_liq(i,j,current), i,js+j-1
    end do
   end do
 

  if (.not.do_leapfrog) then
   surf_liq(:,:,previous) = surf_liq(:,:,current)
   run_res(:,:,previous) = run_res(:,:,current)
   meth_table(:,:,previous) = meth_table(:,:,current)
   height(:,:,previous) = height(:,:,current)
  end if
 
  if (first_hydro_call) first_hydro_call = .false.


end subroutine tam_hydrology_driver


subroutine perform_leap(param,change,previous,current,future)
  real, intent(in) , dimension(:,:)       :: change
  real, intent(inout), dimension(:,:,:)   :: param
  integer, intent(in)   :: previous,current,future
  
  real, dimension(size(param,1),size(param,2)) :: filt
  integer :: i,j

   filt = param(:,:,previous) - 2.0 * param(:,:,current)

   if (previous .eq. current) then
     param(:,:,future) = param(:,:,previous) + change
     param(:,:,current) = param(:,:,current) + robert_surf_liq * RAW_surf_liq * &
        (param(:,:,previous) - 2.0*param(:,:,current) + param(:,:,future))
   else
     param(:,:,current) = param(:,:,current) + robert_surf_liq * RAW_surf_liq * &
        (param(:,:,previous) - 2.0*param(:,:,current))
     param(:,:,future) = param(:,:,previous) + change
     param(:,:,current) = param(:,:,current) + robert_surf_liq * RAW_surf_liq * &
        param(:,:,future)
   end if
   param(:,:,future) = param(:,:,future) + robert_surf_liq * &
        (filt + param(:,:,future)) * (RAW_surf_liq - 1.0)

end subroutine perform_leap


subroutine topo_runoff(runoff, removal, res,infil,dis,liquid,rr,height,rain, evap,lon, porosity,dt, do_methane_table)

  real, intent(in),    dimension(:,:)   :: lon, liquid, rr, height, rain, evap  !all in kg/m2 except height (m)!
  real, intent(inout), dimension(:,:)   :: runoff, res, removal, infil, dis
  real, intent(in)                      :: dt,porosity
  logical, intent(in)                   :: do_methane_table
 
  integer :: i,j,ix,jx,ixs,ixn,ilws,ilwn,iles,ilen,ilw,ile,jln,jls !indices 
  integer :: px, py !processor indices
  integer, dimension(nx+2) :: ii   !polar circle function

  real, dimension(is:ie,js:je)  :: seep,sd_map  !seep is removal, sd_map gives steepest descent index at each point
  real, dimension(is:ie,jsd:jed) :: dld,drd,final_runoff,drd_north,drd_south,dresd,final_res,dresd_n,dresd_s  !data liq/runoff/run_res domains
  real, dimension(ny) :: wts_lat   !latitude weights

  real, dimension(is:ie,jsd:jed) :: cp !cell potential
  real :: d1,d2,d3,d4,d5,d6,d7,d8,d9,c1,c2,c3,c4 !slopes and distances 
  real :: curr_liq,ue,ur,mt,rrc  !cell liq,updated evap,rain,table

  integer    :: steep_des  !grid number of steepest descent
  real, dimension(1:9)   :: descents  !descent array


!       -------------------------
!       |       |       |       |
!       |  1    |   2   |   3   |
!       |-------|--c1---|-------|
!       |      c4       |       |
!       |  4    |  i,j c2   5   |   d9 is always 0
!       |-------|--c3---|-------|
!       |       |       |       |
!       |  6    |   7   |   8   |
!       |-------|-------|-------|


!-----------------------------------------------------------------------
! Set up polar circle function
!-----------------------------------------------------------------------
  
  do i = 2,nx+1
    ii(i) = i-1+nx/2
    if (ii(i) .gt. nx) ii(i) = ii(i) - nx
  end do
  ii(1) = ii(nx+1)
  ii(nx+2) = ii(2)

!-----------------------------------------------------------------------
! Communicate with neighboring cells to get current liquid
!-----------------------------------------------------------------------  
      
  px = size(runoff,1)
  py = size(runoff,2)
  drd = 0.0; dld = 0.0; final_runoff = 0.0
  drd_north = 0.0; drd_south = 0.0
  dresd = 0.0; dresd_n = 0.0; dresd_s = 0.0; final_res = 0.0; 
  seep = 0.0; infil = 0.0; dis = 0.0
 
  !Infiltrate first!
  if (infil_first) then
   do j = 1,py
     do i = 1,px
  
       if (liquid(i,j) .gt. infil_thresh) then 
             if ((do_methane_table .and. (height(i,j) .lt. dtd(i+is-1,j+js-1))) .or. (do_global_infil .and. .not.do_methane_table)) then 
                infil(i,j) = min(liquid(i,j) - infil_thresh,infil_rate*dt*liq_dens)
             end if 
       end if 
     end do
   end do
  end if


  do j = js,je
    do i = is,ie

      dld(i,j) = ( liquid(i-is+1,j-js+1) - infil(i-is+1,j-js+1) )/liq_dens  !Data Liquid Domain (in meters)

    end do
  end do

  call get_wts_lat(wts_lat)
  call mpp_update_domains( dld, topo_domain)
 
 !Cell potential = Topography + Surface liquid
  cp = dld+dtd !in meters


!-----------------------------------------------------------------------
! ONE MAJOR LOOP through grid points to assign steepest descent and move flow!!
!-----------------------------------------------------------------------
  
  do j = 1,py
    do i = 1,px

     !Move from processor indices to grid indices
        jx = js+j-1 
        ix = is+i-1

!-----------------------------------------------------------------------
!  Boundary conditions
!-----------------------------------------------------------------------
     !E-W boundary conditions 
       if(ix == 1) then
          ilw = px
        else
          ilw = ix-1
        end if

        if (ix == px) then
          ile = 1
        else
          ile = ix+1
        end if 
     !these will change if on polar circle
        ixs = ix  
        ixn = ix  
        ilws = ilw 
        iles = ile 
        ilwn = ilw 
        ilen = ile 
     !N-S boundary conditions 
        if(jx == 1) then    !south pole
          jls = 1
          jln = jx+1
          ixs = ii(ix+1) 
          ilws = ii(ix+2) 
          iles = ii(ix)
        else if (jx == ny) then   !north pole
          jls = jx-1
          jln = ny
          ixn = ii(ix+1)  
          ilwn = ii(ix+2) 
          ilen = ii(ix) 
        else
          jls = jx-1
          jln = jx+1
        end if
 
!-----------------------------------------------------------------------
!  Calculate slopes and descents
!-----------------------------------------------------------------------       

        c1 = (lat_halo(ix,jln) - lat_halo(ix,jx))*Radius
        if (ile == 1) then
           c2 = ((2*pi - lon(ix,j))+lon(ile,j))*Radius*cos(lat_halo(ix,jx))
        else 
           c2 = (lon(ile,j) - lon(ix,j))*Radius*cos(lat_halo(ix,jx))
        end if
        c3 = (lat_halo(ix,jx) - lat_halo(ix,jls))*Radius
        if (ilw == px) then
           c4 = ((2*pi - lon(ilw,j))+lon(ix,j))*Radius*cos(lat_halo(ix,jx))
        else
           c4 = abs(lon(ilw,j) - lon(ix,j))*Radius*cos(lat_halo(ix,jx))
        end if
        if (c1 .eq. 0.0) c1 = 2.0*(pi/2.0-lat_halo(ix,jln))*Radius
        if (c3 .eq. 0.0) c3 = 2.0*(pi/2.0-lat_halo(ix,jls))*Radius

        d1 = (cp(ilwn,jln) - cp(ix,jx))/sqrt(c1**2 + c4**2)
        d2 = (cp(ixn,jln) - cp(ix,jx))/c1
        d3 = (cp(ilen,jln) - cp(ix,jx))/sqrt(c1**2 + c2**2)
        d4 = (cp(ilw,jx) - cp(ix,jx))/c4
        d5 = (cp(ile,jx) - cp(ix,jx))/c2
        d6 = (cp(ilws,jls) - cp(ix,jx))/sqrt(c3**2 + c4**2)
        d7 = (cp(ixs,jls) - cp(ix,jx))/c3
        d8 = (cp(iles,jls) - cp(ix,jx))/sqrt(c3**2 + c2**2)
        d9 = 0.0
   
      
        descents = (/ d9,d1,d2,d3,d4,d5,d6,d7,d8 /)
        steep_des = minloc(descents, DIM=1)-1  !if steep_des = 0, then (ix,jx) is a depression
        if (abs(descents(steep_des+1)) .lt. 1e-8) steep_des = 0   !arbitrary slope minmium
        sd_map(ix,jx) = steep_des 
!-----------------------------------------------------------------------
!  Parameterize runoff speed and move the flow
!-----------------------------------------------------------------------       

        curr_liq = dld(ix,jx)*liq_dens

        ue = evap(i,j)      !updated evap 
        ur = rain(i,j)      !updated rain
        mt = height(i,j)    !methane table 
        rrc = rr(i,j)       !runoff reservoir current 

          !PC flow
          !Each call adds the flow to the steepest descent grid and takes away
          !the flow in (ix,jx) in the array 'dresd'. If the steepest descent grid
          !is on another processor, the flow there is recorded at (ix,jx) in the
          !'dresd_s'/'dresd_n' array. For example, if flow moves from
          !(10,17) to (10,16) across a processor (where je = 16, jed = 17), then an amount is subtracted
          !from (10,17) in 'dresd' and added to (10,17) in 'drd_south'
          if (steep_des .eq. 1) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,sqrt(c1**2+c4**2),curr_liq,ix,jx,ilwn,jln,mt,ur,ue,d1,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 2) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,c1,curr_liq,ix,jx,ixn,jln,mt,ur,ue,d2,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 3) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,sqrt(c1**2+c2**2),curr_liq,ix,jx,ilen,jln,mt,ur,ue,d3,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 4) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,c4,curr_liq,ix,jx,ilw,jx,mt,ur,ue,d4,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 5) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,c2,curr_liq,ix,jx,ile,jx,mt,ur,ue,d5,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 6) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,sqrt(c3**2+c4**2),curr_liq,ix,jx,ilws,jls,mt,ur,ue,d6,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 7) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,c3,curr_liq,ix,jx,ixs,jls,mt,ur,ue,d7,dt,dresd,dresd_n,dresd_s,seep,rrc)
          if (steep_des .eq. 8) call gosurfflow(drd,drd_north,drd_south,cp,wts_lat,sqrt(c3**2 + c2**2),curr_liq,ix,jx,iles,jls,mt,ur,ue,d8,dt,dresd,dresd_n,dresd_s,seep,rrc)
    end do
  end do 
 
!-----------------------------------------------------------------------
! Add all accumulated runoff
!-----------------------------------------------------------------------       

  final_runoff = final_runoff + drd                 

  !Fill in halo points, can now access flow moved in from other processors
  call mpp_update_domains( drd_south, topo_domain)
  call mpp_update_domains( drd_north, topo_domain)

  !On southern edge of processor, add flow that moved north
  final_runoff(:,js) = final_runoff(:,js) + drd_north(:,jsd) 
  !Likewise on northern edge. From prior example, at (10,16), we can access
  !drd_south(10,17), which tells us how much ran into (10,16)
  final_runoff(:,je) = final_runoff(:,je) + drd_south(:,jed) 

!----------runoff reservoir movement-----
 
  final_res = final_res + dresd
  call mpp_update_domains(dresd_s,topo_domain)
  call mpp_update_domains(dresd_n,topo_domain)
  final_res(:,js) = final_res(:,js) + dresd_n(:,jsd)
  final_res(:,je) = final_res(:,je) + dresd_s(:,jed)




!-----------------------------------------------------------------------
!  Loop through again to get dis/infil based on runoff amounts
!-----------------------------------------------------------------------
  do j = 1,py
    do i = 1,px
       ix = is+i-1 
       jx = js+j-1
         
         curr_liq = dld(ix,jx)*liq_dens          
         rrc = rr(i,j)
          
          if (do_methane_table .and. (height(i,j) .ge. dtd(ix,jx))) then
             infil(i,j) = 0.0     !ensure no infiltration
             if (sd_map(ix,jx) .eq. 0) then   !Local depression of saturated liquid
               dis(i,j) = final_res(ix,jx) + rrc    !add the amount added this loop (dresd(ix,jx) should be positive at local depression) AND amount in there already
               final_res(ix,jx) = -rrc               !subtract current run res amount from res
             else !Shoreline or near-shoreline interior of sea
               dis(i,j) = max(final_res(ix,jx),0.0)  !anything added from surrouding cells will move to surface 
               final_res(ix,jx) = final_res(ix,jx) - dis(i,j)    
               !whatever is in this cell's rr will be moved to a lower cell
             end if 
          end if
         
          !Local depression or perched liquid 
          if ((sd_map(ix,jx) .eq. 0 .and. .not.do_methane_table) .or. ((sd_map(ix,jx) .eq. 0) .and. do_methane_table .and. (height(i,j) .lt. dtd(ix,jx)))) then 
             dis(i,j) = final_res(ix,jx) + rrc  !add the amount added this loop AND the amount in there already
             final_res(ix,jx) =  -rrc  
          end if

          if (.not.do_global_infil) infil(i,j) = 0.0
          
    end do
  end do 

!-----------------------------------------------------------------------
! Assign to runoff diag
!-----------------------------------------------------------------------       

  do j = 1,py
    do i = 1,px
        jx = js+j-1
        ix = is+i-1
        runoff(i,j) = final_runoff(ix,jx)
        res(i,j) = final_res(ix,jx)
        removal(i,j) = seep(ix,jx)

        if (.not.infil_first) then
          if (liquid(i,j)+dis(i,j) .gt. infil_thresh) then
             if ((do_methane_table .and. (height(i,j) .lt. dtd(i+is-1,j+js-1))) .or. (do_global_infil .and. .not.do_methane_table)) then 
                infil(i,j) = min(liquid(i,j) + dis(i,j) - infil_thresh,infil_rate*dt*liq_dens)
             end if 
          end if 

         if (.not.do_global_infil) infil(i,j) = 0.0
        end if

    end do
  end do




end subroutine topo_runoff


subroutine gosurfflow(drd,drd_north,drd_south,cp,wts,length,curr_liq,ix,jx,rx,ry,height,rain,ue,slope,dt, &
                      dresd,dresd_n,dresd_s,seep,rrc)

  real, intent(in)                              :: length,curr_liq,height,rain,ue,slope,dt,rrc
  real, intent(in), dimension(is:ie,jsd:jed)    :: cp
  real, intent(in), dimension(ny)               :: wts
  real, intent(inout), dimension(is:ie,jsd:jed) :: drd,drd_north,drd_south,dresd,dresd_n,dresd_s
  real, intent(inout), dimension(is:ie,js:je)           :: seep
  
  integer, intent(in)                   :: ix,jx,rx,ry

  real  :: tc,cell_run,max_run,goflow,runflow

  !Maximum runoff from one cell to another
  max_run = 0.5*(cp(ix,jx) - cp(rx,ry))*liq_dens*wts(ry)/wts(jx)  !everything in kg/m2
  if (max_run .lt. 0.0) then 
      max_run = 0.0
  end if

  if (linear_catch) then  !see Watt and Chow (1985) and Horvath et al. (2016)
      tc = 3600.*3.23*0.000326*(basin_length/abs(slope)**0.5)**0.79  !3.23 is scaling due to lower g, rho, and viscosity of liquid methane
      !tc = 3600.*0.41344*((length)/abs(slope)**0.5)**0.79  !this was old way, probably wrong but nicely slow
      cell_run = 2.*dt*curr_liq/tc
 
   !Print out rates of interest in spectral.out
    if (first_hydro_call) print *, 'tc,linear_catch,length,L/S,rain,evap,jx: ',tc,cell_run,length, length/sqrt(abs(slope)),rain,ue,jx
    if (rain*1000*86400/dt/liq_dens .gt. 100.0) print *, 'intense storm > 100 mm/day, rain (mm/day), runoff,tc : ',rain*1000*86400/dt/liq_dens,cell_run,tc,jx

  else if (input_speed) then
      cell_run = run_speed*dt*liq_dens*abs(slope)
  else
      cell_run = 0.0
  end if
 
 !If liquid is still above threshold after infiltration, then go for the flow!! 
 if ((curr_liq .gt. infil_thresh) .and. (abs(slope) .gt. 1.e-8)) then   !arbitrary liquid and slope minima
  if ((cell_run .lt. curr_liq) .and. (cell_run .lt. max_run)) then
      goflow = cell_run
  else if ((cell_run .lt. curr_liq) .and. (cell_run .ge. max_run)) then  
      goflow = max_run
  else if ((cell_run .gt. curr_liq) .and. (curr_liq .lt. max_run)) then
      goflow = curr_liq*.95  !.95 to avoid numerical issues with going too close to 0
  else if ((cell_run .gt. curr_liq) .and. (curr_liq .ge. max_run)) then
      goflow = max_run
  end if
 else
   goflow = 0.0
 end if

 runflow = goflow + rrc 

  if (ry .eq. jsd) then
     drd_south(rx,jx) = drd_south(rx,jx) + goflow*wts(jx)/wts(ry)   !account for latitude
     dresd_s(rx,jx) = dresd_s(rx,jx) + runflow*wts(jx)/wts(ry)
  else if (ry .eq. jed) then
     drd_north(rx,jx) = drd_north(rx,jx) + goflow*wts(jx)/wts(ry)
     dresd_n(rx,jx) = dresd_n(rx,jx) + runflow*wts(jx)/wts(ry)
  else
     drd(rx,ry) = drd(rx,ry) + goflow*wts(jx)/wts(ry)
     dresd(rx,ry) = dresd(rx,ry) + runflow*wts(jx)/wts(ry)
  end if 
  drd(ix,jx) = drd(ix,jx) - goflow
  dresd(ix,jx) = dresd(ix,jx) - rrc
  seep(ix,jx) = goflow   !To be subtracted from surface reservoir...


end subroutine gosurfflow 

SUBROUTINE do_subsurface_flow(diff,sat_flow,subflow,table,height,liquid,lon,porosity,dt)
  
	real, intent(in),    dimension(:,:)   :: lon,table,height,liquid  !liquid is in kg/m2!
	real, intent(inout), dimension(:,:)   :: diff,sat_flow,subflow
	real, intent(in)                      :: dt,porosity
 
  	integer :: i,j,k,ix,jx,ixs,ixn,ilw,ilws,ilwn,ile,iles,ilen,jln,jls
  	integer :: px,py
  	integer, dimension(nx+2) :: ii 

  	real, dimension(is:ie,jsd:jed) 	:: dld,drd,drd_north,drd_south,final_flow !drd's and final_flow are volumes!
  	real, dimension(ny) 			:: wts_lat
  	real, allocatable, dimension(:) :: blon,blat

  	real, dimension(is:ie,jsd:jed)	:: cp 			!cell potential
  	real :: d1,d2,d3,d4,d5,d6,d7,d8,d9,c1,c2,c3,c4	!slopes and distances 
  	real :: cell_run_n,cell_run_s,cell_run_e,cell_run_w,x,total_cell_run,updat_total_cell_run  ! cell flows, L, total volume removed from cell
        real :: cell_area, cell_dh

  	integer    				:: kmax, num_slopes 
  	real, dimension(1:9)  	:: descents
  	real, dimension(1:4)  	:: descents_grid,dfrac,dfrac_dir,dfrac_temp

  	!"Real" diffusion
  	real, dimension(is:ie,js:je) 					:: dt_table,table_diffusion
  	real, dimension(size(liquid,1),size(liquid,2))	:: diag_height
  	
  	px = size(subflow,1)	!number of longitudes on processor
  	py = size(subflow,2)	!number of latitudes on processor	
	
	allocate(blon(is:ie+1), blat(js:je+1))
	call get_grid_boundaries(blon, blat)
    drd = 0.0; drd_north = 0.0; drd_south = 0.0; final_flow = 0.0; dld = 0.0
  	diff = 0.0				!change in methane table height
  	diag_height = height	!overall methane table height 
    
  	if (do_diff) then	!If doing "real" diffusion (with spherical harmonics)
    	
    	dt_table = 0.0
    	table_diffusion = 0.0 
  		
    	call diffuse_table(dt_table,diag_height,dt,hyd_cond,table_diffusion,porosity)
    	
    	do i=is,ie
    		do j=js,je
        		if (diag_height(i-is+1,j-js+1) .ge. dtd(i,j)) then
          			if ((diag_height(i-is+1,j-js+1) + dt_table(i,j)) .ge. dtd(i,j)) then
            			diff(i-is+1,j-js+1) = dt_table(i,j)*porosity                         
            			sat_flow(i-is+1,j-js+1) = dt_table(i,j)*porosity*liq_dens
            			subflow(i-is+1,j-js+1) = 0.0
          			else
            			diff(i-is+1,j-js+1) = (diag_height(i-is+1,j-js+1) + dt_table(i,j) + (diag_height(i-is+1,j-js+1)-dtd(i,j))*(1.0-porosity)) - diag_height(i-is+1,j-js+1)
            			sat_flow(i-is+1,j-js+1) = -(diag_height(i-is+1,j-js+1) - dtd(i,j))*porosity*liq_dens
            			subflow(i-is+1,j-js+1) = (dt_table(i,j) + diag_height(i-is+1,j-js+1) - dtd(i,j))*porosity*liq_dens
          			end if
        		else
          			if ((diag_height(i-is+1,j-js+1) + dt_table(i,j)) .ge. dtd(i,j)) then
            			diff(i-is+1,j-js+1) = (dtd(i,j) + (diag_height(i-is+1,j-js+1)-dtd(i,j)+dt_table(i,j))*porosity) - diag_height(i-is+1,j-js+1) + liquid(i-is+1,j-js+1)/liq_dens!need to add surface liquid
            			sat_flow(i-is+1,j-js+1) = (diag_height(i-is+1,j-js+1)-dtd(i,j)+dt_table(i,j))*porosity*liq_dens
            			subflow(i-is+1,j-js+1) = (dtd(i,j) - diag_height(i-is+1,j-js+1))*porosity*liq_dens
            			sat_flow(i-is+1,j-js+1) = 0.0
            			subflow(i-is+1,j-js+1) = dt_table(i,j)*porosity*liq_dens
          			end if
       		 	end if
     		end do
    	end do

	elseif (do_darcy) then	!Continue with discretized Darcy flow 
  		
  		do i = 2,nx+1
    		ii(i) = i-1+nx/2
    		if (ii(i) .gt. nx) ii(i) = ii(i) - nx
  		end do
  		ii(1) = ii(nx+1)
  		ii(nx+2) = ii(2)
    	
  		do i=is,ie
    		do j=js,je
				dld(i,j) = diag_height(i-is+1,j-js+1)
    		end do
  		end do

  		call get_wts_lat(wts_lat)
                print *,'dld dim 2 x js: ',js*size(dld,2)
		!print *,'topo_domain dim 2 x js: ',js*size(topo_domain,2)
                call mpp_update_domains(dld, topo_domain) !mmm segfaults here when run on multiple cores
		
		cp = dld	!cell potentials, in meters
		
		do j = 1,py
			do i = 1,px
        		jx = js+j-1
        		ix = is+i-1
        		if(ix == 1) then
          			ilw = px
        		else
          			ilw = ix-1
        		end if
        		if (ix == px) then
          			ile = 1
        		else
          			ile = ix+1
        		end if 
 
     			!these will change if on polar circle
        		ixs = ix  
        		ixn = ix  
        		ilws = ilw 
        		iles = ile 
        		ilwn = ilw 
        		ilen = ile 
     			!N-S boundary conditions 
        		if(jx == 1) then    !south pole
          			jls = 1
          			jln = jx+1
          			ixs = ii(ix+1) 
          			ilws = ii(ix+2) 
          			iles = ii(ix)
        		else if (jx == ny) then   !north pole
          			jls = jx-1
          			jln = ny
          			ixn = ii(ix+1)  
          			ilwn = ii(ix+2) 
          			ilen = ii(ix) 
        		else
          			jls = jx-1
          			jln = jx+1
        		end if
       
        		c1 = (lat_halo(ix,jln) - lat_halo(ix,jx))*Radius
        		if (ile == 1) then
           			c2 = ((2.*pi - lon(ix,j))+lon(ile,j))*Radius*cos(lat_halo(ix,jx))
        		else 
           			c2 = (lon(ile,j) - lon(ix,j))*Radius*cos(lat_halo(ix,jx))
        		end if
       			c3 = (lat_halo(ix,jx) - lat_halo(ix,jls))*Radius
        		if (ilw == px) then
           			c4 = ((2.*pi - lon(ilw,j))+lon(ix,j))*Radius*cos(lat_halo(ix,jx))
        		else
           			c4 = abs(lon(ilw,j) - lon(ix,j))*Radius*cos(lat_halo(ix,jx))
        		end if
        		if (c1 .eq. 0.0) c1 = 2.0*(pi/2.0-lat_halo(ix,jln))*Radius
        		if (c3 .eq. 0.0) c3 = 2.0*(pi/2.0-lat_halo(ix,jls))*Radius


 				!PC descents
        		d1 = (cp(ilwn,jln) - cp(ix,jx))/sqrt(c1**2 + c4**2)
        		d2 = (cp(ixn,jln) - cp(ix,jx))/c1
        		d3 = (cp(ilen,jln) - cp(ix,jx))/sqrt(c1**2 + c2**2)
        		d4 = (cp(ilw,jx) - cp(ix,jx))/c4
        		d5 = (cp(ile,jx) - cp(ix,jx))/c2
        		d6 = (cp(ilws,jls) - cp(ix,jx))/sqrt(c3**2 + c4**2)
        		d7 = (cp(ixs,jls) - cp(ix,jx))/c3
        		d8 = (cp(iles,jls) - cp(ix,jx))/sqrt(c3**2 + c2**2)
        		d9 = 0.0
       
        		descents = (/ d9,d1,d2,d3,d4,d5,d6,d7,d8 /)
  				descents_grid = (/d2,d4,d5,d7 /)

          		dfrac = 0.0
          		num_slopes = 0
          		do k = 1,4		!collect cells of lower potential
            		if (descents_grid(k) < 0) then 
              			dfrac_temp(k) = abs(descents_grid(k))
              			num_slopes = num_slopes + 1
            		else 
              			dfrac_temp(k) = 0.
            		end if
          		end do

         		!Sort dfrac and dfrac_dir so the cell to steepest descent is called first
          		do k = 1,4
            		dfrac(k)		= maxval(dfrac_temp)
            		kmax 			= maxloc(dfrac_temp,DIM=1)
            		dfrac_dir(k)	= kmax
            		dfrac_temp(kmax)= 0.
          		end do
!-----------------------------------------------------------------------------------------
!8/21/2018 JML
!This code discretizes Darcy flow
!The equation of interest is based on flow of liquid volume that is proportional to the 
!	head gradient: porosity * dh * area = K * B * x * dh/dl * dt
!	where K is conductivity, B is depth across which flow is operating, and x is the 
!	horizontal length of contact between cells.
!In the following, cell_run = porosity * dh * area = volume of liquid that moves out
!-----------------------------------------------------------------------------------------
                    cell_run_n = 0.0; cell_run_s = 0.0
                    cell_run_w = 0.0; cell_run_e = 0.0
 			do while (num_slopes > 0) 
           			if (dfrac_dir(num_slopes) .eq. 1) then		!northward flow
           				x = 2.*pi*radius * cos(blat(jx+1)) / nx
           				cell_run_n = hyd_cond * cp(ix,jx) * x * (cp(ix,jx)-cp(ixn,jln))/c1 * dt
                                        call checkmaxflow(cell_run_n,cp,wts_lat,ix,jx,ixn,jln,porosity)
           			end if
                    if (dfrac_dir(num_slopes) .eq. 2) then      !westward flow
           		x = pi*radius / ny
                        cell_run_w = hyd_cond * cp(ix,jx) * x * (cp(ix,jx)-cp(ilw,jx))/c4 * dt
                        call checkmaxflow(cell_run_w,cp,wts_lat,ix,jx,ilw,jx,porosity)
                    end if 
                    if (dfrac_dir(num_slopes) .eq. 3) then		!eastward flow
                    	x = pi*radius / ny
                        cell_run_e = hyd_cond * cp(ix,jx) * x * (cp(ix,jx)-cp(ile,jx))/c2 * dt
                        call checkmaxflow(cell_run_e,cp,wts_lat,ix,jx,ile,jx,porosity)
                    end if 
                	if (dfrac_dir(num_slopes) .eq. 4) then		!southward flow
                    	x = 2.*pi*radius * cos(blat(jx)) / nx
                        cell_run_s = hyd_cond * cp(ix,jx) * x * (cp(ix,jx)-cp(ixs,jls))/c3 * dt
                        call checkmaxflow(cell_run_s,cp,wts_lat,ix,jx,ixs,jls,porosity)
                    end if
           			num_slopes = num_slopes - 1
           		end do
          		
                        total_cell_run = cell_run_n+cell_run_w+cell_run_e+cell_run_s
                        updat_total_cell_run = total_cell_run
                        call removesubflow(drd,cp,wts_lat,total_cell_run,ix,jx,porosity,updat_total_cell_run)  !removes total volume from ix,jx
           		if (cell_run_n .ne. 0.0) cell_run_n = cell_run_n*updat_total_cell_run/total_cell_run
           		if (cell_run_w .ne. 0.0) cell_run_w = cell_run_w*updat_total_cell_run/total_cell_run
           		if (cell_run_e .ne. 0.0) cell_run_e = cell_run_e*updat_total_cell_run/total_cell_run
           		if (cell_run_s .ne. 0.0) cell_run_s = cell_run_s*updat_total_cell_run/total_cell_run
                        call gosubflow(drd,drd_north,drd_south,cp,wts_lat,cell_run_n,ix,jx,ixn,jln,porosity)
           		call gosubflow(drd,drd_north,drd_south,cp,wts_lat,cell_run_w,ix,jx,ilw,jx,porosity)
           		call gosubflow(drd,drd_north,drd_south,cp,wts_lat,cell_run_e,ix,jx,ile,jx,porosity)
           		call gosubflow(drd,drd_north,drd_south,cp,wts_lat,cell_run_s,ix,jx,ixs,jls,porosity)
                       
			end do
		end do
                !combine all volume removals and additions for each cell into 'final_flow' array
  		final_flow = final_flow + drd 
  		call mpp_update_domains( drd_south, topo_domain)
  		call mpp_update_domains( drd_north, topo_domain) !mmm segfaults here when run on multiple cores
		final_flow(:,js) = final_flow(:,js) + drd_north(:,jsd) 
		final_flow(:,je) = final_flow(:,je) + drd_south(:,jed)
		
		!Diagnostics
                cell_dh = 0.0; cell_area = 0.0
		do j = 1,py
			do i = 1,px
			jx = js+j-1
        		ix = is+i-1
        		
                        cell_area = 2.*pi*radius*radius * wts_lat(jx) / nx	!horizontal area of cell (m**2)

                        !calculate cell change in m from changes in liquid volume, adjusting for porosity
                        if (cp(ix,jx) .gt. dtd(ix,jx)) then					!Table above surface
                            if (final_flow(ix,jx) .lt. 0.0) then                                    !decreasing table
                                if (-final_flow(ix,jx) .lt. (cp(ix,jx)-dtd(ix,jx))*cell_area) then
                                        cell_dh = final_flow(ix,jx)/cell_area
                                else
				        cell_dh = -((-final_flow(ix,jx) - (cp(ix,jx)-dtd(ix,jx))*cell_area)/cell_area/porosity + (cp(ix,jx)-dtd(ix,jx)))
                                end if
                            else                                                                    !increasing table
                                cell_dh = final_flow(ix,jx)/cell_area
                            end if
                        else                                                                    !Table below surface
                            if (final_flow(ix,jx) .lt. 0.0) then                                    !decreasing table
                                cell_dh = final_flow(ix,jx)/cell_area/porosity
                            else                                                                    !increasing table
                                if (final_flow(ix,jx) .gt. (dtd(ix,jx)-cp(ix,jx))*cell_area*porosity) then
                                        cell_dh = (dtd(ix,jx)-cp(ix,jx)) + (final_flow(ix,jx) - (dtd(ix,jx)-cp(ix,jx))*cell_area*porosity)/cell_area
                                else
                                        cell_dh = final_flow(ix,jx)/cell_area/porosity
        	                end if
                            end if
                        end if	

                        diff(i,j) = cell_dh	!change of overall methane table height (m)


        		if (diag_height(i,j) .ge. dtd(ix,jx)) then
        			if ((diag_height(i,j) + cell_dh) .ge. dtd(ix,jx)) then
        				sat_flow(i,j) = cell_dh*liq_dens
        				subflow(i,j)  = 0.0
        			else
        				sat_flow(i,j) = -(diag_height(i,j) - dtd(ix,jx))*liq_dens
        				subflow(i,j)  = (cell_dh + diag_height(i,j) - dtd(ix,jx))*porosity*liq_dens
        			end if
        		else
        			if ((diag_height(i,j) + cell_dh) .ge. dtd(ix,jx)) then
        				sat_flow(i,j) = (diag_height(i,j) - dtd(ix,jx) + cell_dh)*liq_dens
        				subflow(i,j)  = (dtd(ix,jx) - diag_height(i,j))*porosity*liq_dens
        			else
        				sat_flow(i,j) = 0.0
    					subflow(i,j)  = cell_dh*porosity*liq_dens
    				end if
    			end if
			end do
		end do
	else 
   		subflow = 0.0; sat_flow = 0.0; diff = 0.0
  	end if !end of if(do_diff) statement at very top
 
	deallocate(blon,blat)

END SUBROUTINE do_subsurface_flow

SUBROUTINE checkmaxflow(cell_run,table,wts_lat,ix,jx,rx,ry,porosity)

	real, intent(in)                            	:: porosity
	real, intent(inout)                            	:: cell_run
	real, intent(in), dimension(is:ie,jsd:jed)   	:: table
	real, intent(in), dimension(ny)					:: wts_lat
	integer, intent(in)                   			:: ix,jx,rx,ry
  
	real	:: cell_area,max_dh,max_vol

	cell_area = 2.*pi*radius*radius * wts_lat(jx) / nx	!horizontal area of cell (m**2)


        max_dh = .5*(table(ix,jx) - table(rx,ry))
	if (table(ix,jx) .gt. dtd(ix,jx)) then					!Table above surface
                if (table(rx,ry) .gt. dtd(rx,ry)) then
                        max_vol = max_dh*cell_area
                else
                        if ((table(ix,jx)-max_dh) .gt. dtd(ix,jx)) then
                                max_vol = max_dh*cell_area
                        else
                                max_vol = (table(ix,jx)-dtd(ix,jx))*cell_area + (max_dh - (table(ix,jx)-dtd(ix,jx)))*cell_area*porosity
                        end if
                end if
        else                                                                    !Table below surface
                max_vol = max_dh*cell_area*porosity
        end if

        if (cell_run .gt. max_vol) cell_run = max_vol
		
END SUBROUTINE checkmaxflow

SUBROUTINE removesubflow(drd,table,wts_lat,cell_run_temp,ix,jx,porosity,updat_run)

	real, intent(in)                            	:: porosity,cell_run_temp
	real, intent(inout)                            	:: updat_run
	real, intent(in), dimension(is:ie,jsd:jed)   	:: table
	real, intent(in), dimension(ny)					:: wts_lat
	real, intent(inout), dimension(is:ie,jsd:jed)	:: drd
	integer, intent(in)                   			:: ix,jx
  
	real	:: cell_run,cell_area,avail_vol

        cell_run = cell_run_temp
	if (cell_run .lt. 0.0) cell_run = 0.0

	cell_area = 2.*pi*radius*radius * wts_lat(jx) / nx	!horizontal area of cell (m**2)


	if (table(ix,jx) .gt. dtd(ix,jx)) then					!Table above surface
                avail_vol = dtd(ix,jx)*cell_area*porosity + (table(ix,jx)-dtd(ix,jx))*cell_area
        else
                avail_vol = table(ix,jx)*cell_area*porosity
        end if

        if (cell_run .gt. avail_vol) cell_run = avail_vol
        updat_run = cell_run
        drd(ix,jx) = drd(ix,jx) - cell_run
		
END SUBROUTINE removesubflow

SUBROUTINE gosubflow(drd,drd_north,drd_south,table,wts_lat,cell_run_temp,ix,jx,rx,ry,porosity)

	real, intent(in)                            	:: porosity,cell_run_temp
	real, intent(in), dimension(is:ie,jsd:jed)   	:: table
	real, intent(in), dimension(ny)					:: wts_lat
	real, intent(inout), dimension(is:ie,jsd:jed)	:: drd,drd_north,drd_south
	integer, intent(in)                   			:: ix,jx,rx,ry
  
	real	:: cell_run

        cell_run = cell_run_temp
	if (cell_run .lt. 0.0) cell_run = 0.0

       if (cell_run .ne. 0.0) then
	if (ry .eq. jsd) then
		drd_south(rx,jx) = drd_south(rx,jx) + cell_run
	else if (ry .eq. jed) then
		drd_north(rx,jx) = drd_north(rx,jx) + cell_run
	else
	        drd(rx,ry) = drd(rx,ry) + cell_run
	end if
       end if
  
END SUBROUTINE gosubflow 


subroutine diffuse_table(dt_table,table,dt,damping_coeff_table,table_diffusion,porosity)

   real, intent(inout), dimension(is:ie,js:je) :: dt_table,table_diffusion
   real, intent(in),    dimension(is:ie,js:je) :: table
   real, intent(in)                            :: dt,damping_coeff_table,porosity

   real, dimension(size(table,1),size(table,2)) :: damp_coeff_table_grid,dt_table_test,table_test
   complex, dimension(ms:me,ns:ne)             :: dt_table_sph,table_sph,table_diffusion_sph,damp_coeff_table_sph,dt_table_test_sph

   real, dimension(size(table_sph,1),size(table_sph,2)) :: coeff
   real                                                 :: damping_order_table = 1.0
   real, allocatable, dimension(:,:)                    :: table_damping, eigen

   integer  :: num_fourier, num_spherical
   integer :: i,j

   call get_num_fourier(num_fourier)
   call get_num_spherical(num_spherical)

   allocate(table_damping (0:num_fourier, 0:num_spherical))
   allocate(eigen         (0:num_fourier, 0:num_spherical))
  

   call get_eigen_laplacian(eigen)
  
   call trans_grid_to_spherical(table,table_sph)
   call trans_grid_to_spherical(dt_table,dt_table_sph)

   table_damping(:,:) = damping_coeff_table * (eigen(:,:)**damping_order_table) 
   coeff  = 1.0/(1.0 + table_damping(ms:me,ns:ne)*dt)
   table_diffusion_sph = dt_table_sph
   dt_table_sph = coeff*(dt_table_sph - table_damping(ms:me,ns:ne)*table_sph*dt)
   table_diffusion_sph = table_diffusion_sph*(-1.) + dt_table_sph
   call trans_spherical_to_grid(table_diffusion_sph,table_diffusion)
   call trans_spherical_to_grid(dt_table_sph ,dt_table)  

!-------------------------------------------------------------------------------------
!   Attempt to do diffusion with spatially variable coefficient (transmissivity)
!-------------------------------------------------------------------------------------
     !This is SEVERELY untested; but comment off the above block starting at 
     !'table_damping(:,:) and uncomment the following block:

!   table_damping(:,:) = (eigen(:,:)**damping_order_table)
!   table_diffusion_sph = dt_table_sph
!   dt_table_sph = table_damping(ms:me,ns:ne)*table_sph  
!   table_diffusion_sph = table_diffusion_sph*(-1.) + dt_table_sph
!   call trans_spherical_to_grid(dt_table_sph,dt_table_test)
!   do i=1,size(table,1)
!     do j=1,size(table,2)
!         dt_table_test(i,j) = trans(is+i-1,js+j-1)*dt_table_test(i,j)
!     end do
!   end do
!   call trans_grid_to_spherical(dt_table_test,dt_table_sph)
!   table_damping(:,:) = (eigen(:,:)**damping_order_table) 
!   dt_table_sph = -table_damping(ms:me,ns:ne)*dt_table_sph*dt
!   call trans_spherical_to_grid(dt_table_sph ,dt_table)  


end subroutine diffuse_table

end module tam_hydrology_mod
