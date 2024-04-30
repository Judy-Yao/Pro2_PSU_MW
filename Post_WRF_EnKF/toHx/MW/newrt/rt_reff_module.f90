module rt_reff_module
  use rt_state_module
  use qn2re_wsm6
  use qn2re_gfdlfv3
  use qn2re_thompson08
  use qn2re_nssl

  implicit none
  private
  public :: calc_reff_driver
contains

  !------------------------------------------------------
  ! Driver to calculate effective radius 
  !------------------------------------------------------ 
  subroutine calc_reff_driver(method, state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    character(len=*),  intent(in) :: method 
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)
    
    reff_cloud(:) = 0.0
    reff_ice  (:) = 0.0
    reff_rain (:) = 0.0
    reff_snow (:) = 0.0
    reff_graup(:) = 0.0

    if (method == 'default') then
      call calc_reff_default     (state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    else if (method == 'mp_physics') then
      if (state%mp_physics_name == 'WSM6') then
        call calc_reff_wsm6      (state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      else if (state%mp_physics_name == 'Goddard') then  
        call calc_reff_goddard   (state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      else if (state%mp_physics_name == 'Thompson08') then  
        call calc_reff_Thompson08(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      else if (state%mp_physics_name == 'Morrison') then
        call calc_reff_morrison  (state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      else if (state%mp_physics_name == 'GFDL') then  
        call calc_reff_gfdl      (state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      else if (state%mp_physics_name == 'NSSL') then
        call calc_reff_nssl      (state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      end if
    end if
  end subroutine calc_reff_driver
  ! -----------------------


  !------------------------------------------------------
  ! Default: spherical 
  !------------------------------------------------------ 
  subroutine calc_reff_default(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = 1000.0
    reff_snow (:) = 500.0
    reff_graup(:) = 1500.0
  end subroutine calc_reff_default
  ! -----------------------


  !------------------------------------------------------
  ! WSM6
  !------------------------------------------------------ 
  subroutine calc_reff_wsm6(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = qn2re_WSM6_rain( state%qrain (x,y,:), rho_air) * 1e6
    reff_snow (:) = qn2re_WSM6_snow( state%qsnow (x,y,:), rho_air, state%temperature(x,y,:)) * 1e6
    reff_graup(:) = qn2re_WSM6_graup(state%qgraup(x,y,:), rho_air) * 1e6
  end subroutine calc_reff_wsm6
  ! -----------------------


  !------------------------------------------------------
  ! Goddard
  !------------------------------------------------------ 
  subroutine calc_reff_Goddard(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    print *, 'Goddard not implemeted'
    stop
    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = 1000.0
    reff_snow (:) = 500.0
    reff_graup(:) = 1500.0
  end subroutine calc_reff_Goddard
  ! -----------------------


  !------------------------------------------------------
  ! Thompson 
  !------------------------------------------------------ 
  subroutine calc_reff_Thompson08(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = qn2re_Thompson08_rain( state%qrain (x,y,:), rho_air, state%nrain(x,y,:)) * 1e6
    reff_snow (:) = qn2re_Thompson08_snow( state%qsnow (x,y,:), rho_air, state%temperature(x,y,:))* 1e6
    reff_graup(:) = qn2re_Thompson08_graup(state%qgraup(x,y,:), rho_air ) * 1e6
  end subroutine calc_reff_Thompson08
  ! -----------------------


  !------------------------------------------------------
  ! Morrison
  !------------------------------------------------------ 
  subroutine calc_reff_morrison(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain,reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    print *, 'Morrison not implemeted'
    stop
    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = 1000.0
    reff_snow (:) = 500.0
    reff_graup(:) = 1500.0
  end subroutine calc_reff_morrison
  ! -----------------------


  !------------------------------------------------------
  ! GFDL
  !------------------------------------------------------ 
  subroutine calc_reff_GFDL(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = qn2re_GFDLFV3_rain( state%qrain (x,y,:), rho_air) * 1e6
    reff_snow (:) = qn2re_GFDLFV3_snow( state%qsnow (x,y,:), rho_air) * 1e6
    reff_graup(:) = qn2re_GFDLFV3_graup(state%qgraup(x,y,:), rho_air) * 1e6
  end subroutine calc_reff_GFDL
  ! -----------------------


  !------------------------------------------------------
  ! NSSL
  !------------------------------------------------------ 
  subroutine calc_reff_nssl(state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x,y
    real,              intent(in) :: rho_air(:)
    real,              intent(out):: reff_cloud(:), reff_ice(:), reff_rain(:), reff_snow(:), reff_graup(:)

    reff_cloud(:) = 16.8
    reff_ice  (:) = 25.0
    reff_rain (:) = qn2re_nssl_rain( state%qrain (x,y,:), rho_air, state%nrain(x,y,:)) * 1e6
    reff_snow (:) = qn2re_nssl_snow( state%qsnow (x,y,:), rho_air, state%nsnow(x,y,:)) * 1e6
    reff_graup(:) = qn2re_nssl_graup(state%qgraup(x,y,:), rho_air, state%ngraup(x,y,:)) * 1e6
  end subroutine calc_reff_nssl
  ! -----------------------


end module rt_reff_module
