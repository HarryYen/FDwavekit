!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routine for defining velocity/attenuation structure
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_vmodel_user

  use m_std
  use m_global
  use m_geomap
  use m_fdtool
  use m_readini
  use m_seawater
  use linear_interp
  implicit none
  private
  save

  public :: vmodel_user

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !!
  !! This is a user-specific routine to define original veloicty model.
  !!
  !! Input:
  !!    io_prm                          :: I/O number of parameter file (file has been opened already)
  !!    i0,i1, j0,j1, k0,k1             :: model area by indices in i-,j-,and k- directions
  !!    xc(i0:i1), yc(j0:j1), zc(k0:k1) :: Cartesian coordinate location
  !!    vcut                            :: cut-off velocity specified by input parameter
  !!
  !! Output:
  !!    rho(k0:k1, i0:i1, j0:j1)        :: mass density (usually in g/cm^3)
  !!    lam(k0:k1, i0:i1, j0:j1)        :: Lame's parameter (usually in (g/cm^3) * (km/s)^2)
  !!    mu (k0:k1, i0:i1, j0:j1)        :: Lame's parameter (usually in (g/cm^3) * (km/s)^2)
  !!    qp (k0:k1, i0:i1, j0:j1)        :: Attenuation QP
  !!    qs (k0:k1, i0:i1, j0:j1)        :: Attenuation QS
  !!    bd (i0:i1, j0:j1, 0:NBD)        :: Boundary depths
  !!
  !! Note:
  !! bd(:,:,0) are treated as topography shape for output.
  !!    this is only for output and visualization. topography in the simulation will be automatically detected by medium params.
  !! bd(:,:,1:NBD) may contain internal boundary depths. The boundary number can be specified as source depth or station depth. 
  !!
  !<
  !! ----
  subroutine vmodel_user( io_prm, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                         !< i-region
    integer,  intent(in)  :: j0, j1                         !< j-region
    integer,  intent(in)  :: k0, k1                         !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )                  !< x-coordinate location
    real(SP), intent(in)  :: yc  ( j0:j1 )                  !< y-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )                  !< z-coordinate location
    real(SP), intent(in)  :: vcut                           !< cut-off minimum velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1, j0:j1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s)**2 ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s)**2 ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1, j0:j1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1, j0:j1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, j0:j1, 0:NBD )    !< Boundary depths
    !! --

    integer  :: i, j, k
    real(SP) :: vp0, vs0, rho0, qp0, qs0, topo0
    real(SP) :: Vp, Vs
    real(SP) :: vp1, vs1
    real(SP) :: dum
    logical  :: use_munk, earth_flattening, zero_topo
    logical  :: user_Vp, user_rho
    real(SP) :: zs(k0:k1) ! spherical depth for earth_flattening
    real(SP) :: Cv(k0:k1) ! velocity scaling coefficient for earth_flattening
    real(SP) :: clon, clat, phi, new_lat, new_lon
    character(len=60) :: vp_model_dir, vs_model_dir, rho_model_dir, topo_dir
    real(SP) :: topo_user
    integer :: line_tomo, line_topo, error
    integer :: nx, ny
    real(SP) :: dz
    !! ----

    !!
    !! The following dummy code is an example how to discribe the routine. 
    !!

    !!
    !! subroutine readini() can access parameters defined in the input file.
    !! Any original parameters can be added in the input file. 
    !!
    call readini( io_prm, 'vp0',    vp0, 5.0 )
    call readini( io_prm, 'vs0',    vs0, vp0/sqrt(3.0) )
    call readini( io_prm, 'rho0',   rho0, 2.7 )
    call readini( io_prm, 'qp0',    qp0, 1000000.0 )
    call readini( io_prm, 'qs0',    qs0, 1000000.0 )
    call readini( io_prm, 'topo0', topo0, 0.0 )
    call readini( io_prm, 'clon', clon, 0.0 )
    call readini( io_prm, 'clat', clat, 0.0 )
    call readini( io_prm, 'phi', phi, 0.0 )
    call readini( io_prm, 'topo_dir', topo_dir, '/home/data/data4/OpenSWPC-5.2.0/src/3D_topo/3dtopo_re.xyz')
    call readini( io_prm, 'nx', nx, 1000 )
    call readini( io_prm, 'ny', ny, 1000 )
    call readini( io_prm, 'dz', dz, 0.1 )
    !call interp(Model_p, Model_s, x, y, z)
    
    !! seawater
    call readini( io_prm, 'munk_profile', use_munk, .false. )
    call seawater__init( use_munk )
    !! earth-flattening tranformation
    !! if this option is true, zs(:) array is nonlinearly mapped from evenly-spaced 
    !! zc(:). Use zs(:) to set velocity models. Please note that P and S wave velocities
    !! should be multiplied Cv(k) which depends on depth. 
    call readini( io_prm, 'earth_flattening', earth_flattening, .false. )
    if( earth_flattening ) then
      do k=k0, k1
        zs(k) = R_EARTH - R_EARTH * exp( - zc(k) / R_EARTH )
        Cv(k) = exp( zc(k) / R_EARTH)
      end do
    else
      zs(:) = zc(:)
      Cv(:) = 1.0
    end if

    
    call readini( io_prm, 'user_Vp', user_Vp, .false. )
    call readini( io_prm, 'user_rho', user_rho, .false. )

    if( user_Vp ) then
      call readini( io_prm, 'vp_model_dir', vp_model_dir, '/home/harry' )
      open(2019, file=vp_model_dir, status = 'OLD', form = 'formatted',&
           access = 'direct', recl=6)
    end if
    
    if( user_rho ) then  
      call readini( io_prm, 'rho_model_dir', rho_model_dir, '/home/harry' )
      open(2020, file=rho_model_dir, status = 'OLD', form = 'formatted',&
           access = 'direct', recl=6)
    end if

    call readini( io_prm, 'vs_model_dir', vs_model_dir, '/home/harry')
    open(2021, file=vs_model_dir, status = 'OLD', form = 'formatted',&
         access = 'direct', recl=6)
    
    !! read topo file by harry
    call readini( io_prm, 'zero_topo', zero_topo, .true. )
    if ( .not. ( zero_topo ) ) then
      open(2022, file=topo_dir, status = 'OLD', form = 'formatted',&
      access = 'direct', recl=11)
    end if    

    !! The medium parameter must be set from given region (i0:i1, j0:j1, k0:k1)
    !! Note that the order of indices is k->i->j, for improving performance
    !!
    do j = j0, j1
      do i = i0, i1

        !! define topography shape here
        if ( zero_topo ) then
          bd(i,j,0) = topo0
        else  
          line_topo = (j + 2) * (nx + 6) + (i + 3)
          read(2022, fmt="(F10.3)", rec=line_topo, IOSTAT=error)topo_user
          bd(i,j,0) = topo_user/1000.
        end if 
        
        do k = k0, k1
          
          if( zs( k ) > bd(i,j,0) ) then
            line_tomo = (j+2)*(nx+6)*(nz+6) + (i+2)*(nz+6) + (k+3)
            
            read(2021,fmt="(F5.3)", rec=line_tomo, IOSTAT=error)Vs
            if ( user_Vp ) then
              read(2019,fmt="(F5.3)", rec=line_tomo, IOSTAT=error)Vp
            else     
              Vp = Vs * sqrt(3.0)
              !Vp = 0.9409+2.0947*Vs-0.8206*(Vs**2)+0.2683*(Vs**3)-0.0251*(Vs**4)
            end if
            !! elastic medium
            vp1 = Cv(k) * Vp
            vs1 = Cv(k) * Vs
            
            if ( user_rho ) then
              read(2020,fmt="(F5.3)", rec=line_tomo, IOSTAT=error)rho(k,i,j)
            else 
              rho(k,i,j) = 1.6612*Vp-0.4721*Vp**2+0.0671*Vp**3-0.0043*Vp**4+0.000106*Vp**5
            end if

            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = qp0
            qs (k,i,j) = qs0

            !! added by Hung-Yu 2022/09/19
            !! also for ocean column
               
            if( bd(i,j,0) > 0 .AND. 0 < zs( k ) .AND. zs(k) <= dz ) then    
              vp1 = Cv(k) * seawater__vel( zc(k) )
              vs1 = 0.0
              rho(k,i,j) = 1.0
              mu (k,i,j) = rho(k,i,j) * vs1 * vs1
              lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
              qp (k,i,j) = 1000000.0 ! effectively no attenuation in ocean column
              qs (k,i,j) = 1000000.0
            end if          
            

          else if ( zs (k) >= 0.0 ) then

            !!
            !! ocean column
            !!
            !! The code treat the uppermost layer as ocean column if P-wave velocity is finite and S-wave velocity is zero
            !!
            vp1 = Cv(k) * seawater__vel( zc(k) )
            vs1 = 0.0
            rho(k,i,j) = 1.0
            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = 1000000.0 ! effectively no attenuation in ocean column
            qs (k,i,j) = 1000000.0

          else

            !!
            !! air column
            !!
            !! The air column must have zero P- & S-wave velocity (i.e., mu=lam=0)
            !! Please use non-zero but very small density (e.g., 0.001) for avoiding zero division with satisfying boundary cond.
            !! Since waves do not penetrate to the air column, qp and qs does not affect. Just set dummy.
            !!
            vp1 = 0.0
            vs1 = 0.0
            rho(k,i,j) = 0.001
            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
            qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column

          end if
        end do
      end do
    end do

    !! dummy value
    bd(:,:,1:NBD) = -9999

    ! substitute to a dummy variable for avoiding compiler warnings
    dum = xc(i0)
    dum = yc(j0)
    dum = zc(k0)
    dum = vcut

  end subroutine vmodel_user
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_vmodel_user
!! ----------------------------------------------------------------------------------------------------------------------------- !!
