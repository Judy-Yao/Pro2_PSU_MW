!=======================================================================================
subroutine xb_to_surface(inputfile,proj,xa,ix,jx,kx,nv,iob,xland,lu_index,znu,znw,p_top,times,xb)
use constants
use namelist_define
use obs_define
use netcdf
use wrf_tools
use map_utils
implicit none
character(len=10), intent(in)           :: inputfile
character(len=19), intent(in)           :: times
character(len=10) :: obstype
type(proj_info), intent(in)             :: proj                   ! 1st guest map info
real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
integer, intent(in)                     :: ix, jx, kx, nv, iob
real, dimension(ix, jx ), intent(in)    :: xland, lu_index
real,  intent(in)                       :: p_top
real, dimension(kx), intent(in)         :: znu
real, dimension(kx+1), intent(in)       :: znw
real, intent(out)                       :: xb
real, dimension(ix, jx, kx+1)           :: ph, phb
real, dimension(ix, jx, kx  )           :: pt, qv, qc, qr, pb
real, dimension(ix, jx      )           :: mu, mub, t2m, th2m, q2m, u10m, v10m
real, dimension(ix+1, jx, kx)           :: u
real, dimension(ix, jx+1, kx)           :: v
real, dimension(2, 2, kx+1)             :: p
real, dimension(kx)                     :: pres, ptt, qvt, ht
real                                    :: mu1, mub1, long, grid_u, grid_v
integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk, obs_ii,obs_jj
integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_var, i_u, i_v
integer                                 :: i_t2, i_th2, i_q2, i_u10, i_v10
real                                    :: dx, dxm, dy, dym, dz, dzm
real                                    :: psfc, tsfc, psfcm, u10, v10, t2, q2, th2
real, dimension(ix, jx      )           :: rough
obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

call roughness_from_landuse ( 'USGS', times, ix, jx, lu_index, rough ) 

!- calculate q, pressure profile on (obs_ii, obs_jj)
i_qv = 0 
i_qc = 0 
i_qr = 0 
i_mu = 0
i_mub = 0
i_ph = 0
i_phb = 0
i_pt = 0
i_u = 0
i_v = 0
i_t2 = 0
i_th2 = 0
i_q2 = 0
i_u10 = 0
i_v10 = 0
do m = 1, nv
   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
   if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
   if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
   if ( enkfvar(m) == 'MU        ' ) i_mu=m
   if ( enkfvar(m) == 'MUB       ' ) i_mub=m
   if ( enkfvar(m) == 'PH        ' ) i_ph=m
   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
   if ( enkfvar(m) == 'T         ' ) i_pt=m
   if ( enkfvar(m) == 'U         ' ) i_u=m
   if ( enkfvar(m) == 'V         ' ) i_v=m
   if ( enkfvar(m) == 'T2        ' ) i_t2=m
   if ( enkfvar(m) == 'TH2       ' ) i_th2=m
   if ( enkfvar(m) == 'Q2        ' ) i_q2=m
   if ( enkfvar(m) == 'U10       ' ) i_u10=m
   if ( enkfvar(m) == 'V10       ' ) i_v10=m
enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_qc>0) qc (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qc )
if(i_qr>0) qr (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qr )
if(i_mu>0)  mu (i1:i1+1, j1:j1+1)   = xa( 1:2, 1:2, 1, i_mu )
if(i_mub>0) mub (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_mub )
if(i_t2>0)  t2m (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_t2 )
if(i_th2>0) th2m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_th2 )
if(i_q2>0)  q2m (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_q2 )
if(i_u10>0) u10m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_u10 )
if(i_v10>0) v10m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_v10 )
if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1,    mu  )
if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1,    mub )
if ( i_t2== 0  ) call get_variable2d(inputfile, 'T2        ', ix, jx, 1,    t2m )
if ( i_th2== 0 ) call get_variable2d(inputfile, 'TH2       ', ix, jx, 1,    th2m)
if ( i_q2== 0  ) call get_variable2d(inputfile, 'Q2        ', ix, jx, 1,    q2m )
if ( i_u10== 0 ) call get_variable2d(inputfile, 'U10       ', ix, jx, 1,    u10m)
if ( i_v10== 0 ) call get_variable2d(inputfile, 'V10       ', ix, jx, 1,    v10m)

qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)

!...... get total qv at obs' position from horizontal interpolation
qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
!...... get mu,mub at obs' position from horizontal interpolation
mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))
!...... calculate pressure profile from qv
call cal_press_from_q( kx, znu, znw, qvt, mu1, mub1, p_top, pres )

!- calculate t (not theta) and height profile
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_phb )
if(i_pt>0) pt (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_pt )
if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, pt )

!...... get geopotential height profile around obs%position, then horizontal interpolated to obs's position
p(1:2,1:2,1:3) = (ph(i1:i1+1, j1:j1+1, 1:3) + phb(i1:i1+1, j1:j1+1, 1:3))/g
p(1,1,1:3) = dym*(dx*p(2,1,1:3) + dxm*p(1,1,1:3)) + dy*(dx*p(2,2,1:3) + dxm*p(1,2,1:3))
ht(1:2) = 0.5*(p(1,1,1:2)+p(1,1,2:3))
!...... get ptt(theta)  at obs' position from horizontal interpolation
ptt(1) = dym*(dx*pt(i1+1,j1,  1) + dxm*pt(i1,j1,  1)) + dy*(dx*pt(i1+1,j1+1,1) + dxm*pt(i1,j1+1,1))
ptt(1) = theta_to_temp(ptt(1)+to, pres(1))
 
!-------------------------------------------------
!- calculate surface p and ground t 
   psfc = pres(1) + (1. - znu(1))*mub1
   tsfc = ptt(1) + 0.0065*(ht(1)-p(1,1,1))
!
!-------------------------------------------------
!...... get model u and v
if(i_u>0) u(i1:i1+2,j1:j1+1,1)=xa(1:3,1:2,1,i_u)
if(i_v>0) v(i1:i1+1,j1:j1+2,1)=xa(1:2,1:3,1,i_v)
if ( i_u == 0 ) call get_variable3d(inputfile, 'U         ', ix+1, jx, kx, 1, u )
if ( i_v == 0 ) call get_variable3d(inputfile, 'V         ', ix, jx+1, kx, 1, v )

!...... horizontal interp for U
if ( obs_ii-i1 == 0.500 ) then
   grid_u = u(i1+1,j1,1)*(j1+1-obs_jj) + u(i1+1,j1+1,1)*(obs_jj-j1)
else if ( obs_ii-i1 > 0.500 ) then
   grid_u = (j1+1-obs_jj)*( u(i1+1,j1  ,1)*(i1+1.5-obs_ii)+u(i1+2,j1  ,1)*(obs_ii-i1-0.5) ) + &
             (obs_jj-j1)  *( u(i1+1,j1+1,1)*(i1+1.5-obs_ii)+u(i1+2,j1+1,1)*(obs_ii-i1-0.5) )
else if ( obs_ii-i1 < 0.500 ) then
   grid_u = (j1+1-obs_jj)*( u(i1,j1  ,1)*(i1+0.5-obs_ii)+u(i1+1,j1  ,1)*(obs_ii-i1+0.5) ) + &
             (obs_jj-j1)  *( u(i1,j1+1,1)*(i1+0.5-obs_ii)+u(i1+1,j1+1,1)*(obs_ii-i1+0.5) )
endif

!...... horizontal interp for V
if ( obs_jj-j1 == 0.500 ) then
   grid_v = v(i1,j1+1,1)*(i1+1-obs_ii) + v(i1+1,j1+1,1)*(obs_ii-i1)
else if ( obs_jj-j1 > 0.500 ) then
   grid_v = (i1+1-obs_ii)*( v(i1  ,j1+1,1)*(j1+1.5-obs_jj) + v(i1  ,j1+2,1)*(obs_jj-j1-0.5) ) + &
             (obs_ii-i1)  *( v(i1+1,j1+1,1)*(j1+1.5-obs_jj) + v(i1+1,j1+2,1)*(obs_jj-j1-0.5) )
else if ( obs_jj-j1 < 0.500 ) then
   grid_v = (i1+1-obs_ii)*( v(i1  ,j1,1)*(j1+0.5-obs_jj) + v(i1  ,j1+1,1)*(obs_jj-j1+0.5) ) + &
             (obs_ii-i1)  *( v(i1+1,j1,1)*(j1+0.5-obs_jj) + v(i1+1,j1+1,1)*(obs_jj-j1+0.5) )
endif

!-------------------------------------------------
!- calculate 10m wind, 2m t and q
!call sfc_wtq( psfc, tsfc, pres(1), ptt(1), qvt(1), grid_u, grid_v,          &
!              pres(2), ptt(2), qvt(2), ht(1), rough(i1,j1), xland(i1,j1),   &
!              u10, v10, t2, q2 )  
t2  = dym*(dx*t2m (i1+1,j1) + dxm*t2m (i1,j1)) + dy*(dx*t2m (i1+1,j1+1) + dxm*t2m (i1, j1+1))
th2 = dym*(dx*th2m(i1+1,j1) + dxm*th2m(i1,j1)) + dy*(dx*th2m(i1+1,j1+1) + dxm*th2m(i1, j1+1))
q2  = dym*(dx*q2m (i1+1,j1) + dxm*q2m (i1,j1)) + dy*(dx*q2m (i1+1,j1+1) + dxm*q2m (i1, j1+1))
u10 = dym*(dx*v10m(i1+1,j1) + dxm*u10m(i1,j1)) + dy*(dx*u10m(i1+1,j1+1) + dxm*u10m(i1, j1+1))
v10 = dym*(dx*u10m(i1+1,j1) + dxm*v10m(i1,j1)) + dy*(dx*v10m(i1+1,j1+1) + dxm*v10m(i1, j1+1))


!-------------------------------------------------
!- Correct surface pressure
call da_sfc_pre ( psfcm, psfc, t2, q2, p(1,1,1), obs%sta(iob,1), obs%sta(iob,3), obs%sta(iob,4)/1000.)

!if( print_detail > 100 )write(*,'(3x,a,5f)')'model u =', u(i1,j1  ,1), u(i1+1,j1  ,1), u(i1+2,j1  ,1)  
!if( print_detail > 100 )write(*,'(3x,a,5f)')'model u =', u(i1,j1+1,1), u(i1+1,j1+1,1), u(i1+2,j1+1,1)  
!if( print_detail > 100 )write(*,'(3x,a,5f)')'grid_u, grid_v, u10, v10 =',grid_u, grid_v, u10, v10
!if( print_detail > 100 )write(*,'(3x,a,4f)')'station elevation, p, t, q =', obs%sta(iob,1:4)
!if( print_detail > 100 )write(*,'(3x,a,5f)')'t2, q2, model_terrain =', t2, q2, p(1,1,1)
!if( print_detail > 100 )write(*,'(3x,a,5f)')'psfc and corrected psfc =', psfc, psfcm
!-------------------------------------------------
!- get xb
!...... surface pressure
if ( obstype(10:10) == 'P' ) then
     xb = psfcm
     xb = psfc
else if ( obstype(10:10) == 'U' ) then
     xb = u10
else if ( obstype(10:10) == 'V' ) then
     xb = v10
else if ( obstype(10:10) == 'T' ) then
     xb = t2
else if ( obstype(10:10) == 'Q' ) then
     xb = q2*1000.
else if ( obstype(10:10) == 'D' ) then
     xb = mixrat_to_tdew(q2, psfcm)
else if ( obstype(10:10) == 'R' ) then
     xb = rel_humidity(q2, t2, psfcm)
else if ( obstype(10:10) == 'S' ) then   !wind speed
     xb = sqrt(u10**2.+v10**2.)
endif

!if( print_detail > 100 )write(*,'(3x,a,3f)')'xb_to_surface '//obstype//' obs xb, obs_ii, obs_jj :', xb, obs_ii, obs_jj

end subroutine xb_to_surface

!=======================================================================================
   subroutine xb_to_sounding (inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,xb,itruth,isimulated )
   use constants
   use namelist_define
   use obs_define
   use netcdf
   use wrf_tools
   use map_utils
   use mpi_module

   implicit none

   character(len=10), intent(in)           :: inputfile
   type(proj_info), intent(in)             :: proj                   ! 1st guest map info
   real, dimension(3,3,kx+1,nv), intent(in)  :: xa                     ! 1st guest
   integer, intent(in)                     :: ix, jx, kx, nv, iob  
   real, dimension(ix, jx ), intent(in)    :: xlong
   real, dimension(kx+1), intent(in)       :: znw
   real, dimension(kx), intent(in)         :: znu
   real,  intent(in)                       :: p_top
   integer, intent(in)                     :: itruth, isimulated

   real :: obs_ii, obs_jj, obs_kk, obs_pp !
   real, intent(out)   :: xb
   character(len=10)   :: obstype

   real, dimension(ix, jx, kx+1)           :: ph, phb
   real, dimension(ix, jx, kx  )           :: pt, pb, qv, qc, qr
   real, dimension(ix, jx      )           :: mu, mub
   real, dimension(ix+1, jx, kx)           :: u
   real, dimension(ix, jx+1, kx)           :: v
   real, dimension(kx+1)                   :: znw0
   real, dimension(kx)                     :: znu0

   real, dimension(2, 2, kx+1)             :: p
   real, dimension(kx+1)                   :: work, worku, workv

   real, dimension(kx)                     :: pres, ptt, qvt
   real                                    :: mu1, mub1, long, grid_u, grid_v, true_u, true_v
   real                                    :: dx, dxm, dy, dym, dz, dzm
   integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk
   integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_var, i_u, i_v

   obstype=obs%type(iob)
   obs_ii=obs%position(iob,1)
   obs_jj=obs%position(iob,2)
   obs_kk=obs%position(iob,3)
   obs_pp=obs%position(iob,4)
   i1 = int( obs_ii )
   j1 = int( obs_jj ) 
   dx  = obs_ii-real(i1)
   dxm = real(i1+1)-obs_ii
   dy  = obs_jj-real(j1)
   dym = real(j1+1)-obs_jj

!.. calculate pressure from geopotential height
    if ( itruth == 0 )then
         call get_variable1d(inputfile, 'ZNW       ', kx+1, 1, znw0)
         call get_variable1d(inputfile, 'ZNU       ', kx  , 1, znu0)
    else
         znw0 = znw 
         znu0 = znu
    endif 

!.. Calculate obs%position(iob,3)
    if ( isimulated == 0 .or. obstype(10:10) == 'T' .or. obstype(10:10) == 'D' .or.  &
                              obstype(10:10) == 'R' .or. obstype(10:10) == 'Q' ) then
   
!..    get data from background
       i_ph = 0
       i_phb = 0
       i_mu = 0
       i_mub = 0
       i_pt = 0
       i_qv = 0
       do m = 1, nv
          if ( enkfvar(m) == 'PH        ' ) i_ph=m
          if ( enkfvar(m) == 'PHB       ' ) i_phb=m
          if ( enkfvar(m) == 'T         ' ) i_pt=m
          if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
          if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
          if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
          if ( enkfvar(m) == 'MU        ' ) i_mu=m
          if ( enkfvar(m) == 'MUB       ' ) i_mub=m
       enddo
       if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
       if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
       if(i_pt>0) pt(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pt)
       if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
       if(i_qc>0) qc(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qc)
       if(i_qr>0) qr(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qr)
       if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
       if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)

       if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
       if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
       if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, pt )
       if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
       if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
       if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
       if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1,    mu )
       if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1,    mub)

       qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)

!...... get mut mu, mub at obs' position from horizontal interpolation
       mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
       mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))

!........ get qvt (qvapor) at obs' position from horizontal interpolation
!........ get ptt(theta)  at obs' position from horizontal interpolation
       qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
       ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))

       p(1:2,1:2,1:kx+1) = ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1)
       p(1,1,1:kx+1) = dym*(dx*p(2,1,1:kx+1) + dxm*p(1,1,1:kx+1)) + dy*(dx*p(2,2,1:kx+1) + dxm*p(1,2,1:kx+1))

       call eta_to_pres(znw0(1:kx+1), mu1+mub1, qvt(1:kx), p(1,1,1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))
!   if( print_detail==3 )write(*,'(3x,a,4f)')'obs location =', obs%position(iob,1:4)
!   if( print_detail==3 )write(*,'(3x,a,f)')'          mu =', mu(1,1)
!   if( print_detail==3 )write(*,'(3x,a,10f)')'          p  =', p(1,1,1:10)
!       call cal_press_from_q( kx, znu0, znw0, qvt, mu1, mub1, p_top, pres )
!       if( print_detail > 100 ) write(*,'(a,100f10.0)')'Pressure = ', pres
!       if( print_detail > 100 ) write(*,'(a,100f10.4)')'qvt = ', qvt

       if ( isimulated == 0 ) then
          if ( obstype(1:1) == 'H' ) then
!........... get geopotential height profile around obs%position, then horizontal interpolated to obs's position
             call to_zk(obs%position(iob,4), p(1,1,1:kx)/g, obs%position(iob,3), kx) 
          else if ( obstype(1:1) == 'P' ) then
             call to_zk(obs%position(iob,4), pres(1:kx), obs%position(iob,3), kx)
          endif
!          if ( obs%position(iob,3) .ge. real(kx-1) ) obs%position(iob,3) = real(kx-1)
          if ( obs%position(iob,3) .lt. 1. ) obs%position(iob,3) = 1.
          obs_kk = obs%position(iob,3)
       endif
    endif

!.. Get xb from background
   if ( obs_kk .gt. real(kx-1) .or. obs_kk .lt. 1. ) then
      xb = -999999.
      return
   endif

   k1  = int( obs_kk )
   dz  = obs_kk-real(k1)
   dzm = real(k1+1)-obs_kk

!.. U, V
    if ( obstype(10:10) == 'U' .or. obstype(10:10) == 'V' .or. obstype(10:10) == 'S' ) then
       i_u = 0
       i_v = 0
       do m = 1, nv
          if ( enkfvar(m) == 'U         ' ) i_u=m
          if ( enkfvar(m) == 'V         ' ) i_v=m
       enddo
       u(i1:i1+2,j1:j1+1,k1:k1+1)=xa(1:3,1:2,k1:k1+1,i_u)
       v(i1:i1+1,j1:j1+2,k1:k1+1)=xa(1:2,1:3,k1:k1+1,i_v)

       worku = -88888.
       workv = -88888.
       if ( obs_ii-i1 == 0.500 ) then
          worku(k1:k1+1) = dym*u(i1+1,j1,k1:k1+1) + dy*u(i1+1,j1+1,k1:k1+1)
       else if ( obs_ii-i1 > 0.500 ) then
          worku(k1:k1+1) = dym*( u(i1+1,j1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1,k1:k1+1)*(dx-0.5) ) + &
                           dy *( u(i1+1,j1+1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1+1,k1:k1+1)*(dx-0.5) )
       else if ( obs_ii-i1 < 0.500 ) then
          worku(k1:k1+1) = dym*( u(i1,j1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1,k1:k1+1)*(dx+0.5) ) + &
                           dy *( u(i1,j1+1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1+1,k1:k1+1)*(dx+0.5) )
       endif
       if ( obs_jj-j1 == 0.500 ) then
          workv(k1:k1+1) = v(i1,j1+1,k1:k1+1)*dxm + v(i1+1,j1+1,k1:k1+1)*dx
       else if ( obs_jj-j1 > 0.500 ) then
          workv(k1:k1+1) = dxm*( v(i1,j1+1,k1:k1+1)*(dym+0.5) + v(i1,j1+2,k1:k1+1)*(dy-0.5) ) + &
                           dx *( v(i1+1,j1+1,k1:k1+1)*(dym+0.5) + v(i1+1,j1+2,k1:k1+1)*(dy-0.5) )
       else if ( obs_jj-j1 < 0.500 ) then
          workv(k1:k1+1) = dxm*( v(i1,j1,k1:k1+1)*(dym-0.5) + v(i1,j1+1,k1:k1+1)*(dy+0.5) ) + &
                           dx *( v(i1+1,j1,k1:k1+1)*(dym-0.5) + v(i1+1,j1+1,k1:k1+1)*(dy+0.5) )
       endif

!       if( print_detail > 100 ) write(*,*)'i1,j1,k1, v=',i1,j1,k1,v(i1,j1,k1),v(i1,j1+1,k1),v(i1,j1+2,k1),v(i1+1,j1,k1)
!       if( print_detail > 100 ) write(*,'(a,100f10.2)')'U profile =', worku(k1:k1+1)
!       if( print_detail > 100 ) write(*,'(a,100f10.2)')'V profile =', workv(k1:k1+1)

       if ( obs_kk .le. 1. ) then
          grid_u = worku(k1)
          grid_v = workv(k1)
       else
          grid_u = dzm*worku(k1)+dz*worku(k1+1)
          grid_v = dzm*workv(k1)+dz*workv(k1+1)
       endif
!       if( print_detail > 100 ) write(*,'(a,f10.2)')'grid U =', grid_u, 'grid V =', grid_v

!------011108 changed obs. wind to gridwind
       long = ( xlong(i1  ,j1)*dym + xlong(i1  ,j1+1)*dy ) * dxm +          & 
              ( xlong(i1+1,j1)*dym + xlong(i1+1,j1+1)*dy ) * dx
       call gridwind_to_truewind( long, proj, grid_u, grid_v, true_u, true_v )
    
       if ( obstype(10:10) == 'U' ) then
           xb = grid_u
       else if ( obstype(10:10) == 'V' ) then
           xb = grid_v
       else if ( obstype(10:10) == 'S' ) then   !wind speed
           xb = sqrt(grid_u**2.+grid_v**2.)
       endif

!       if( print_detail > 100 ) then
!           write(*,'(12f10.3)')obs_kk, grid_u, grid_v, dz, dzm, long, worku(k1:k1+1), workv(k1:k1+1), true_u, true_v
!       endif

!.. T
    else if ( obstype(10:10) == 'T' ) then
       do k = k1, k1+1
          work(k) = theta_to_temp(ptt(k)+to, pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(1)
       else
          xb = dzm*work(k1)+dz*work(k1+1)
       endif

!.. TD(if used TD, not RH)
    else if ( obstype(10:10) == 'D' ) then
       do k = k1, k1+1
          work(k) = mixrat_to_tdew(qvt(k), pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(1)
       else
          xb = dzm*work(k1) + dz*work(k1+1)
       endif

!.. RH
    else if ( obstype(10:10) == 'R' )then
       do k = k1, k1+1
          work(k) = theta_to_temp(ptt(k)+to, pres(k))
       enddo
       do k = k1, k1+1
          work(k) = rel_humidity(qvt(k), work(k), pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(1)
       else
          xb = dzm*work(k1) + dz*work(k1+1)
       endif
       if ( xb > 100. )xb = 100.

!.. Q
    else if ( obstype(10:10) == 'Q' )then
       if ( obs_kk .le. 1. ) then
          xb = qvt(k1)*1000.
       else
          xb = (dzm*qvt(k1)+dz*qvt(k1+1))*1000.
       endif
       if ( xb .le. 0.0 ) xb = -99999.

!.. HG
    else if ( obstype(10:10) == 'H' )then
       call destag_zstag(znu0, znw0, kx, p(1,1,1:kx+1), work(1:kx))
       if ( obs_kk .le. 1. ) then
          xb = work(1)/g
       else
          xb = (dzm*work(k1) + dz*work(k1+1))/g
       endif 

    endif

!   if( print_detail > 100 )write(*,'(a,3f8.2, f10.1)')'xb_to_sounding '//obstype//' obs position :', obs_ii, obs_jj, obs%position(iob,3:4)
    
   end subroutine xb_to_sounding
!
!=======================================================================================
   subroutine xb_to_idealsound(inputfile,xa,ix,jx,kx,nv,iob,xb)
   use constants
   use namelist_define
   use obs_define 
   use netcdf 
   use wrf_tools
   use map_utils 

   implicit none

   character(len=10), intent(in)           :: inputfile
   real, dimension(3,3,kx+1,nv), intent(in) :: xa 
   integer, intent(in) :: ix, jx, kx, nv, iob
   real :: obs_ii, obs_jj, obs_kk 
   character(len=10) :: obstype
   real, intent(out) :: xb
   integer :: i1, j1, k1, i, j, k, m, ii, jj, kk

   obstype=obs%type(iob)
   obs_ii=obs%position(iob,1)
   obs_jj=obs%position(iob,2)
   obs_kk=obs%position(iob,3)
   i1 = int( obs_ii )
   j1 = int( obs_jj )
   k1 = int( obs_kk )
    
   if (obstype .eq. 'idealU    ' ) then
      do m = 1, nv
         if ( enkfvar(m) == 'U         ' ) then
            xb = xa(i1,j1,k1,m)
         endif
      enddo
   else if (obstype .eq. 'idealV    ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'V         ' ) then
            xb = xa(i1,j1,k1,m)
         endif
      enddo 
   else if (obstype .eq. 'idealPT   ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'T         ' ) then
            xb = xa(i1,j1,k1,m) 
         endif
      enddo
   else if (obstype .eq. 'idealQV   ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'QVAPOR    ' ) then
            xb = xa(i1,j1,k1,m)*1000.
         endif
      enddo
   else if (obstype .eq. 'idealPH   ' ) then
      xb = 0.
      do m = 1, nv
         if ( enkfvar(m) == 'PH        ' .or. enkfvar(m) == 'PHB       ' ) then
            xb = xb + xa(i1,j1,k1,m)
         endif
      enddo
   endif

   end subroutine xb_to_idealsound

!=======================================================================================
  subroutine xb_to_rv(inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znw,xb,kkflag)
  use constants 
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils 
  use netcdf
  use wrf_tools
  use radar
  implicit none
  character (len=10), intent(in) :: inputfile
  type(proj_info), intent(in) :: proj
  integer, intent(in)         :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in)  :: xa
  real, dimension(ix, jx ), intent(in) :: xlong
  real, dimension(kx+1), intent(in)    :: znw
  real, intent(out)                    :: xb
  integer, intent(in)                  :: kkflag
  real  :: obs_ii, obs_jj, obs_kk, obs_hh, radar_ii, radar_jj, radar_elv

  real, dimension(ix+1, jx, kx)        :: u
  real, dimension(ix, jx+1, kx)        :: v
  real, dimension(ix, jx, kx+1)        :: w, ph, phb
  real, dimension(ix, jx, kx  )        :: t, qv, p, pb
  real, dimension(ix, jx      )        :: mu, mub, terrain
    
  integer :: m, ii, jj, kk, i, j, k, i1, j1
  integer :: i_u, i_v, i_w, i_t, i_qv, i_mu, i_ph, i_phb, i_mub, i_pb

  obs_ii=obs%position(iob,1)
  obs_jj=obs%position(iob,2)
  obs_kk=obs%position(iob,3)
  obs_hh=obs%position(iob,4)

  radar_ii=obs%sta(iob,1)
  radar_jj=obs%sta(iob,2)
  radar_elv=obs%sta(iob,4)

  i1 = int( obs_ii )
  j1 = int( obs_jj )
    
  i_u  = 0 ; i_v   = 0 ; i_w  = 0 ; i_t = 0 ; i_qv = 0 ; i_mu = 0; i_mub = 0;
  i_ph = 0 ; i_phb = 0 ; i_pb = 0 ;

  do m=1,nv
     if ( enkfvar(m) == 'U         ' ) i_u=m
     if ( enkfvar(m) == 'V         ' ) i_v=m
     if ( enkfvar(m) == 'W         ' ) i_w=m
  enddo
  if(i_u>0) u(i1:i1+2,j1:j1+1,1:kx)=xa(1:3,1:2,1:kx,i_u)
  if(i_v>0) v(i1:i1+1,j1:j1+2,1:kx)=xa(1:2,1:3,1:kx,i_v)
  if(i_w>0) w(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_w)
  if ( i_u == 0 ) call get_variable3d(inputfile, 'U         ', ix+1, jx, kx, 1, u )
  if ( i_v == 0 ) call get_variable3d(inputfile, 'V         ', ix, jx+1, kx, 1, v )
  if ( i_w == 0 ) call get_variable3d(inputfile, 'W         ', ix, jx, kx+1, 1, w )
  
  if ( kkflag == 0 ) then
     do m = 1, nv
        if ( enkfvar(m) == 'T         ' ) i_t=m
        if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
        if ( enkfvar(m) == 'PH        ' ) i_ph=m
        if ( enkfvar(m) == 'PHB       ' ) i_phb=m
        if ( enkfvar(m) == 'PB        ' ) i_pb=m
        if ( enkfvar(m) == 'MU        ' ) i_mu=m
        if ( enkfvar(m) == 'MUB       ' ) i_mub=m
     end do
     if(i_t>0) t(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_t)
     if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
     if(i_pb>0) pb(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pb)
     if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
     if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
     if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
     if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)
     if( i_t <1  ) call get_variable3d(inputfile, 'T         ', ix, jx, kx, 1, t )
     if( i_qv <1 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx, 1, qv )
     if( i_pb <1 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx, 1, pb )
     if( i_ph <1 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
     if( i_phb <1) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb )
     if( i_mu <1 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1, mu )
     if( i_mub <1) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1, mub )
     do j = j1, j1+1
     do i = i1, i1+1
        call cal_ph( kx, znw, t(i,j,1:kx), qv(i,j,1:kx), pb(i,j,1:kx), mu(i,j), mub(i,j), phb(i,j,1:kx+1), ph(i,j,1:kx+1) )
     enddo
     enddo

!     if (print_detail > 1000) write(*,'(a,70f8.0)') 'phb+ph = ',(phb(i1,j1,1:kx+1)+ph(i1,j1,1:kx+1))/g
!     if (print_detail > 1000) write(*,'(a,f8.0)') 'obs_hh =', obs_hh
     call calc_radar_data_point_kk ( ix, jx, kx, ph, phb, obs_ii, obs_jj, obs_hh, obs%position(iob,3) )
     obs_kk=obs%position(iob,3)
!     if (print_detail > 1000) write(*,'(a, 4f8.2)')'obs_ijhk =', obs_ii, obs_jj, obs_hh, obs_kk
  endif

  call calculate_rv ( ix, jx, kx, u, v, w, xlong, proj, obs_ii, obs_jj, obs_kk, obs_hh, radar_ii, radar_jj, radar_elv, xb )
!  write(*,*)'i,j,k =', obs_ii, obs_jj, obs_kk
!  write(*,*)'u,v,w =', u(int(obs_ii),int(obs_jj),int(obs_kk)), v(int(obs_ii),int(obs_jj),int(obs_kk)), w(int(obs_ii),int(obs_jj),int(obs_kk))

  end subroutine xb_to_rv

!=======================================================================================
subroutine xb_to_pw(inputfile,xa,ix,jx,kx,nv,iob,znu,znw,xb)
  use constants
  use namelist_define
  use obs_define
  use netcdf
  use wrf_tools
  implicit none
  character(len=10), intent(in)           :: inputfile
  character(len=10) :: obstype
  integer, intent(in)                     :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, intent(out)                       :: xb
  real, dimension(2, 2)                   :: pw,temp,vtemp,zdiff
  real, dimension(ix+1, jx, kx  )         :: u
  real, dimension(ix, jx+1, kx  )         :: v
  real, dimension(ix, jx, kx  )           :: t,p,pb,z,qv
  real, dimension(ix, jx, kx+1)           :: ph,phb
  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
  integer                                 :: i_t, i_u, i_v, i_ph, i_phb, i_qv, i_p, i_pb
  real                                    :: obs_ii, obs_jj, dx,dxm,dy,dym

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

i_qv = 0 
i_ph = 0
i_phb = 0
i_t = 0
i_p = 0
i_pb = 0
do m = 1, nv
   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
   if ( enkfvar(m) == 'P         ' ) i_p=m
   if ( enkfvar(m) == 'PB        ' ) i_pb=m
   if ( enkfvar(m) == 'PH        ' ) i_ph=m
   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
   if ( enkfvar(m) == 'T         ' ) i_t=m
enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_phb )
if(i_p>0) p (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_p )
if(i_pb>0) pb (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_pb )
if(i_t>0) t (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_t )
     
if ( i_t  == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
if ( i_p  == 0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

ph(i1:i1+1, j1:j1+1, 1:kx+1) = (ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1))/g
p(i1:i1+1, j1:j1+1, 1:kx)  = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)
do j=j1,j1+1
do i=i1,i1+1
  z(i, j, 1:kx) = (ph(i, j, 1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i, j, 2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
enddo
enddo
do j=j1,j1+1
do i=i1,i1+1
do k=1,kx
  if( p(i,j,k)/100. >= 1200. ) write(*,'(a,3i4,2f8.0)')'P error at:',i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
enddo
enddo
enddo
pw = 0.
do k=1,kx-1
  temp=(t(i1:i1+1,j1:j1+1,k)+300.)*(p(i1:i1+1,j1:j1+1,k)/Pr)**(rd/cp)
  vtemp=(1+0.61*qv(i1:i1+1,j1:j1+1,k))*temp
  zdiff=z(i1:i1+1,j1:j1+1,k+1)-z(i1:i1+1,j1:j1+1,k)
  pw=pw+p(i1:i1+1,j1:j1+1,k)/(rd*vtemp)*zdiff*qv(i1:i1+1,j1:j1+1,k)
enddo
xb = dym*(dx*pw(2,1) + dxm*pw(1,1)) + dy*(dx*pw(2,2) + dxm*pw(1,2))
xb = xb/10

end subroutine xb_to_pw

!=======================================================================================
subroutine xb_to_slp(inputfile,xa,ix,jx,kx,nv,iob,znu,znw,xb )

!---------------------
! slp_center subroutine calculates sea level pressure and find the hurricane center
!---------------------
  use constants
  use namelist_define
  use obs_define
  use netcdf
  use wrf_tools

  implicit none

  character(len=10), intent(in)           :: inputfile
  character(len=10) :: obstype
  integer, intent(in)                     :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, intent(out)                       :: xb

  real, dimension(2, 2)                   :: slp
  real, dimension(ix+1, jx, kx  )         :: u
  real, dimension(ix, jx+1, kx  )         :: v
  real, dimension(ix, jx, kx  )           :: t
  real, dimension(ix, jx, kx+1)           :: ph
  real, dimension(ix, jx, kx+1)           :: phb
  real, dimension(ix, jx, kx  )           :: p
  real, dimension(ix, jx, kx  )           :: z
  real, dimension(ix, jx, kx  )           :: pb
  real, dimension(ix, jx, kx  )           :: qv, qc, qr
  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
  integer                                 :: i_t, i_u, i_v, i_ph, i_phb, i_qv, i_qc, i_qr, i_p, i_pb  ! variables flag
  real                                    :: obs_ii, obs_jj, dx,dxm,dy,dym

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

i_qv = 0 
i_qc = 0 
i_qr = 0 
i_ph = 0
i_phb = 0
i_t = 0
i_p = 0
i_pb = 0
do m = 1, nv
   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
   if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
   if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
   if ( enkfvar(m) == 'P         ' ) i_p=m
   if ( enkfvar(m) == 'PB        ' ) i_pb=m
   if ( enkfvar(m) == 'PH        ' ) i_ph=m
   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
   if ( enkfvar(m) == 'T         ' ) i_t=m
enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_qc>0) qc (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qc )
if(i_qr>0) qr (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qr )
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_phb )
if(i_p>0) p (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_p )
if(i_pb>0) pb (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_pb )
if(i_t>0) t (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_t )
     
if ( i_t  == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
if ( i_p  == 0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)
ph(i1:i1+1, j1:j1+1, 1:kx+1) = (ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1))/g
p(i1:i1+1, j1:j1+1, 1:kx)  = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)

do j=j1,j1+1
do i=i1,i1+1
  z(i, j, 1:kx) = (ph(i, j, 1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i, j, 2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
enddo
enddo
do j=j1,j1+1
do i=i1,i1+1
do k=1,kx
  if( p(i,j,k)/100. >= 1200. ) write(*,'(a,3i4,2f8.0)')'P error at:', i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
enddo
enddo
enddo
   call compute_seaprs(2, 2, kx, z(i1:i1+1,j1:j1+1,:), t(i1:i1+1,j1:j1+1,:), p(i1:i1+1,j1:j1+1,:), qv(i1:i1+1,j1:j1+1,:), slp, 0)

   xb = dym*(dx*slp(2,1) + dxm*slp(1,1)) + dy*(dx*slp(2,2) + dxm*slp(1,2))
   xb = xb*100.

end subroutine xb_to_slp






!=======================================================================================
subroutine xb_to_radiance(inputfile,proj,ix,jx,kx,xlong,xlat,landmask,iob_radmin,iob_radmax,xb_tb,use_cloud)

!---------------------
! radiance subroutine calculates brightness temperature for satellite channels
!---------------------

  USE constants
  USE netcdf
  USE mpi_module
  USE CRTM_Module
  use namelist_define
  use obs_define
  use wrf_tools


  implicit none

  integer, intent(in)                      :: ix, jx, kx
  integer, intent(out)                     :: iob_radmin,iob_radmax
  character(len=10), intent(in)            :: inputfile
  type(proj_info), intent(in)              :: proj                   ! 1st guestmap info
  logical, intent(in)                      :: use_cloud
  real, dimension(obs%num), intent(out)    :: xb_tb
  real, dimension(ix, jx ), intent(in)     :: xlong
  real, dimension(ix, jx ), intent(in)     :: xlat
  real, dimension(ix, jx ), intent(in)     :: landmask
  integer                                  :: iob,irad
  real                                     :: obs_ii, obs_jj, dx,dxm,dy,dym

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
  REAL, PARAMETER :: Re=6378000.0
  !====================
  !setup for geostationary satellites
   REAL, PARAMETER :: sat_h=35780000.0
   REAL :: sat_lon
  !====================
!  INTEGER, intent(in) :: ix = ix  !total number of the x-grid
!  INTEGER, parameter, intent(in) :: jx = jx  !total number of the y-grid
!  INTEGER, parameter, intent(in) :: kx = kx        !level range
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1 
!  INTEGER, PARAMETER :: N_LAYERS    = kx
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
!  INTEGER, PARAMETER :: N_CLOUDS    = kx*5
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis

  ! Variables
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  CHARACTER(256) :: FILE_NAME
  CHARACTER(256) :: obstype
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: l, m, irec, yend, ystart, nyi
  integer :: ncid,ncrcode
  character(LEN=16) :: var_name
  character(LEN=3)  :: file_ens
  integer :: x,y,tt,v,z,n,reci,ens,n_ec,num_radgrid
  INTEGER :: ncl,icl,k1,k2
  integer :: lat_radiance(ix*jx)  ! latitude
  integer :: lon_radiance(ix*jx) ! longitude
  real :: lat(ix,jx)   ! in radian
  real :: lon(ix,jx)   ! in radian
  real :: p(ix,jx,kx)
  real :: pb(ix,jx,kx)
  real :: pres(ix,jx,kx)
  real :: ph(ix,jx,kx+1)
  real :: phb(ix,jx,kx+1)
  real :: delz(kx)
  real :: t(ix,jx,kx)
  real :: tk(ix,jx,kx)
  real :: qvapor(ix,jx,kx)
  real :: qcloud(ix,jx,kx)
  real :: qrain(ix,jx,kx)
  real :: qice(ix,jx,kx)
  real :: qsnow(ix,jx,kx)
  real :: qgraup(ix,jx,kx)
  real :: psfc(ix,jx)
  real :: hgt(ix,jx)
  real :: tsk(ix,jx)
  real, allocatable, dimension(:,:,:) :: Tbsend, Tb

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)
  ! ============================================================================

  ! ============================================================================
  ! 1.5. **** make a loop to get the number of satellite-radiance-iob ****
  !
  num_radgrid = 0
  check_cycle:do iob=1,obs%num
    obstype = obs%type(iob)
    if ( obstype(1:8) == 'Radiance' ) then
     if(num_radgrid == 0) then
      num_radgrid = num_radgrid + 1
      lon_radiance(num_radgrid) = int(obs%position(iob,1))
      lat_radiance(num_radgrid) = int(obs%position(iob,2))
      iob_radmin = iob
      iob_radmax = iob
     else
      iob_radmax = iob
      do irad = 1,num_radgrid
      if((lon_radiance(irad)==int(obs%position(iob,1))).and.(lat_radiance(irad)==int(obs%position(iob,2))))cycle check_cycle
      enddo
      num_radgrid = num_radgrid + 1
      lon_radiance(num_radgrid) = int(obs%position(iob,1))
      lat_radiance(num_radgrid) = int(obs%position(iob,2))
     endif
    endif
  enddo check_cycle


  ! ============================================================================
  ! --------------
  CALL CRTM_Version( Version )
  !if(my_proc_id==0)  write(*,*) "CRTM ver.",TRIM(Version) 
  ! Get sensor id from user
  ! -----------------------
  !It assumes that all the Radiance data is same sattelite as the first data.
  Sensor_Id = trim(adjustl(obs%sat(iob_radmin)))
 
  if (Sensor_Id == 'ahi_h8' .or. Sensor_Id == 'ahi_himawari8') then
     sat_lon = 140.0/180.0*3.14159
  else if (Sensor_Id == 'abi_gr' .or. Sensor_Id == 'imgr_g13') then
     sat_lon = -75.2/180.0*3.14159
  endif
  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. This initializes the CRTM for the sensors
  !     predefined in the example SENSOR_ID parameter.
  !     NOTE: The coefficient data file path is hard-
  !           wired for this example.
  ! --------------------------------------------------
  !if(my_proc_id==0) WRITE( *,'(/5x,"Initializing the CRTM...")' )
  if (my_proc_id==0 .AND. (inputfile == 'fort.80011' .OR. inputfile == 'fort.80071')) then
    Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
                              ChannelInfo  , &  ! Output
                              IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                              IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                              File_Path='coefficients/',&
                              Quiet=.false.)
  else
    Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
                              ChannelInfo  , &  ! Output
                              IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                              IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                              File_Path='coefficients/',&
                              Quiet=.true.)
  endif
  
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  ! Specify channel 14 for GOES-R ABI
  if (Sensor_Id == 'abi_gr' .or. Sensor_Id == 'ahi_h8' ) then
    Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset =(/8,9,10/) )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error initializing CRTM'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
  endif 
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  allocate(Tb(ix,jx,n_Channels))
  allocate(Tbsend(ix,jx,n_Channels))
  Tb = 0.
  Tbsend = 0.

  ! ============================================================================

  ! ============================================================================
  ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 3a. Allocate the ARRAYS
  ! -----------------------
  ! Note that only those structure arrays with a channel
  ! dimension are allocated here because we've parameterized
  ! the number of profiles in the N_PROFILES parameter.
  !
  ! Users can make the 
  ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  ! 3b. Allocate the STRUCTURES
  ! ---------------------------
  ! The input FORWARD structure
  CALL CRTM_Atmosphere_Create( Atm, kx, N_ABSORBERS, kx*5, N_AEROSOLS)
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! Fill the Atm structure array.
  ! NOTE: This is an example program for illustrative purposes only.
  !       Typically, one would not assign the data as shown below,
  !       but rather read it from file
  !
  ! 4a1. Loading Atmosphere and Surface input
  ! --------------------------------
  call get_variable3d(inputfile,'P',ix,jx,kx,1,p)
  call get_variable3d(inputfile,'PB',ix,jx,kx,1,pb)
  call get_variable3d(inputfile,'PH',ix,jx,kx+1,1,ph)
  call get_variable3d(inputfile,'PHB',ix,jx,kx+1,1,phb)
  call get_variable3d(inputfile,'T',ix,jx,kx,1,t)
  call get_variable3d(inputfile,'QVAPOR',ix,jx,kx,1,qvapor)
  call get_variable2d(inputfile,'PSFC',ix,jx,1,psfc)
  call get_variable2d(inputfile,'TSK',ix,jx,1,tsk)
  call get_variable2d(inputfile,'HGT',ix,jx,1,hgt)
  if (use_cloud) then
     call get_variable3d(inputfile,'QCLOUD',ix,jx,kx,1,qcloud)
     call get_variable3d(inputfile,'QRAIN',ix,jx,kx,1,qrain)
     call get_variable3d(inputfile,'QICE',ix,jx,kx,1,qice)
     call get_variable3d(inputfile,'QSNOW',ix,jx,kx,1,qsnow)
     call get_variable3d(inputfile,'QGRAUP',ix,jx,kx,1,qgraup)
  else
     qcloud = 0.
     qrain = 0.
     qice = 0.
     qsnow = 0.
     qgraup = 0.
  endif

  lat = xlat/180.0*3.14159
  lon = xlong/180.0*3.14159
  pres = P + PB
  tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )
  where(qvapor.lt.0.0) qvapor=1.0e-8
  where(qcloud.lt.0.0) qcloud=0.0
  where(qice.lt.0.0) qice=0.0
  where(qrain.lt.0.0) qrain=0.0
  where(qsnow.lt.0.0) qsnow=0.0
  where(qgraup.lt.0.0) qgraup=0.0

  ! 4a2. Parallerization with grids
  ! --------------------------------
  !--- preparation for the x,y-loop
  if(mod(num_radgrid,nprocs).eq.0) then
     nyi=num_radgrid/nprocs
  else
     nyi=num_radgrid/nprocs+1
  endif
  ystart=my_proc_id*nyi+1
  yend=min(num_radgrid,(my_proc_id+1)*nyi)

  do iob = ystart, yend
     obs_ii=lon_radiance(iob)
     obs_jj=lat_radiance(iob)
     x = int( obs_ii )
     y = int( obs_jj )

  ! 4a3. Converting WRF data for CRTM structure
  ! --------------------------------
  !--- converting the data to CRTM structure

  !*******************************************************************************
  ! satellite information
  !*******************************************************************************

  sat_dis=sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
  SCAN_ANGLE=180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))**2))
  ZENITH_ANGLE=SCAN_ANGLE+180.0/3.14159*acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))

  !*******************************************************************************
  ! load WRF data into CRTM structures
  !*******************************************************************************
  !--- calcurating delz
  do z=1,kx
   if(z.eq.1) then
    delz(z) = (PH(x,y,z+1) + PHB(x,y,z+1)) / 9.806 - hgt(x,y)
   else
    !delz(z) = ((PH(x,y,z+1) + PHB(x,y,z+1))-(PH(x,y,z) + PHB(x,y,z)))/2/9.806  ! corrected 9 July 2020 after group meeting
    delz(z) = ((PH(x,y,z+1) + PHB(x,y,z+1))-(PH(x,y,z) + PHB(x,y,z)))/9.806
   endif
  enddo
  if (delz(1) <= 0.) delz(1) = delz(2)
  !---Atmospheric Profile
  atm(1)%Climatology         = TROPICAL
  atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
  atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)
  atm(1)%Level_Pressure(0) = (pres(x,y,kx)*3.0/2.0 - pres(x,y,kx-1)/2.0)/100.0  ! convert from Pa to hPA
  do z=kx,1,-1
    if(z.eq.1) then
      atm(1)%Level_Pressure(kx-z+1) = max(psfc(x,y), pres(x,y,1)*3.0/2.0-pres(x,y,2)/2.0)/100.0
    else
      atm(1)%Level_Pressure(kx-z+1) = ((pres(x,y,z-1)+pres(x,y,z))/2.0)/100.0  ! convert from Pa to hPA
    endif
    atm(1)%Pressure(kx-z+1)       = pres(x,y,z) / 100.0
    atm(1)%Temperature(kx-z+1)    = tk(x,y,z)
    atm(1)%Absorber(kx-z+1,1)     = qvapor(x,y,z)*1000.0
  enddo
  atm(1)%Absorber(:,2) = 5.0E-02 
  !---Cloud Profile
  do z=1,kx*5
   atm(1)%Cloud(z)%Type = 0
   atm(1)%Cloud(z)%Effective_Radius = 0.0
   atm(1)%Cloud(z)%Water_Content = 0.0
  enddo
  ncl = 0
  icl = 0
  !--calculating # of clouds (cloud and rain)
  do z=kx,1,-1
    if(qcloud(x,y,z).gt.0.0) then
      ncl = ncl + 1
    endif
    if(qrain(x,y,z).gt.0.0) then
      ncl = ncl + 1
    endif
    if(qice(x,y,z).gt.0.0) then
      ncl = ncl + 1
    endif
    if(qsnow(x,y,z).gt.0.0) then
      ncl = ncl + 1
    endif
    if(qgraup(x,y,z).gt.0.0) then
      ncl = ncl + 1
    endif
  enddo
  !--Data for cloud
  atm(1)%n_Clouds         = ncl
  IF ( atm(1)%n_Clouds > 0 ) THEN
  do z=kx,1,-1
    if(qcloud(x,y,z).gt.0.0) then
      icl = icl + 1
      k1 = kx-z+1
      k2 = kx-z+1
      atm(1)%Cloud(icl)%Type = WATER_CLOUD
      atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 16.8_fp
      atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
          qcloud(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
    endif
  enddo
  do z=kx,1,-1
    if(qrain(x,y,z).gt.0.0) then
      icl = icl + 1
      k1 = kx-z+1
      k2 = kx-z+1
      atm(1)%Cloud(icl)%Type = RAIN_CLOUD
      atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1000.0_fp
      atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
          qrain(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
    endif
  enddo
  do z=kx,1,-1
    if(qice(x,y,z).gt.0.0) then
      icl = icl + 1
      k1 = kx-z+1
      k2 = kx-z+1
      atm(1)%Cloud(icl)%Type = ICE_CLOUD
      atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 25.0_fp
      atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
          qice(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
    endif
  enddo
  do z=kx,1,-1
    if(qsnow(x,y,z).gt.0.0) then
      icl = icl + 1
      k1 = kx-z+1
      k2 = kx-z+1
      atm(1)%Cloud(icl)%Type = SNOW_CLOUD
      atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 750.0_fp
      atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
          qsnow(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
    endif
  enddo
  do z=kx,1,-1
    if(qgraup(x,y,z).gt.0.0) then
      icl = icl + 1
      k1 = kx-z+1
      k2 = kx-z+1
      atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
      atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1500.0_fp
      atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
          qgraup(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
    endif
  enddo
  ENDIF

  !---Surface data
  if(landmask(x,y).eq.1.0) then
   sfc(1)%Water_Coverage = 0.0_fp
   sfc(1)%Land_Coverage = 1.0_fp
   sfc(1)%Land_Temperature = tsk(x,y)
   sfc(1)%Soil_Temperature = tsk(x,y)
  else
   sfc(1)%Water_Coverage = 1.0_fp
   sfc(1)%Land_Coverage = 0.0_fp
   sfc(1)%Water_Type = 1  ! Sea water
   sfc(1)%Water_Temperature = tsk(x,y)
  endif


  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! 4c. Use the SOI radiative transfer algorithm
  ! --------------------------------------------
  Options%RT_Algorithm_ID = RT_SOI
  ! ============================================================================

  ! ============================================================================
  ! 5. **** CALL THE CRTM FORWARD MODEL ****
  !
  Error_Status = CRTM_Forward( Atm        , &
                               Sfc        , &
                               Geometry   , &
                               ChannelInfo, &
                               RTSolution , &
                               Options = Options )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================



  ! ============================================================================
  ! 6. **** Collecting output ****
  !
  ! User should read the user guide or the source code of the routine
  ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
  ! select the needed variables for outputs.  These variables are contained
  ! in the structure RTSolution.
  !
  !DO m = 1, N_PROFILES
  !  WRITE( *,'(//7x,"Profile ",i0," output for ",a )') n, TRIM(Sensor_Id)
  !  DO l = 1, n_Channels
  !    WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
  !    CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
  !  END DO
  !END DO

  !---for file output, edited 2014.9.26
  do l = 1, n_Channels
      Tbsend(x,y,l) = real(RTSolution(l,1)%Brightness_Temperature)
  enddo

  !--- end of iob(x,y)-loop
  enddo

  CALL MPI_Allreduce(Tbsend,Tb,ix*jx*n_Channels,MPI_REAL,MPI_SUM,comm,ierr)

  ! ============================================================================

  ! ============================================================================
  !6.5  **** writing the output ****
  !
  if(my_proc_id==0) then
  do iob = iob_radmin, iob_radmax
     obs_ii=obs%position(iob,1)
     obs_jj=obs%position(iob,2)
     x = int( obs_ii )
     y = int( obs_jj )
     if (Sensor_Id == 'abi_gr' .or. Sensor_Id == 'ahi_h8' ) then
         if (obs%ch(iob) .eq. 8) xb_tb(iob) = Tb(x,y,1) !6.19um
         if (obs%ch(iob) .eq. 9) xb_tb(iob) = Tb(x,y,2) !6.95um
         if (obs%ch(iob) .eq. 10) xb_tb(iob) = Tb(x,y,3) !7.34um
         if (obs%ch(iob) .eq. 14) write(*,*)'change channel setting for ch14' !xb_tb(iob) = Tb(x,y,4) !11.2um
     elseif (Sensor_Id == 'imgr_g13' ) then
         if (obs%ch(iob) .eq. 3) xb_tb(iob) = Tb(x,y,2) !6.19um
         if (obs%ch(iob) .eq. 4) xb_tb(iob) = Tb(x,y,3) !11.2um
     else
         do l = 1, n_Channels
            if (obs%ch(iob) .eq. RTSolution(l,1)%Sensor_Channel) xb_tb(iob) = Tb(x,y,l) 
         enddo
     endif
  enddo
  !--initializing the Tbsend fields for Bcast
  Tbsend = 0.0
  endif
  if(my_proc_id==0) &
   WRITE(*,'(a10," Tb=",f6.2,"~",f6.2)')inputfile,minval(xb_tb(iob_radmin:iob_radmax)),maxval(xb_tb(iob_radmin:iob_radmax))

  ! ============================================================================
  !  **** initializing all Tb and Tbsend fields ****
  !
  Tb = 0.0
  CALL MPI_BCAST(Tbsend,ix*jx*n_Channels,MPI_REAL,0,comm,ierr)

  ! ============================================================================
  ! 7. **** DESTROY THE CRTM ****
  !
  deallocate(Tb)
  deallocate(Tbsend)
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================

end subroutine xb_to_radiance










!=======================================================================================
subroutine xb_to_microwave(inputfile,proj,ix,jx,kx,xlong,xlat,landmask,times,iob_radmin,iob_radmax,xb_tb,write_file_override,use_cloud)

!---------------------
! radiance subroutine calculates brightness temperature for satellite channels
!---------------------

  USE constants
  USE netcdf
  USE mpi_module
  USE CRTM_Module
  use namelist_define
  use obs_define
  use wrf_tools


  implicit none

  integer, intent(in)                      :: ix, jx, kx
  integer, intent(out)                     :: iob_radmin,iob_radmax
  character(len=10), intent(in)            :: inputfile
  type(proj_info), intent(in)              :: proj                   ! 1st guestmap info
  logical, intent(in)                      :: use_cloud
  real, dimension(obs%num), intent(out)    :: xb_tb
  real, dimension(ix, jx ), intent(in)     :: xlong
  real, dimension(ix, jx ), intent(in)     :: xlat
  real, dimension(ix, jx ), intent(in)     :: landmask
  character(len=80), intent(in)            :: times
  logical, intent(in)                      :: write_file_override
  integer                                  :: iob,isensor,n_sensors_used,ich,ich2,ich3
  real                                     :: obs_ii, obs_jj, dx,dxm,dy,dym

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
  INTEGER, PARAMETER :: N_SUPPORTED_SENSOR_CLASSES = 8
  CHARACTER(*), PARAMETER, DIMENSION(N_SUPPORTED_SENSOR_CLASSES) :: SENSOR_CLASS_NAMES = &
      (/'ssmis','ssmi','gmi_gpm_hf','gmi_gpm_lf','saphir_meghat','amsr2_gcom-w1','atms_npp','mhs'/)
  INTEGER, DIMENSION(N_SUPPORTED_SENSOR_CLASSES) :: N_CHANNELS = &
      (/24, 7, 13, 13, 6, 14, 22,5/)
  !INTEGER, PARAMETER :: N_MAX_CH_CLASSES = 7
  !REAL(fp), PARAMETER, DIMENSION(N_SUPPORTED_SENSOR_CLASSES) :: SCAN_ANGLE = &
  !    (/53.1,52.8/)
  !REAL(fp), PARAMETER, DIMENSION(N_SUPPORTED_SENSOR_CLASSES) :: ZENITH_ANGLE = &
  !    (/45.0,48.5/)
  !INTEGER, DIMENSION(N_MAX_CHS,N_SUPPORTED_SENSOR_CLASSES) :: SENSOR_SATCONV_EFOV_CLASS = &
  !  (/ (/5,5,5,5,5,5,5,3,3,3,3,1,1,1,2,2,4,4,6,6,6,6,6,5/),   &
  !     (/1,1,2,2,3,4,4,5,5,6,6,7,7,0,0,0,0,0,0,0,0,0,0,0/) /)
  ! index 1 is EFOV along-scan and index 2 is EFOV perpendicular to scan
  ! gmi EFOV from Petty and Bennartz 2017 (corrigendum)
  !REAL, DIMENSION(2,N_MAX_CH_CLASSES,N_SUPPORTED_SENSOR_CLASSES)  :: SENSOR_SATCONV_EFOV = &
  !  (/ (/ (/44.8,73.6/),(/27.5,45.0/),(/13.2,15.5/),(/13.2,15.5/),(/17.6,27.3/),(/17.6,27.3/),(/00.0,00.0/) /),   &
  !     (/ (/19.8,32.1/),(/11.7,18.1/),(/10.5,16.0/),(/10.3,15.6/),(/ 6.4, 7.2/),(/ 5.8, 6.3/),(/ 5.6, 5.8/) /) /)


  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
  REAL, PARAMETER :: Re=6378000.0
  REAL, PARAMETER :: MY_PI = 4*ATAN(1.)
  REAL, PARAMETER :: RADS_PER_DEGREE = MY_PI / 180.
  !====================
!  INTEGER, intent(in) :: ix = ix  !total number of the x-grid
!  INTEGER, parameter, intent(in) :: jx = jx  !total number of the y-grid
!  INTEGER, parameter, intent(in) :: kx = kx        !level range
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1 
!  INTEGER, PARAMETER :: N_LAYERS    = kx
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
!  INTEGER, PARAMETER :: N_CLOUDS    = kx*5
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
!  INTEGER, PARAMETER :: N_STREAMS = -1  ! let CRTM decide automatically
!  INTEGER, PARAMETER :: N_STREAMS = 4
!  INTEGER, PARAMETER :: N_STREAMS = 6
!  INTEGER, PARAMETER :: N_STREAMS = 8
  INTEGER, PARAMETER :: N_STREAMS = 16
  INTEGER, PARAMETER :: x_coarse = 1, y_coarse = 1  ! how many wrf gridpoints should be averaged together in a CRTM profile

  ! Variables
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_ID
  CHARACTER(256) :: CRTM_OUT_DIR
  CHARACTER(256) :: FILE_NAME
  CHARACTER(4)   :: YEAR_STR
  CHARACTER(2)   :: MONTH_STR, DAY_STR
  CHARACTER(3)   :: DOMAIN
  CHARACTER(256) :: obstype
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: sum_n_channels            ! total number of channels across all included sensors
  INTEGER :: n_Channels_all            ! total number of channels with observations across all sensors
  INTEGER :: n_Channels_crtm           ! total number of channels with observations across all sensors that need CRTM calculations
  INTEGER :: l, m, irec, grand_count, yend, ystart, xend, xstart
  integer :: ncid,ncrcode,fid
  character(LEN=16) :: var_name
  character(LEN=3)  :: file_ens
  integer :: x,y,i,j,tt,v,z,n,reci,ens,n_ec,num_radgrid
  INTEGER :: ncl,icl,k1,k2
  real :: sigma, search_radius
  real :: cputime
  integer, allocatable, dimension(:) :: lat_radiance_min, lat_radiance_max ! latitude
  integer, allocatable, dimension(:) :: lon_radiance_min, lon_radiance_max ! longitude
  integer :: search_lat_min, search_lat_max
  integer :: search_lon_min, search_lon_max
  real :: beam_conv_simple, beam_conv_gaussian_simple
  real(fp) :: my_scan_angle, my_zenith_angle, my_azimuth_angle, my_azimuth_sine, my_azimuth_cosine, my_azimuth_angle_est
  integer :: my_angle_sum_count, scan_angle_search_radius
  INTEGER :: YEAR, MONTH, DAY
  

  logical, DIMENSION(:), allocatable :: is_channel_used, is_channel_crtm   ! logical array, length sum_n_channels

  integer, DIMENSION(:), allocatable :: channel_numbers_crtm, sensor_index_crtm, channel_search   ! length n_Channels_crtm

  integer, DIMENSION(:,:), allocatable :: channel_numbers_all_eachSensor   ! size n_sensors_used X n_Channels_all (using a total of only n_Channels_crtm elements)
  integer, allocatable, dimension(:) :: n_channels_all_eachSensor  ! length number of sensors
  INTEGER :: n_Channels_all_mySensor   ! total number of channels with observations for the sensor of interest
  integer, allocatable, dimension(:) :: channel_numbers_all_mySensor  ! length n_channels_crtm_mySensor
  INTEGER :: mySensor_Tb_offset  ! offset of channel numbers in the Tb array for sensor of interest

  integer, DIMENSION(:,:), allocatable :: channel_numbers_crtm_eachSensor   ! size n_sensors_used X n_Channels_crtm (using a total of only n_Channels_crtm elements)
  integer, allocatable, dimension(:) :: n_channels_crtm_eachSensor  ! length number of sensors
  integer :: n_Channels_crtm_mySensor  ! total number of channels with observations for the sensor of interest that need CRTM calculations
  integer, allocatable, dimension(:) :: channel_numbers_crtm_mySensor  ! length n_channels_crtm_mySensor

  integer, allocatable, dimension(:) :: sensors

  logical :: crtm_out_files
  !integer :: numx_crtm, numy_crtm
  integer, allocatable, dimension(:) :: numx_crtm, numy_crtm  ! length number of sensors
  real :: lat(ix,jx)   ! in radian
  real :: lon(ix,jx)   ! in radian
  real :: p(ix,jx,kx)
  real :: pb(ix,jx,kx)
  real :: pres(ix,jx,kx)
  real :: ph(ix,jx,kx+1)
  real :: phb(ix,jx,kx+1)
  real :: delz(kx)
  real :: t(ix,jx,kx)
  real :: tk(ix,jx,kx)
  real :: qvapor(ix,jx,kx)
  real :: qcloud(ix,jx,kx)
  real :: qrain(ix,jx,kx)
  real :: qice(ix,jx,kx)
  real :: qsnow(ix,jx,kx)
  real :: qgraup(ix,jx,kx)
  real :: nrain(ix,jx,kx)
  real :: psfc(ix,jx)
  real :: hgt(ix,jx)
  real :: tsk(ix,jx)
  real :: xland(ix,jx)
  real :: u10(ix,jx)
  real :: v10(ix,jx)
  real :: windspeed(ix,jx)
  real :: westwind(ix,jx)
  real :: winddir(ix,jx)
  real, allocatable, dimension(:,:,:) :: Tbsend, Tb_crtm, Tb  ! size # sensors X spatial grid
  real, allocatable, dimension(:,:,:) :: scan_angle, zenith_angle, azimuth_angle, azimuth_sine, azimuth_cosine  ! size # sensors X spatial grid
  real :: azimuth_angle_theta
  integer, allocatable, dimension(:,:,:) :: angle_sum_count  ! size # sensors X spatial grid
  logical :: is_sensor_used(N_SUPPORTED_SENSOR_CLASSES)
  integer, allocatable, dimension(:) :: sensors_used, sensor_ch_offset  ! length number of sensors
  CHARACTER(256), allocatable, dimension(:) :: Sensor_IDs_used
  CHARACTER(256) :: CRTM_Init_Sensor_ID
  integer, allocatable, dimension(:) :: sensor_indices   ! length number of obs
  integer :: my_sensor_index
  integer :: max_n_channels_all_sensor
  integer :: scan_angle_search_number
  integer :: scan_angle_search_lonmin, scan_angle_search_lonmax
  integer :: scan_angle_search_latmin, scan_angle_search_latmax

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)

  integer :: mp_scheme = 0  ! wsm6, gfdlfv3, thompson08
  logical :: l_use_default_re = .false.
  integer, parameter :: N_CLOUDS=5

  ! ============================================================================

  ! ============================================================================
  ! 1.5. **** make a loop to get the number of satellite-microwave-iob ****
  !      **** and lat-lon bounds for CRTM calculations                 ****
  !

  if(my_proc_id==0) then
    call cpu_time(cputime)
    write(*,*) 'entered xb.f at ', cputime
  endif
  

  call open_file(inputfile, nf_nowrite, fid)
  ncrcode = nf_get_att_real(fid, nf_global,'DX', dx)
  ncrcode = nf_get_att_real(fid, nf_global,'DY', dy)
  ! convert to km
  dx = dx / 1000.
  dy = dy / 1000.
  call close_file(fid)


  ! determine the number of unique sensors of observation data in the file
  ! determine the total number of microwave observations
  is_sensor_used = .FALSE.
  num_radgrid = 0
  check_cycle1:do iob=1,obs%num
    obstype = obs%type(iob)
    if ( obstype(1:9) == 'Microwave' ) then
      if(num_radgrid == 0) then
        iob_radmin = iob
      endif
      iob_radmax = iob
      num_radgrid = num_radgrid + 1

      ! Get sensor id from observation object
      Sensor_ID = trim(adjustl(obs%sat(iob)))
      do isensor = 1,N_SUPPORTED_SENSOR_CLASSES
        if(INDEX(trim(Sensor_ID), trim(SENSOR_CLASS_NAMES(isensor)))>0) exit
      enddo
  
      is_sensor_used(isensor) = .TRUE.

    endif
  enddo check_cycle1
  n_sensors_used = count(is_sensor_used)
  if (my_proc_id==0) write(*,*) 'n_sensors_used:', n_sensors_used

  ! determine number of channels across all sensors with observations
  allocate(sensors_used(n_sensors_used), STAT = Allocate_Status)
  allocate(sensor_ch_offset(n_sensors_used), STAT = Allocate_Status)
  sensor_ch_offset = 0
  sum_n_channels = 0
  max_n_channels_all_sensor = 0
  i = 1
  do isensor = 1,N_SUPPORTED_SENSOR_CLASSES
    if (is_sensor_used(isensor)) then
      sum_n_channels = sum_n_channels + N_CHANNELS(isensor)
      max_n_channels_all_sensor = max(max_n_channels_all_sensor, N_CHANNELS(isensor))
      sensors_used(i) = isensor   ! will result in sensor indices of increasing value from 1 to number of sensors
      if (i < n_sensors_used) then
        sensor_ch_offset(i+1) = sum(sensor_ch_offset(1:i)) + N_CHANNELS(isensor)
        i = i + 1
      else
        exit
      endif
    endif
  enddo
  if(my_proc_id==0) write(*,*) 'sensors_used:',sensors_used


  ! allocate arrays. 
  ! the channels are arranged by the order of sesnors in SENSOR_CLASS_NAMES

  allocate(lon_radiance_min(n_sensors_used), STAT = Allocate_Status)
  allocate(lon_radiance_max(n_sensors_used), STAT = Allocate_Status)
  allocate(lat_radiance_min(n_sensors_used), STAT = Allocate_Status)
  allocate(lat_radiance_max(n_sensors_used), STAT = Allocate_Status)
  lon_radiance_min = ix
  lon_radiance_max = 0
  lat_radiance_min = jx
  lat_radiance_max = 0

  allocate(scan_angle(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(zenith_angle(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(azimuth_angle(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(azimuth_sine(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(azimuth_cosine(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(angle_sum_count(n_sensors_used,ix,jx), STAT = Allocate_Status)
  scan_angle = 0.d0
  zenith_angle = 0.d0
  azimuth_angle = 0.d0
  azimuth_sine = 0.d0
  azimuth_cosine = 0.d0
  angle_sum_count = 0

  allocate(is_channel_used(sum_n_channels), STAT = Allocate_Status)
  allocate(is_channel_crtm(sum_n_channels), STAT = Allocate_Status)
  is_channel_used = .FALSE.
  is_channel_crtm = .FALSE.

  allocate(sensor_indices(num_radgrid), STAT = Allocate_Status)

  allocate(Sensor_IDs_used(n_sensors_used), STAT = Allocate_Status)
 
  ! determine locations of observations for each sensor
  ! get Sensor_ID for each sensor used
  check_cycle2:do i=1,num_radgrid
    iob = i + iob_radmin - 1

    ! Get sensor id index for each observation, save them for easy reference later
    Sensor_ID = trim(adjustl(obs%sat(iob)))
    do isensor = 1,N_SUPPORTED_SENSOR_CLASSES 
      if(INDEX(trim(Sensor_ID), trim(SENSOR_CLASS_NAMES(isensor)))>0) exit
    enddo
    do my_sensor_index = 1,n_sensors_used
      if(sensors_used(my_sensor_index) .EQ. isensor) exit
    enddo
    sensor_indices(i) = my_sensor_index
    Sensor_IDs_used(my_sensor_index) = trim(Sensor_ID)

    ! identify that this channel will be used
    is_channel_used(obs%ch(iob) + sensor_ch_offset(my_sensor_index)) = .TRUE.
    
    ! determine satellite beam convolution area for the observation
    !if (my_proc_id == 0) write(*,*) 'about to calculate sigma:\n', iob, obs%efov_aScan(iob), obs%efov_cScan(iob)
    sigma = (0.5 * ( ( ( ( obs%efov_aScan(iob) + obs%efov_cScan(iob) ) /2 ) / 1.18) ) )
    search_radius = 1 + (2 * sigma)

    search_lon_min = max(1, nint(obs%position(iob,1)) - ceiling(1+(search_radius / dx)))
    search_lon_max = min(ix,nint(obs%position(iob,1)) + ceiling(1+(search_radius / dx)))
    search_lat_min = max(1, nint(obs%position(iob,2)) - ceiling(1+(search_radius / dy)))
    search_lat_max = min(jx,nint(obs%position(iob,2)) + ceiling(1+(search_radius / dy)))
    !if(my_proc_id==0 )  write(*,*) "search_lon: ", search_lon_min, ' to ', search_lon_max
    !if(my_proc_id==0 )  write(*,*) "search_lat: ", search_lat_min, ' to ', search_lat_max
    !if(my_proc_id==0 )  write(*,*) "obs%scan_angle(iob): ", obs%scan_angle(iob)
    !if(my_proc_id==0 )  write(*,*) "obs%zenith_angle(iob): ", obs%zenith_angle(iob)

    ! include specifications for this observation in determining the optimal scan and zenith angles for CRTM simulations
    scan_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    scan_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + abs(obs%scan_angle(iob))
    !if(my_proc_id==0 )  write(*,*) "scan_angle sum example: ", scan_angle(my_sensor_index,search_lon_min,search_lat_min)

    zenith_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    zenith_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + abs(obs%zenith_angle(iob))

    !if(my_proc_id==0 )  write(*,*) "about to read obs%azimuth_angle "
    azimuth_angle_theta = 90 - obs%azimuth_angle(iob)
    !if(my_proc_id==0 )  write(*,*) "azimuth_angle_theta: ", azimuth_angle_theta
    if (azimuth_angle_theta < 0) azimuth_angle_theta = azimuth_angle_theta + 360
    azimuth_angle_theta = azimuth_angle_theta * RADS_PER_DEGREE

    azimuth_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    azimuth_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + azimuth_angle_theta

    azimuth_sine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    azimuth_sine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + sin(azimuth_angle_theta)

    azimuth_cosine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    azimuth_cosine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + cos(azimuth_angle_theta)

    angle_sum_count(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    angle_sum_count(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + 1

    if (.NOT. write_file_override) then
      lon_radiance_min(my_sensor_index) = min(lon_radiance_min(my_sensor_index), nint(obs%position(iob,1)) - ceiling(1+(search_radius / dx)))
      lon_radiance_max(my_sensor_index) = max(lon_radiance_max(my_sensor_index), nint(obs%position(iob,1)) + ceiling(1+(search_radius / dx)))
      lat_radiance_min(my_sensor_index) = min(lat_radiance_min(my_sensor_index), nint(obs%position(iob,2)) - ceiling(1+(search_radius / dy)))
      lat_radiance_max(my_sensor_index) = max(lat_radiance_max(my_sensor_index), nint(obs%position(iob,2)) + ceiling(1+(search_radius / dy)))
    endif
  enddo check_cycle2
  if (my_proc_id==0 .AND. inputfile == 'fort.80011') write(*,*) "Sensor_IDs_used: ", Sensor_IDs_used

  if (write_file_override) then
    lon_radiance_min = 1
    lon_radiance_max = ix
    lat_radiance_min = 1
    lat_radiance_max = jx
  else
    do isensor = 1,n_sensors_used
      lon_radiance_min(isensor) = max(1,lon_radiance_min(isensor))
      lon_radiance_max(isensor) = min(ix,lon_radiance_max(isensor))
      lat_radiance_min(isensor) = max(1,lat_radiance_min(isensor))
      lat_radiance_max(isensor) = min(jx,lat_radiance_max(isensor))
    enddo
  endif
  !if(my_proc_id==0)  write(*,*) "iob_rad range: ", iob_radmin, iob_radmax
  !if(my_proc_id==0)  write(*,*) "lon_radiance: ", lon_radiance_min, lon_radiance_max
  !if(my_proc_id==0)  write(*,*) "lat_radiance: ", lat_radiance_min, lat_radiance_max

  allocate(numx_crtm(n_sensors_used), STAT = Allocate_Status)
  allocate(numy_crtm(n_sensors_used), STAT = Allocate_Status)
  numx_crtm = lon_radiance_max - lon_radiance_min + 1
  if(my_proc_id==0)  write(*,*) "numx_crtm: ", numx_crtm
  numy_crtm = lat_radiance_max - lat_radiance_min + 1
  if(my_proc_id==0)  write(*,*) "numy_crtm: ", numy_crtm
  
  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. Determine which if any channels that 
  !     brightness temperatures have already
  !     been calculated, and load those
  ! ------------------------------------------

  n_Channels_all = count(is_channel_used)
  !if(my_proc_id==0)  write(*,*) "is_channel_used: ", is_channel_used
  !if(my_proc_id==0)  WRITE(*,*) 'there will be ', n_Channels_all, ' assimilated in total'
  ALLOCATE(n_channels_all_eachSensor(n_sensors_used), STAT = Allocate_status)
  ALLOCATE(channel_numbers_all_eachSensor(n_sensors_used,max_n_channels_all_sensor), STAT = Allocate_status)
  n_channels_all_eachSensor = 0

  ! BTs of all cahnnels from all sensors is sotred in a single array.
  ! allocate array for entire WRF domain, but use only locations relevant
  ! to the observations from each sensor.
  allocate(Tb(ix,jx,n_Channels_all))
  Tb = -8888.

  ! get array of channels used and the sensor they are with, and
  ! load brightness temperature file if it exists. if no file exists,
  ! then identify the channel for CRTM simulation
  CALL GET_ENVIRONMENT_VARIABLE('CRTM_OUT_DIR',CRTM_OUT_DIR)
  !if(my_proc_id==0)  WRITE(*,*) 'CRTM_OUT_DIR: ', trim(CRTM_OUT_DIR)
  ich2 = 0
  crtm_out_files = .FALSE.
  
  do isensor = 1,n_sensors_used
    ich3 = 0
    do ich = 1,N_CHANNELS(sensors_used(isensor))
      !if(my_proc_id==0)  write(*,*) "considering channel ", ich, " of sensor ", sensors_used(isensor)
      if (is_channel_used(sensor_ch_offset(isensor) + ich)) then
        ich2 = ich2 + 1
        !if(my_proc_id==0)  write(*,*) "  this channel is used"

        ich3 = ich3 + 1
        channel_numbers_all_eachSensor(isensor,ich3) = ich

        !if(my_proc_id==0)  write(*,*) 'FILE_NAME should have sensor ', isensor, ' name ', trim(SENSOR_CLASS_NAMES(sensors_used(isensor)))
        write(FILE_NAME,'(a,a,a,a,i2.2,a,a)') trim(CRTM_OUT_DIR), '/Tb_', trim(SENSOR_CLASS_NAMES(sensors_used(isensor))), '_ch', ich,'_', trim(inputfile)
        !if(my_proc_id==0)  WRITE(*,*) '  searching for file ', trim(FILE_NAME)

        INQUIRE(FILE=trim(FILE_NAME),EXIST=crtm_out_files)
        if (crtm_out_files) then ! read file of brightness temperatures
          if(my_proc_id==0)  WRITE(*,*) '  Loading BT data from this file: ', trim(FILE_NAME)
          if (numx_crtm(isensor) .NE. ix .OR. numy_crtm(isensor) .NE. jx) then
            lon_radiance_min(isensor) = 1
            lon_radiance_max(isensor) = ix
            lat_radiance_min(isensor) = 1
            lat_radiance_max(isensor) = jx
            numx_crtm(isensor) = ix
            numy_crtm(isensor) = jx

            if(my_proc_id==0)  write(*,*) "  new numx_crtm(isensor) and numy_crtm(isensor): ", numx_crtm(isensor), numy_crtm(isensor)
          endif

          !if(my_proc_id==0)  WRITE(*,*) '  reading ', numx_crtm*numy_crtm, ' BTs'
          open(ens,file=trim(FILE_NAME),access='direct',recl=4)
          irec = 0
          do j = 1, numy_crtm(isensor)
          !if(my_proc_id==0)  write(*,*) "reading column ", j
          do i = 1, numx_crtm(isensor)
            !if(my_proc_id==0 .AND. j==numy_crtm)  write(*,*) "reading row ", i
            irec= irec +1
            read( ens, rec=irec) Tb(i,j,ich3)
            !if(my_proc_id==0 .AND. i==numx_crtm)  write(*,*) "  read Tb: ", Tb(i,j,ich2)
          enddo
          enddo
          close (ens)
        else
          is_channel_crtm(sensor_ch_offset(isensor) + ich) = .TRUE.
          !if(my_proc_id==0)  write(*,*) "during searching, this many crtm channels: ", count(is_channel_crtm)
        endif
      endif
      !if(my_proc_id==0)  write(*,*) "ich = ", ich, ", if ich ", ich2, " equals n_Channels_all", n_Channels_all, " then we exit"
    enddo
    n_channels_all_eachSensor(isensor) = ich3
    if (ich2 .eq. n_Channels_all) exit  ! stop when all the channels have been found
  enddo
  deallocate(is_channel_used)
  !do isensor = 1,n_sensors_used
  !  if(my_proc_id==0)  WRITE(*,*) 'channel_numbers_all sensor ', isensor, ' : ', channel_numbers_all_eachSensor(isensor,:)
  !enddo
  

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------

  !if(my_proc_id==0)  write(*,*) "done loading any Tb files"
  !if(my_proc_id==0)  write(*,*) "is_channel_crtm: ", is_channel_crtm
  !if(my_proc_id==0)  write(*,*) "about to assign this to n_Channels_crtm: ", count(is_channel_crtm)
  n_Channels_crtm = count(is_channel_crtm)

  ! determine sensor indices and channel numbers for crtm
  ALLOCATE(channel_numbers_crtm(n_Channels_crtm), STAT = Allocate_Status)
  ALLOCATE(sensor_index_crtm(n_Channels_crtm), STAT = Allocate_Status)
  ALLOCATE(channel_numbers_crtm_eachSensor(n_sensors_used,n_Channels_crtm), STAT = Allocate_status)
  ALLOCATE(n_channels_crtm_eachSensor(n_sensors_used), STAT = Allocate_status)
  channel_numbers_crtm_eachSensor = 0
  ich2 = 0
  do isensor = 1,n_sensors_used
    ich3 = 0
    do ich = 1,N_CHANNELS(sensors_used(isensor))
      !if(my_proc_id==0)  write(*,*) "is channel ", ich, " of sensor ", sensors_used(isensor), " for crtm?"
      if (is_channel_crtm(sensor_ch_offset(isensor) + ich)) then
        !if(my_proc_id==0)  write(*,*) "yes"
        ich2 = ich2 + 1
        channel_numbers_crtm(ich2) = ich
        sensor_index_crtm(ich2) = isensor

        ich3 = ich3 + 1
        channel_numbers_crtm_eachSensor(isensor,ich3) = ich
      endif
    enddo
    n_channels_crtm_eachSensor(isensor) = ich3
    if (ich2 .eq. n_Channels_crtm) exit  ! stop when all the channels have been found
  enddo
  deallocate(is_channel_crtm)
  !if(my_proc_id==0)  WRITE(*,*) 'channel_numbers_crtm: ', channel_numbers_crtm
  !if(my_proc_id==0)  WRITE(*,*) 'sensor_index_crtm: ', sensor_index_crtm

  ! ============================================================================

  ! ============================================================================
  ! 3. **** GET INPUT DATA ****
  !
  ! For fill the Atm structure array.
  !
  ! 4a1. Loading Atmosphere and Surface input
  ! --------------------------------
  call get_variable3d(inputfile,'P',ix,jx,kx,1,p)
  call get_variable3d(inputfile,'PB',ix,jx,kx,1,pb)
  call get_variable3d(inputfile,'PH',ix,jx,kx+1,1,ph)
  call get_variable3d(inputfile,'PHB',ix,jx,kx+1,1,phb)
  call get_variable3d(inputfile,'T',ix,jx,kx,1,t)
  call get_variable3d(inputfile,'QVAPOR',ix,jx,kx,1,qvapor)
  if(my_proc_id==0)  WRITE(*,*) 'use cloud: ', use_cloud
  if (use_cloud) then
     call get_variable3d(inputfile,'QCLOUD',ix,jx,kx,1,qcloud)
     call get_variable3d(inputfile,'QRAIN',ix,jx,kx,1,qrain)
     call get_variable3d(inputfile,'QICE',ix,jx,kx,1,qice)
     call get_variable3d(inputfile,'QSNOW',ix,jx,kx,1,qsnow)
     call get_variable3d(inputfile,'QGRAUP',ix,jx,kx,1,qgraup)
  else
     qcloud = 0.
     qrain = 0.
     nrain = 0.
     qice = 0.
     qsnow = 0.
     qgraup = 0.
  endif
  call get_variable2d(inputfile,'PSFC',ix,jx,1,psfc)
  call get_variable2d(inputfile,'TSK',ix,jx,1,tsk)
  call get_variable2d(inputfile,'XLAND',ix,jx,1,xland)
  call get_variable2d(inputfile,'HGT',ix,jx,1,hgt)
  call get_variable2d(inputfile,'U10',ix,jx,1,u10)
  call get_variable2d(inputfile,'V10',ix,jx,1,v10)

  lat = xlat/180.0*MY_PI
  lon = xlong/180.0*MY_PI
  pres = P + PB
  tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )
  where(qvapor.lt.0.0) qvapor=1.0e-8
  where(qcloud.lt.0.0) qcloud=0.0
  where(qice.lt.0.0) qice=0.0
  where(qrain.lt.0.0) qrain=0.0
  where(nrain.lt.0.0) nrain=0.0
  where(qsnow.lt.0.0) qsnow=0.0
  where(qgraup.lt.0.0) qgraup=0.0
  windspeed = sqrt(U10**2 + V10**2)
  where (ISNAN(windspeed)) windspeed = 0.0d0
  westwind = (U10 .lt. 0)
  ! corrected for CRTM version 2.3.0
  ! winddir = -180*(westwind) + ( 90 - ( atan(V10/U10)/RADS_PER_DEGREE) )
  winddir = 180*(1 - westwind) + (90 - ( atan(V10/U10)/RADS_PER_DEGREE) )
  where (ISNAN(winddir)) winddir = 0.0d0



  ! ============================================================================

  ! ============================================================================
  ! 4. **** INITIALIZE CRTM AND ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 4a. This initializes the CRTM for the sensors
  !     predefined in the example SENSOR_ID parameter.
  !     NOTE: The coefficient data file path is hard-
  !           wired for this example.
  ! --------------------------------------------------
  !if(my_proc_id==0) WRITE( *,'(/5x,"Initializing the CRTM...")' )

  CALL CRTM_Version( Version )

  do isensor = 1,n_sensors_used

    ! determine if CRTM gets run for any channels of this sensor
    n_channels_crtm_mySensor = n_channels_crtm_eachSensor(isensor)
    if (n_channels_crtm_mySensor > 0) then  

    if(my_proc_id==0) then
      call cpu_time(cputime)
      WRITE(*,*) 'starting CRTM for sensor ', isensor,' at ', cputime
      WRITE(*,*) 'sensor ', isensor,' is ', trim(Sensor_IDs_used(isensor))
    endif

    ! If sensor class is either gpm_gmi_lf or gmi_gpm_hf (low- or high-frequency),
    ! then initialize the CRTM with 'gpm_gmi'. Otherwise, use the actual Sensor_ID
    ! specified by the observation.
    ! 'gpm_gmi_lf' is channels 1-9, while 'gpm_gmi_hg' is channels 10-13.
    ! These two sets of channels have different scan and zenith angles. If observations
    ! from both of these two sets of channels are assimilated, then treating these sets of
    ! channels as different classes of sensors is necessary to have the calculated average 
    ! scan and zenith angles of observations in the vicinity of a given model grid point 
    ! to be correct for each observation. However, the CRTM still needs to be told 'gmi_gpm'
    ! for either.
    if(INDEX(trim(Sensor_IDs_used(isensor)), 'gmi')>0) then
      CRTM_Init_Sensor_ID = 'gmi_gpm'
    else
      CRTM_Init_Sensor_ID = trim(Sensor_IDs_used(isensor))
    endif
  
    if(my_proc_id==0)  write(*,*) 
    if(my_proc_id==0)  write(*,*) "CRTM ver.",TRIM(Version) 
    if (my_proc_id==0 .AND. (inputfile == 'fort.80011' .OR. inputfile == 'fort.80071')) then
      Error_Status = CRTM_Init( (/CRTM_Init_Sensor_ID/), &  ! Input... must be an array, hencethe (/../)
                                ChannelInfo  , &  ! Output
                                IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                                IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                                File_Path='coefficients/', &
                                CloudCoeff_File_rain  = 'WSM6_RainLUT_-109z-1.bin',&
                                CloudCoeff_File_snow  = 'WSM6_SnowLUT_-109z-1.bin',&
                                CloudCoeff_File_graup = 'WSM6_GraupelLUT_-109z-1.bin',&
                                Quiet=.false.)
    else
      Error_Status = CRTM_Init( (/CRTM_Init_Sensor_ID/), &  ! Input... must be an array, hencethe (/../)
                                ChannelInfo  , &  ! Output
                                IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                                IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                                File_Path='coefficients/', &
                                CloudCoeff_File_rain  = 'WSM6_RainLUT_-109z-1.bin',&
                                CloudCoeff_File_snow  = 'WSM6_SnowLUT_-109z-1.bin',&
                                CloudCoeff_File_graup = 'WSM6_GraupelLUT_-109z-1.bin',&
                                Quiet=.true.)
    end if

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error initializing CRTM'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF

    !if(my_proc_id==0)  write(*,*) "n_Channels_crtm: ", n_Channels_crtm
    !if(my_proc_id==0)  write(*,*) "channel_numbers_crtm: ", channel_numbers_crtm
    !if(my_proc_id==0)  write(*,*) "ChannelInfo n_Channels after CRTM_Init: "
    !if(my_proc_id==0)  call CRTM_ChannelInfo_Inspect(ChannelInfo(1))

    n_channels_crtm_mySensor = n_channels_crtm_eachSensor(isensor)
    ALLOCATE(channel_numbers_crtm_mySensor(n_channels_crtm_mySensor), STAT=Allocate_Status )
    channel_numbers_crtm_mySensor = channel_numbers_crtm_eachSensor(isensor,1:n_channels_crtm_mySensor)
    !if(my_proc_id==0)  write(*,*) "n_channels_crtm_mySensor: ", n_channels_crtm_mySensor
    !if(my_proc_id==0)  write(*,*) "channel_numbers_crtm_mySensor: ", channel_numbers_crtm_eachSensor(isensor,1:n_channels_crtm_mySensor)
    
    !if (my_proc_id==0) write(*,*) 'about to call CRTM_ChannelInfo_Subset'
    Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset = channel_numbers_crtm_mySensor )
    !if (my_proc_id==0) write(*,*) 'successful with CRTM_ChannelInfo_Subset'
    !if(my_proc_id==0)  write(*,*) "ChannelInfo_Subset went fine: ", Error_Status
    !if(my_proc_id==0)  call CRTM_ChannelInfo_Inspect(ChannelInfo(1))

    ! 4b. Allocate the ARRAYS
    ! -----------------------
    ! Note that only those structure arrays with a channel
    ! dimension are allocated here because we've parameterized
    ! the number of profiles in the N_PROFILES parameter.
    !
    ! Users can make the 
    ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
    ALLOCATE( RTSolution( n_channels_crtm_mySensor, N_PROFILES ), STAT=Allocate_Status )
    IF ( Allocate_Status /= 0 ) THEN
      Message = 'Error allocating structure arrays'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
    !if(my_proc_id==0)  write(*,*) "RTSolution allocated"

    ! 4c. Allocate the STRUCTURES
    ! ---------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( Atm, kx, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
      Message = 'Error allocating CRTM Atmosphere structures'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
    !if(my_proc_id==0)  write(*,*) "CRTM_Atmosphere created"

    !if(my_proc_id==0)  write(*,*) "Allocating Tb arrays" 
    !if(my_proc_id==0)  write(*,*) "sizes: ", numx_crtm(isensor), numy_crtm(isensor), n_channels_crtm_mySensor
  
    allocate(Tbsend(numx_crtm(isensor),numy_crtm(isensor),n_channels_crtm_mySensor))
    allocate(Tb_crtm(numx_crtm(isensor),numy_crtm(isensor),n_channels_crtm_mySensor))
    Tbsend = 0.
    Tb_crtm = 0.

    ! 4a2. Parallerization with grids
    ! --------------------------------

    ! original system divids domain into 1-pixel strips in the y-dimension
    ! new system divdes domain by single pixels
    !ystart=my_proc_id+1
    !yend=numy_crtm
    !do j=ystart, yend, nprocs
    !   y = (j-1) + lat_radiance_min
    !   !WRITE(*,*) 'proc_id ', my_proc_id, ' j=',j,', y=',y
    !do i=1, numx_crtm
    !   x = (i-1) + lon_radiance_min
    grand_count = 0
    do j=1, numy_crtm(isensor)
       y = (j-1) + lat_radiance_min(isensor)
    do i=1, numx_crtm(isensor)
       x = (i-1) + lon_radiance_min(isensor)
       grand_count = grand_count + 1
       !if(my_proc_id==0)  write(*,'(a,i4,i4,i6)') 'point ', x, y, grand_count
       if (mod(grand_count,nprocs) .EQ. my_proc_id) then ! calcualte if this processor should
       !if (my_proc_id==0) then ! calcualte if this processor should
       !if(my_proc_id==0)  write(*,'(a,i4,i4,i6)') 'proc 0 is calculating', x, y, grand_count
       !write(*,'(a,i3,a,i4,i4,i6)') 'proc ', my_proc_id,' is calculating', x, y, grand_count

    ! 4a3. Converting WRF data for CRTM structure
    ! --------------------------------
    !--- converting the data to CRTM structure

    !*******************************************************************************
    ! satellite information is contained in parameter arrays
    !*******************************************************************************


    ! 4a. GeometryInfo input
    ! ----------------------


    ! if this location was not designated for use by any observation, then use the average value of
    ! nearby locations. increase search "radius" until averaging over at least 2 locations.
    my_scan_angle = scan_angle(isensor,x,y)
    !if(my_proc_id==0)  write(*,*) "my_scan_angle before search: ", my_scan_angle 
    scan_angle_search_radius = 0
    scan_angle_search_number = 0
    if (my_scan_angle < 2*tiny(my_scan_angle)) then  

      !if(my_proc_id==0 .AND. inputfile == 'fort.80011')  write(*,*) "    point", x, y, " needs a scan angle"
      do while (scan_angle_search_number < 2)
        scan_angle_search_radius = scan_angle_search_radius + 1
        scan_angle_search_lonmin = max(1, x-scan_angle_search_radius)
        scan_angle_search_lonmax = min(ix,x+scan_angle_search_radius)
        scan_angle_search_latmin = max(1, y-scan_angle_search_radius)
        scan_angle_search_latmax = min(jx,y+scan_angle_search_radius)
        !if(my_proc_id==0)  write(*,*) "  max angle_sum_count: ", maxval(angle_sum_count(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
        !                                                                        scan_angle_search_latmin:scan_angle_search_latmax))
        scan_angle_search_number = count(angle_sum_count(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                                 scan_angle_search_latmin:scan_angle_search_latmax) .GT. 1)
        !if(my_proc_id==0)  write(*,*) "  scan_angle_search_number: ", scan_angle_search_number
      enddo
      
      my_scan_angle = sum(scan_angle(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                 scan_angle_search_latmin:scan_angle_search_latmax) )
      my_zenith_angle = sum(zenith_angle(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                     scan_angle_search_latmin:scan_angle_search_latmax) )
      my_azimuth_sine = sum(azimuth_sine(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                     scan_angle_search_latmin:scan_angle_search_latmax) )
      my_azimuth_cosine = sum(azimuth_cosine(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                         scan_angle_search_latmin:scan_angle_search_latmax) )
      my_azimuth_angle_est = sum(azimuth_angle(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                           scan_angle_search_latmin:scan_angle_search_latmax) )
      my_angle_sum_count = sum(angle_sum_count(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                           scan_angle_search_latmin:scan_angle_search_latmax) )
      !if(my_proc_id==0)  write(*,*) "  my_scan_angle: ", my_scan_angle
      !if(my_proc_id==0)  write(*,*) "  final scan_angle_search_radius: ", scan_angle_search_radius
    else
      !if(my_proc_id==0)  write(*,*) "  no need to search for locations with scan angle here"
      my_zenith_angle = zenith_angle(isensor,x,y)
      my_azimuth_sine = azimuth_sine(isensor,x,y)
      my_azimuth_cosine = azimuth_cosine(isensor,x,y)
      my_azimuth_angle_est = azimuth_angle(isensor,x,y)
      my_angle_sum_count = angle_sum_count(isensor,x,y)
      !if(my_proc_id==0)  write(*,*) "  my_angle_sum_count: ", my_angle_sum_count
    endif

    my_scan_angle = my_scan_angle / my_angle_sum_count
    my_zenith_angle = my_zenith_angle / my_angle_sum_count

    my_azimuth_angle = atan2(my_azimuth_sine / my_angle_sum_count, my_azimuth_cosine / my_angle_sum_count) / RADS_PER_DEGREE
    my_azimuth_angle = 90 - my_azimuth_angle
    if (my_azimuth_angle < 0) my_azimuth_angle = my_azimuth_angle + 360
 
    my_azimuth_angle_est = (my_azimuth_angle_est / my_angle_sum_count) / RADS_PER_DEGREE
    my_azimuth_angle_est = 90 - my_azimuth_angle_est
    if (my_azimuth_angle_est < 0) my_azimuth_angle_est = my_azimuth_angle_est + 360

    !if(my_proc_id==0)  write(*,*) "  final my_scan_angle: ", my_scan_angle
    !if(my_proc_id==0 .AND. inputfile == 'fort.80071' .AND. (my_scan_angle < 45 .OR. my_scan_angle > 60) )  write(*,*) "    final my_scan_angle: ", my_scan_angle
    !if(my_proc_id==0 .AND. inputfile == 'fort.80011' .AND. my_azimuth_angle > (MY_PI/-1.25) .AND. my_azimuth_angle < (MY_PI/1.25) .AND. j .EQ. 100 )  write(*,*) "    i", i, ", my_azimuth_angle and my_azimuth_angle_est: ", my_azimuth_angle, my_azimuth_angle_est
    if(my_proc_id==0 .AND. inputfile == 'fort.80011' .AND. j .EQ. 1 )  write(*,*) "    i", i, ", my_azimuth_angle and my_azimuth_angle_est: ", my_azimuth_angle, my_azimuth_angle_est
    !if(my_proc_id==0 .AND. inputfile == 'fort.80071' .AND. (my_zenith_angle < 40 .OR. my_zenith_angle > 55) )  write(*,*) "    final my_zenith_angle: ", my_zenith_angle
    if ( my_scan_angle > 90 ) then
      write(*,*) "x",x," y",y, ", adjusted my_scan_angle: ", my_scan_angle
      write(*,*) "  final scan_angle_search_radius: ", scan_angle_search_radius
      write(*,*) "  my_angle_sum_count: ", my_angle_sum_count
      write(*,*) "  adjusted my_zenith_angle: ", my_zenith_angle
      write(*,*) "  scan_angle_search_lon: ", scan_angle_search_lonmin, scan_angle_search_lonmax
      write(*,*) "  scan_angle_search_lat: ", scan_angle_search_latmin, scan_angle_search_latmax
      write(*,*) "  my_proc_id: ", my_proc_id
      write(*,*)
    endif
    !if(my_proc_id==0)  write(*,*) ""

    YEAR_STR = times(1:4)
    MONTH_STR = times(6:7)
    DAY_STR = times(9:10)
    read(YEAR_STR,*) YEAR
    read(MONTH_STR,*) MONTH
    read(DAY_STR,*) DAY

    CALL CRTM_Geometry_SetValue( Geometry, &
                                 Sensor_Zenith_Angle  = my_zenith_angle, &
                                 Sensor_Scan_Angle    = my_scan_angle, &
                                 Sensor_Azimuth_Angle = my_azimuth_angle, &
                                 Year                 = YEAR, &
                                 Month                = MONTH, &
                                 Day                  = DAY )
    ! 4b. Converting WRF data for CRTM structure
    ! --------------------------------

    if ( use_slant_path) then
      call Load_CRTM_Structures_MW_slant( Atm, Sfc, x, y, kx, &
               my_zenith_angle, my_azimuth_angle, dx, .FALSE.)
    else
      call Load_CRTM_Structures_MW( Atm, Sfc, x, y, kx, .FALSE.)
    end if


    ! 4c. Use the SOI radiative transfer algorithm
    ! --------------------------------------------
    Options%RT_Algorithm_ID = RT_SOI


    ! 4d. Specify number of streams
    IF (N_STREAMS .GT. 0) THEN
      Options%Use_N_Streams = .TRUE.
      Options%n_Streams = N_STREAMS
    END IF
    ! ============================================================================

    ! ============================================================================
    ! 5. **** CALL THE CRTM FORWARD MODEL ****
    !

    if (Sfc(1)%Water_Coverage < 0) then
      WRITE(*,*)
      WRITE(*,*) 'x: ', x
      WRITE(*,*) 'y: ', y
      WRITE(*,*) 'xland(x,y): ', xland(x,y)
      WRITE(*,*) 'Water coverage: ', Sfc(1)%Water_Coverage
      WRITE(*,*) 'Land coverage: ', Sfc(1)%Land_Coverage
      WRITE(*,*) 'Pressure z=10: ', atm(1)%Pressure(10)
      WRITE(*,*)
    endif

    Error_Status = CRTM_Forward( Atm        , &
                                 Sfc        , &
                                 Geometry   , &
                                 ChannelInfo, &
                                 RTSolution , &
                                 Options = Options )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Forward Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
    ! ============================================================================



    ! ============================================================================
    ! 6. **** Collecting output ****
    !
    ! User should read the user guide or the source code of the routine
    ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
    ! select the needed variables for outputs.  These variables are contained
    ! in the structure RTSolution.

    !---for file output, edited 2014.9.26
    do l = 1, n_channels_crtm_mySensor
        Tbsend(i,j,l) = real(RTSolution(l,1)%Brightness_Temperature)

        !if(my_proc_id==0 .AND. mod(grand_count,nprocs*300)==my_proc_id)  write(*,'(a,i4,i4,a,i4,a,f7.2)') 'proc 0: ', x, y, ', ch ', l, ' Tb = ',Tbsend(i,j,l)
        !if(my_proc_id==0)  write(*,'(a,i4,i4,a,f7.2)') 'proc 0: ', x, y, ', Tb = ',Tbsend(i,j,l)
        if ((i .EQ. 121 .AND. j .EQ. 1)) then
          WRITE(*,*) '  at x=121, y=1, Tbsend=',Tbsend(i,j,l)
        endif
        if (Tbsend(i,j,l) .NE. Tbsend(i,j,l)) then
          write(*,*) '  Tbsend is NaN at x=',x,' y=',y
          write(*,*) 'Pressures: ', atm(1)%Pressure
          write(*,*) 'Level_pressures: ', atm(1)%Level_Pressure

        endif
        if (Tbsend(i,j,l) .GT. 999) then
          WRITE(*,*) '  at x=',i,'y=',j,'Tbsend=',Tbsend(i,j,l)
          write(*,*) 'Pressures: ', atm(1)%Pressure
          write(*,*) 'Level_pressures: ', atm(1)%Level_Pressure

        endif
    enddo
    !WRITE(*,'(7x,"Profile (",i0,", ",i0,") finished Tb =  ",f6.2)')x,y,Tbsend(x,y,2)

    !--- end of x and y loops
    endif ! calcualte if this processor should
    enddo ! loop over x dimension
    enddo ! loop over y dimension

    CALL MPI_Allreduce(Tbsend,Tb_crtm,numx_crtm(isensor)*numy_crtm(isensor)*n_channels_crtm_mySensor,MPI_REAL,MPI_SUM,comm,ierr)
    !if(my_proc_id==0)  write(*,*) "out of spatial loop"


    ! ============================================================================

    ! ============================================================================
    !6.5  **** satellite beam convolution, writing the output ****
    !

    if(my_proc_id==0) then

    ! write crtm BTs to file
    if (.false.) then
    !if (numx_crtm(isensor) .EQ. ix .AND. numy_crtm(isensor) .EQ. jx) then
      do l = 1, n_channels_crtm_mySensor

        write(*,*) 'FILE_NAME should have sensor ', isensor, ' name ', trim(SENSOR_CLASS_NAMES(sensors_used(isensor)))
        write(FILE_NAME,'(a,a,a,a,i2.2,a,a)') trim(CRTM_OUT_DIR), '/Tb_', trim(SENSOR_CLASS_NAMES(sensors_used(isensor))), &
                                       '_ch', channel_numbers_crtm_mySensor(l), '_', trim(inputfile)
        write(*,*) 'saving TBs to ', trim(FILE_NAME)
        open(ens,file=trim(FILE_NAME),access='direct',recl=4)
        write(*,*) 'saving TBs file open'
        irec = 0
        do j = 1, numy_crtm(isensor)
        do i = 1, numx_crtm(isensor)
          irec= irec +1
          write( ens, rec=irec) Tb_crtm(i,j,l)
          !if(i == 1) write(*,'(i4,i4,a,f7.2)') i, j, ', Tb = ',Tb_crtm(i,j,l)
        enddo
        enddo
        !deallocate(FILE_NAME)
      enddo
       close (ens)
    endif  ! proc_0 writing simulated Tbs to file
    if(my_proc_id==0)  write(*,*) "done writing CRTM simulations to file"
    endif

    endif  ! whether the CRTM gets run for this sensor at all

    !write(*,*) "proc_id", my_proc_id, " is done running CRTM for sensor ", isensor
    n_channels_all_mySensor = n_channels_all_eachSensor(isensor)
    ALLOCATE(channel_numbers_all_mySensor(n_channels_all_mySensor), STAT=Allocate_Status )
    channel_numbers_all_mySensor = channel_numbers_all_eachSensor(isensor,1:n_channels_all_mySensor)
    !write(*,*) "proc_id", my_proc_id, " is done figuring out channel_numbers_all_mySensor for sensor ", isensor

    ! merge file-read and crtm BTs into Tb array (BTs of all channels of all sensors) before looping through all obs for convolution
    ich2 = 1
    mySensor_Tb_offset = 0
    do l = 1, isensor-1
      mySensor_Tb_offset = mySensor_Tb_offset + n_channels_all_eachSensor(l)
      !write(*,*) "proc_id", my_proc_id, " sees that mySensor Tb_offset is  ", mySensor_Tb_offset
    enddo
    do l = 1, n_Channels_all_mySensor
      !write(*,*) "proc_id", my_proc_id, " is considering my sensor channel index ", l
      !write(*,*) "proc_id", my_proc_id, " is considering my sensor channel number ", channel_numbers_all_mySensor(l)
      !write(*,*) "proc_id", my_proc_id, " is considering my sensor CRTM channel index ", ich
      !write(*,*) "proc_id", my_proc_id, " asks if CRTM channel number ", channel_numbers_crtm_mySensor(ich2), " is equal to other channel number"
      if (channel_numbers_all_mySensor(l) .EQ. channel_numbers_crtm_mySensor(ich2) ) then
        !write(*,*) "proc_id", my_proc_id, " found that yes it does"
        !write(*,*) "proc_id", my_proc_id, " asks what l+mySensor_Tb_offset is:", l+mySensor_Tb_offset
        Tb(:,:,l+mySensor_Tb_offset) = Tb_crtm(:,:,ich2)
        !if(my_proc_id==0 )  WRITE(*,*) '  at x=121, y=1, Tb=',Tb(121,1,l+mySensor_Tb_offset)
        ich2 = ich2 + 1
      endif
      if (ich2 .GT. n_channels_crtm_mySensor) exit
    enddo
    !if(my_proc_id==0)  write(*,*) "done merging file and CRTM BTs"

    ! debug
    !if(my_proc_id==0) then
    !  WRITE(*,'(a10,"   Tb=",f6.2,"~",f6.2)')inputfile,minval(Tb),maxval(Tb)
    !  write(*,*) 'all the Tb values for this sensor:'
    !  do l = 1, n_channels_crtm_mySensor
    !    do j = 1, numy_crtm(isensor)
    !    write(*,*) 'row ',j
    !    do i = 1, numx_crtm(isensor), 3
    !      write(*,'(3f7.2)') Tb(i:(i+2),j,l)
    !    enddo
    !    enddo
    !  enddo
    !endif


    if(my_proc_id==0) then
      call cpu_time(cputime)
      WRITE(*,*) 'proc0 starting satellite beam convultion ', cputime
    end if
    if(my_proc_id==1) then
      call cpu_time(cputime)
      WRITE(*,*) 'proc1 starting satellite beam convultion ', cputime
    end if

    if(my_proc_id==0)  write(*,*) iob_radmin, nprocs, iob_radmax, isensor
    !if(my_proc_id==0)  write(*,*) 'length of obs%ype: ', size(obs%type)
  
    ALLOCATE(channel_search(n_Channels_all_mySensor), STAT = Allocate_Status)  ! array used to identify the channel number index from the channel number specified by the obs, for beam convolution
    do iob=iob_radmin+my_proc_id,iob_radmax,nprocs   ! parallelize across all procs by observation
      obstype = obs%type(iob)
      !if(iob < 30 .AND. inputfile == 'fort.80071')  write(*,*) "proc ", my_proc_id, ", iob ", iob, ", obstype ", obstype(1:9), ", sensor_index ", sensor_indices(iob-iob_radmin+1)
      if ( obstype(1:9) .EQ. 'Microwave' .AND. sensor_indices(iob-iob_radmin+1) .EQ. isensor ) then

        obs_ii=obs%position(iob,1)-lon_radiance_min(isensor)+1
        obs_jj=obs%position(iob,2)-lat_radiance_min(isensor)+1
        !if(my_proc_id==0)  write(*,*) "iob ", iob, " obs%position ", obs%position(iob,1), ", obs_ii ", obs_ii
        !if(my_proc_id==0)  write(*,*) "iob ", iob, " obs%position ", obs%position(iob,2), ", obs_jj ", obs_jj
    
        channel_search = abs(obs%ch(iob) - channel_numbers_all_mySensor)
        ich = minloc(channel_search,1) + mySensor_Tb_offset
        !if(my_proc_id==0)  write(*,*) "iob ", iob, " ich ", ich

        !write(*,*) "for iob ", iob, ", ", obs_ii, obs_jj, obs%efov_aScan(iob), obs%efov_cScan(iob)

        !xb_tb(iob) = beam_conv_simple(numx_crtm(isensor), numy_crtm(isensor), Tb(:,:,ich), &
        !                 obs_ii, obs_jj, dx, dy, obs%efov_aScan(iob), obs%efov_cScan(iob) )
        xb_tb(iob) = beam_conv_gaussian_simple(numx_crtm(isensor), numy_crtm(isensor), Tb(:,:,ich), &
                         obs_ii, obs_jj, dx, dy, obs%efov_aScan(iob), obs%efov_cScan(iob) )
        !if(my_proc_id==0)  write(*,*) "iob ", iob, " obs%efov_aScan ", obs%efov_aScan(iob), ", obs%efov_cScan ", obs%efov_cScan(iob)
        !if(iob < 30 .AND. inputfile == 'fort.80071')  write(*,*) "iob ", iob, " numx_crtm ", numx_crtm, ", numy_crtm ", numy_crtm
        !write(*,*) "for iob ", iob, ", xb_tb = ", xb_tb(iob)

        if (xb_tb(iob) > 500 .OR. xb_tb(iob) .LE. 1 .OR. (xb_tb(iob) .NE. xb_tb(iob)) ) then
          WRITE(*,*) 'iob ', iob, ',  obs_ii ', obs_ii,      ',  obs_jj ', obs_jj
          WRITE(*,*) '        xb_tb(iob): ', xb_tb(iob)
        endif
      
      !elseif (obstype(1:9) .NE. 'Microwave') then
      !  WRITE(*,*) 'this iob did not have type Microwave: ', iob
      endif
    enddo
    deallocate(channel_search)

    ! write summary of results
    if(my_proc_id==0) then
      WRITE(*,*) lon_radiance_min(isensor), lon_radiance_max(isensor), lat_radiance_min(isensor), lat_radiance_max(isensor)
      !WRITE(*,*) 'all xb_tb on proc0: ', xb_tb(iob_radmin:iob_radmax)
      WRITE(*,'(a10,"   Tb=",f6.2,"~",f6.2)')inputfile,minval(Tb),maxval(Tb)
      WRITE(*,'(" Tb_conv=",f6.2,"~",f6.2)')minval(xb_tb(iob_radmin:iob_radmax),DIM=1,MASK=xb_tb(iob_radmin:iob_radmax) .GT. 0),maxval(xb_tb(iob_radmin:iob_radmax))
    endif
    do i = 10, nprocs, 10
      if(my_proc_id==i) then
        WRITE(*,'("proc_",i3.3," Tb_conv=",f6.2,"~",f6.2)')my_proc_id,minval(xb_tb(iob_radmin:iob_radmax),DIM=1,MASK=xb_tb(iob_radmin:iob_radmax) .GT. 0),maxval(xb_tb(iob_radmin:iob_radmax))
      endif
    enddo

    ! ============================================================================
    !  **** initializing all Tb and Tbsend fields ****
    !
    deallocate(channel_numbers_all_mySensor)

    if (n_Channels_crtm_mySensor .gt. 0) then
      Tbsend = 0.0
      Tb_crtm = 0.0
      CALL MPI_BCAST(Tbsend,numx_crtm(isensor)*numy_crtm(isensor)*max(n_channels_crtm_mySensor,1),MPI_REAL,0,comm,ierr)

      ! ============================================================================
      ! 7. **** DESTROY THE CRTM ****
      !
      deallocate(Tb_crtm)
      deallocate(Tbsend)
      deallocate(channel_numbers_crtm_mySensor)
      deallocate(RTSolution)
  
      Error_Status = CRTM_Destroy( ChannelInfo )
      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error destroying CRTM'
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP
      END IF
    end if
    ! ============================================================================

  end do ! cycle to the next sensor

  if(my_proc_id==0) then
    call cpu_time(cputime)
    WRITE(*,*) 'finished with xb.f at ', cputime
  endif

  deallocate(Tb)

  deallocate(sensors_used)

  deallocate(lon_radiance_min)
  deallocate(lon_radiance_max)
  deallocate(lat_radiance_min)
  deallocate(lat_radiance_max)

  deallocate(scan_angle)
  deallocate(zenith_angle)
  deallocate(angle_sum_count)
  deallocate(sensor_indices)

  deallocate(Sensor_IDs_used)

  deallocate(numx_crtm)
  deallocate(numy_crtm)

  deallocate(n_channels_all_eachSensor)
  deallocate(channel_numbers_all_eachSensor)
  deallocate(channel_numbers_crtm)
  deallocate(sensor_index_crtm)
  deallocate(channel_numbers_crtm_eachSensor)
  deallocate(n_channels_crtm_eachSensor)


CONTAINS

  INCLUDE 'Load_CRTM_Structures_MW.inc'
  INCLUDE 'Load_CRTM_Structures_MW_slant.inc'

end subroutine xb_to_microwave



