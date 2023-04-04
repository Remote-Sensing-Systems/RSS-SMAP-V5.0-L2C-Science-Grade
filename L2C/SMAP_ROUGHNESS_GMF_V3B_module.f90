module SMAP_ROUGHNESS_GMF_V3B_module
   use external_files_l2c_module
   use dir_paths_module
   use l2c_module_smap

   implicit none
   private
   save


   integer(4), parameter                       :: n_rad=3
   real(4), parameter                          :: sst_ref0=20.0, sss_ref0=35.0, freq_aq=1.413, wcasp=11.0, wextrapol=17.0, teff=290.
   real(4), parameter, dimension(n_rad)        :: tht_ref = (/29.36, 38.44, 46.29/)

! reduce kappascale from 1.4 (Aquarius V5) to 0.7 (SMAP V3/V4/V5)
   real(4), parameter                          :: kappascale=0.7 ! empirical scale factor for delta form factor

! input files
! read in from tausq_file
   integer(2)                                  :: itau(360,180,12,91)
   integer(4), parameter                       :: iu=3


   public :: find_demiss_rough_wspd_SMAP_V3, find_demiss_rough_wdir_SMAP_V3
   public :: find_demiss_wind_SMAP_RESIDUAL
   public :: find_demiss_rough_sstformfactor_SMAP, find_refl_RTM, fd_water_refl
   public :: freq_aq, teff

contains


   subroutine find_demiss_rough_wspd_SMAP_V3 (wspd, sst, tht,    demiss_rough)
! isotropic part of roughness correction
! Aquarius V5.0
! interpolate EIA
      implicit none

      real(4), intent(in)                                 ::  wspd ! [m/s]
      real(4), intent(in)                                 ::  sst  ! [Celsius]
      real(4), intent(in)                                 ::  tht  ! [deg]


      real(4), intent(out), dimension(2)                  ::  demiss_rough ! [*290K]

      real(4), dimension(2)                               ::  dew2, dew3
      integer(4)                                          ::  irad
      real(4)                                             ::  brief
      real(4), dimension(2)                               ::  yy

      brief = (tht-tht_ref(2))/(tht_ref(3)-tht_ref(2)) ! linear interpolation to Aquarius EIA

      irad=2
      call find_demiss_rough_wspd_AQV5 (2, wspd, sst,    dew2)

      irad=3
      call find_demiss_rough_wspd_AQV5 (3, wspd, sst,    dew3)

      yy = dew2*(1.0-brief) + dew3*brief

      demiss_rough = yy

      return
   end subroutine find_demiss_rough_wspd_SMAP_V3



   subroutine find_demiss_wind_SMAP_RESIDUAL(wspd, demiss_wind_res)
      use external_files_l2c_module
! posthoc correction to SMAP isotropic wind emissivity
! derived in C:/SMAP/SSS_algo_V3/roughness_model/residual_DTB_W_T_PHI_WH
! and C:/SMAP/SSS_algo_V3/roughness_model/IDL/analyze_SMAP_dew_residual_3.pro
      implicit none

      real(4), intent(in)                     ::  wspd
      real(4), intent(out), dimension(2)      ::  demiss_wind_res ! [*290K]

      real(4)                     ::  xwspd, brief, yy(2), ww0 ,ww1
      integer(4), parameter       ::  n=500
      real(4), parameter          ::  b0=0.0, del=0.1
      real(8), dimension(2), save ::  dw_res_tab(2,0:n)
      integer(4), save            ::  istart=1
      integer(4)                  ::  ipol, i0, i1

      call get_dir_paths(table_dir, data_dir)
      if (istart == 1) then
         istart=0
         open(unit=3,file=trim(table_dir) // '/' // demiss_res_file,action='read',form='unformatted',access='stream',status='old')
         read(3) dw_res_tab
      endif

      xwspd = wspd
      if (wspd>14.5) xwspd=14.5 ! cut off at 14.5 m/s
      if (wspd<0.0001) xwspd=0.0001

      i0 = floor((xwspd-b0)/del)
      i1= i0+1

      if (i0<0 .or. i0>n .or. i1<0 .or. i1>n) then
         write(*,*) wspd,xwspd,i0,i1, ' error in find_demiss_wind_SMAP_RESIDUAL. pgm stop.'
         stop
      endif

      ww0= b0 + float(i0)*del
      ww1= b0 + float(i1)*del

      brief = (xwspd-ww0)/(ww1-ww0)

      do ipol=1,2
         yy(ipol) = (1.0-brief)*dw_res_tab(ipol,i0) + brief*dw_res_tab(ipol,i1)
      enddo

      demiss_wind_res = yy

      return
   end subroutine find_demiss_wind_SMAP_RESIDUAL



! SMAP V3. Also used for V4 and V5
   subroutine find_demiss_rough_wdir_SMAP_V3(wspd, phir,   demiss_rough_wdir)
      implicit none

      real(4), intent(in)                                 ::  wspd ! [m/s]
      real(4), intent(in)                                 ::  phir ! [deg] mathematical anti-clockwise


      real(4), intent(out), dimension(4)                  ::  demiss_rough_wdir ! [*290K] 1=V, 2=H, 3=S3, 4=S4

      integer(4)                                          ::  ipol

      real(4), dimension(0:2,4)                           ::  aharm  !1=V, 2=H, 3=S3, 4=S4
! the 0th harmonic is set ot 0

      if (phir < -900.) then
         demiss_rough_wdir = 0.0 ! isotropic
         return
      endif

      call fd_SMAP_harmonics(wspd,        aharm)

      do ipol=1,2 ! odd
         demiss_rough_wdir(ipol) = aharm(1,ipol)*cos((phir)*rad) + aharm(2,ipol)*cos((2.0*phir)*rad)
      enddo

      do ipol=3,4 ! even
         demiss_rough_wdir(ipol) = aharm(1,ipol)*sin((phir)*rad) + aharm(2,ipol)*sin((2.0*phir)*rad)
      enddo

      return
   end subroutine find_demiss_rough_wdir_SMAP_V3



   subroutine fd_SMAP_harmonics(wspd,      aharm)
      use external_files_l2c_module
      implicit none

      real(4), intent(in)                                 :: wspd

      real(4), dimension(0:2,4), intent(out)              :: aharm  !1=V, 2=H, 3=S3, 4=S4
                                                                    ! the 0th harmonic is set ot 0

      integer(4), parameter                               :: npoly=5


      real(4), dimension(4)                               :: A0,  A1,  A2  !1=V, 2=H, 3=S3, 4=S4


      real(4)                                             :: ww
      integer(4)                                          :: ipol, iharm
      real(4)                                             :: fval

      real(4), parameter                                  :: wcut=24.5 ! cutoff for 1st and 2nd harm


      integer(4), save                                    :: istart=1
      real(8), dimension(0:2,4,npoly), save               :: acoef
      ! harmonic coefficients for radiometer wind direction signal
      ! V/H/S3/S4
      call get_dir_paths(table_dir, data_dir)
      if (istart==1) then
         istart=0
         open(unit=3,file=trim(table_dir) // '/' // emiss_coeff_harm_file,action='read',&
         &form='unformatted',access='stream',status='old')
         read(3) acoef
         close(3)
      endif

      A0=0.0
      A1=0.0
      A2=0.0

      do ipol=1,4

         ! A0 set to 0
         ! The A0 is taken form the Aquarius V5 GMF
         A0(ipol) = 0.0

         ! A1
         iharm=1
         ww = wspd
         if (wspd >= wcut) ww=wcut ! cutoff at wcut
         fval = ww*acoef(iharm,ipol,1) + (ww**2) *acoef(iharm,ipol,2) +       (ww**3)*acoef(iharm,ipol,3) + &
         &(ww**4)*acoef(iharm,ipol,4) +   (ww**5)*acoef(iharm,ipol,5)
         A1(ipol)  = fval

         ! A2
         iharm=2
         ww = wspd
         if (wspd >= wcut) ww=wcut ! cutoff at wcut
         fval = ww*acoef(iharm,ipol,1) + (ww**2) *acoef(iharm,ipol,2) +       (ww**3)*acoef(iharm,ipol,3) +  &
         &(ww**4)*acoef(iharm,ipol,4) +   (ww**5)*acoef(iharm,ipol,5)
         A2(ipol)  = fval

      enddo !ipol

      do ipol=1,4
         aharm(0,ipol) = A0(ipol)
         aharm(1,ipol) = A1(ipol)
         aharm(2,ipol) = A2(ipol)
      enddo

      return
   end subroutine fd_SMAP_harmonics



! Aquarius V5.0
   subroutine fd_TM_emiss_harmonics(irad,wspd,     aharm, daharm)
      use external_files_l2c_module
      implicit none

      integer(4), intent(in)                              :: irad
      real(4), intent(in)                                 :: wspd

      real(4), dimension(0:2,2), intent(out)              ::  aharm  !1=V, 2=H
      real(4), dimension(0:2,2), intent(out), optional    :: daharm  !1=V, 2=H

      integer(4), parameter                               :: n_rad=3, npoly=5


      real(4), dimension(2)                               ::  A0,  A1,  A2  !1=V, 2=H
      real(4), dimension(2)                               :: dA0, dA1, dA2  !1=V, 2=H


      real(4)                                             :: ww
      integer(4)                                          :: ipol, iharm
      real(4)                                             :: fval, dval

      real(4)                                             :: w0, w1, w2 ! linear extrapolation/cutoff point

      integer(4), save                                    :: istart=1
      real(8), dimension(0:2,2,n_rad,npoly), save         :: acoef      ! harmonic coefficients for radiometer wind direction signal
      real(8), dimension(0:2,2,n_rad), save               :: wspd_max_a ! high wind speed for radiometer wind speed signal


      call get_dir_paths(table_dir, data_dir)
      if (istart==1) then
         istart=0
         open(unit=3,file=trim(table_dir) // '/' // emiss_coeff_harm_file_AQ, form='unformatted',&
         & access='stream', action='read', status='old')
         read(3) acoef
         read(3) wspd_max_a
         ! overwrite
         wspd_max_a = wextrapol
         close(3)
      endif


      A0=0.0
      A1=0.0
      A2=0.0

      do ipol=1,2

         ! A0
         iharm=0
         w0 = wspd_max_a(iharm,ipol,irad)
         ww = wspd
         if (wspd >= w0) ww=w0 ! extrapolation at w0
         fval = &
            ww*acoef(iharm,ipol,irad,1)  +  (ww**2)*acoef(iharm,ipol,irad,2) +       (ww**3)*acoef(iharm,ipol,irad,3) +&
         &(ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)
         dval = &
            acoef(iharm,ipol,irad,1)  + (2.0*ww)*acoef(iharm,ipol,irad,2) + (3.0*(ww**2))*acoef(iharm,ipol,irad,3) +&
         &(4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5)

         if (wspd<=w0) then
            A0(ipol) = fval
         else
            A0(ipol) = fval + dval*(wspd-w0)
         endif

         dA0(ipol) = dval


         ! A1
         iharm=1
         w1 = wspd_max_a(iharm,ipol,irad)
         ww = wspd
         if (wspd >= w1) ww=w1 ! cutoff at w1
         fval = ww*acoef(iharm,ipol,irad,1) + (ww**2) *acoef(iharm,ipol,irad,2) +       (ww**3)*acoef(iharm,ipol,irad,3)      +&
         &(ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)
         dval =    acoef(iharm,ipol,irad,1) + (2.0*ww)*acoef(iharm,ipol,irad,2) + (3.0*(ww**2))*acoef(iharm,ipol,irad,3)      +&
         &(4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5)

         A1(ipol)  = fval
         dA1(ipol)  = dval

         ! A2
         iharm=2
         w2 = wspd_max_a(iharm,ipol,irad)
         ww = wspd
         if (wspd >= w2) ww=w2 ! cutoff at w2
         fval = ww*acoef(iharm,ipol,irad,1) + (ww**2)      *acoef(iharm,ipol,irad,2) +        (ww**3)*acoef(iharm,ipol,irad,3) +&
         &(ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)
         dval =    acoef(iharm,ipol,irad,1) + (2.0*ww)     *acoef(iharm,ipol,irad,2) +  (3.0*(ww**2))*acoef(iharm,ipol,irad,3) +&
         &(4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5)

         A2(ipol)  = fval
         dA2(ipol) = dval

      enddo !ipol


      do ipol=1,2
         aharm(0,ipol) = A0(ipol)
         aharm(1,ipol) = A1(ipol)
         aharm(2,ipol) = A2(ipol)
      enddo

      if (present(daharm)) then
         do ipol=1,2
            daharm(0,ipol) = dA0(ipol)
            daharm(1,ipol) = dA1(ipol)
            daharm(2,ipol) = dA2(ipol)
         enddo
      endif


      return
   end subroutine fd_TM_emiss_harmonics



   subroutine find_demiss_rough_wspd_AQV5 (irad, wspd, sst,    demiss_rough)
      use external_files_l2c_module
! isotropic part of roughness correction
      implicit none

      integer(4), intent(in)                              ::  irad
      real(4), intent(in)                                 ::  wspd ! [m/s]
      real(4), intent(in)                                 ::  sst  ! Celsius


      real(4), intent(out), dimension(2)                  ::  demiss_rough ! [*290K]

      real(4)                                             ::  xsst, xwspd, xtht

      real(4), dimension(2)                               ::  dew_1, dew_2, xem0, yem0

      real(4), dimension(0:2,2)                           ::  aharm  !1=V, 2=H
      real(4), dimension(0:2,2)                           :: daharm  !1=V, 2=H

      real(4), dimension(2)                               ::  delta
      integer(4)                                          ::  ipol

      integer(4)                                          ::  ksst1, ksst2
      real(4)                                             ::  x1, x2, brief, y1, y2

      integer(4), parameter                               ::  msst=35, n_rad=3
      integer(4), save                                    ::  istart=1


      real(4), dimension(2,msst,n_rad), save              ::  dtab
      real(4), parameter                                  ::  sst_step=1.0, sst0=0.0, sstmax=30.0, wmax=wcasp

      call get_dir_paths(table_dir, data_dir)
      if (istart==1) then
         istart=0
         open(unit=3,form='unformatted',access='stream',file=trim(table_dir) // '/' // delta_file,action='read',status='old')
         read(3) dtab
         close(3)
      endif

      call fd_TM_emiss_harmonics(irad,wspd,           aharm, daharm)
      dew_1(1:2) = aharm(0,1:2) ! *290K

      xtht=tht_ref(irad)
      call fdem0_meissner_wentz(freq_aq,xtht,sst,    sss_ref0, xem0)
      call fdem0_meissner_wentz(freq_aq,xtht,sst_ref0,sss_ref0, yem0)

! sst adjustment
      xsst=sst
      if (xsst<sst0+sst_step/2)   xsst=sst0+sst_step/2
      if (xsst>sstmax) xsst=sstmax

      xwspd=wspd
      if (xwspd<0.0)  xwspd=0.0
      if (xwspd>wmax) xwspd=wmax

      call fd_TM_emiss_harmonics(irad,xwspd,          aharm, daharm)
      dew_2(1:2) = aharm(0,1:2)   ! *290K
! *290K dew_2 = dew_1 below wmax. above wmax it is kept constant


      ksst1 = floor((xsst - (sst0+sst_step/2))/sst_step) + 1
      if (ksst1< 1)           ksst1 = 1
      if (ksst1> msst-1)      ksst1 = msst-1
      ksst2 = ksst1 + 1

      do ipol=1,2
         x1 = sst0 + (ksst1-1)*sst_step + sst_step/2
         x2 = x1   + sst_step
         brief = (xsst-x1)/sst_step
         y1 = dtab(ipol,ksst1,irad)
         y2 = dtab(ipol,ksst2,irad)
         delta(ipol) = y1*(1.0-brief) + y2*brief
      enddo ! ipol

! total roughness correction
      do ipol=1,2
         demiss_rough(ipol) = dew_1(ipol)*(xem0(ipol)/yem0(ipol)) + kappascale*delta(ipol)*dew_2(ipol)
      enddo

      return
   end subroutine find_demiss_rough_wspd_AQV5


   subroutine find_demiss_rough_sstformfactor_SMAP (tht, sst,    delta1, delta2)
      use external_files_l2c_module
! SST form factors in roughness correction
      implicit none

      real(4), intent(in)                                 ::  tht
      real(4), intent(in)                                 ::  sst  ! Celsius
      real(4), intent(out), dimension(2)                  ::  delta1, delta2

      real(4)                                             ::  xsst, brief_t

      real(4), dimension(2)                               ::  xem0, yem0, xdelta, yy, ydelta2, ydelta3

      integer(4)                                          ::  ipol, irad

      integer(4)                                          ::  ksst1, ksst2
      real(4)                                             ::  x1, x2, brief, y1, y2

      integer(4), parameter                               ::  msst=35, n_rad=3
      integer(4), save                                    ::  istart=1


      real(4), dimension(2,msst,n_rad), save              ::  dtab
      real(4), parameter                                  ::  sst_step=1.0, sst0=0.0, sstmax=30.0

      call get_dir_paths(table_dir, data_dir)
      if (istart==1) then
         istart=0
         open(unit=3,form='unformatted',access='stream',file=trim(table_dir) // '/' // delta_file,action='read',status='old')
         read(3) dtab
         close(3)
      endif

      call fdem0_meissner_wentz(freq_aq,tht,sst,    sss_ref0, xem0)
      call fdem0_meissner_wentz(freq_aq,tht,sst_ref0,sss_ref0, yem0)

      delta1 = xem0/yem0

      xsst=sst
      if (xsst<sst0+sst_step/2)   xsst=sst0+sst_step/2
      if (xsst>sstmax) xsst=sstmax

      ksst1 = floor((xsst - (sst0+sst_step/2))/sst_step) + 1
      if (ksst1< 1)           ksst1 = 1
      if (ksst1> msst-1)      ksst1 = msst-1
      ksst2 = ksst1 + 1

      irad=2
      do ipol=1,2
         x1 = sst0 + (ksst1-1)*sst_step + sst_step/2
         x2 = x1   + sst_step
         brief = (xsst-x1)/sst_step
         y1 = dtab(ipol,ksst1,irad)
         y2 = dtab(ipol,ksst2,irad)
         xdelta(ipol) = y1*(1.0-brief) + y2*brief
         ydelta2(ipol) =  kappascale*xdelta(ipol)
      enddo ! ipol


      irad=3
      do ipol=1,2
         x1 = sst0 + (ksst1-1)*sst_step + sst_step/2
         x2 = x1   + sst_step
         brief = (xsst-x1)/sst_step
         y1 = dtab(ipol,ksst1,irad)
         y2 = dtab(ipol,ksst2,irad)
         xdelta(ipol) = y1*(1.0-brief) + y2*brief
         ydelta3(ipol) =  kappascale*xdelta(ipol)
      enddo ! ipol


! EIA interpolation
      brief_t = (tht-tht_ref(2))/(tht_ref(3)-tht_ref(2)) ! linear interpolation to Aquarius EIA
      yy = ydelta2*(1.0-brief_t) + ydelta3*brief_t

      delta2 = yy


      return
   end subroutine find_demiss_rough_sstformfactor_SMAP


   subroutine find_refl_RTM(idayjl,xlat,xlon,tht,phir,frac_land,frac_ice,sss,surtep,wind,lst,sm, refl_tot,tbsur)
! this is Frank's routine find_refl_tot
! adds swh as optional parameter
      implicit none

      integer(4), intent( in) :: idayjl
      real(4),    intent( in) :: xlat,xlon,tht,phir,frac_land,frac_ice,sss,surtep,wind,lst,sm

      real(4), dimension(2), intent(out), optional :: refl_tot, tbsur

      real(4)     ::  sst,frac_ice_adj,frac_water
      real(4)     ::  refl_water(2),refl_ice(2),refl_land(2),tausq,sm_adj, xrefl_tot(2), xtbsur(2)


      sst=surtep - 273.15
      if(sst.lt.-2.) sst=-2.
      if(sst.gt.35.) sst=35.

      frac_ice_adj=frac_ice
      if(frac_ice.gt.1.-frac_land) frac_ice_adj=1.-frac_land

      frac_water=1. - frac_land - frac_ice_adj

      refl_water=0.
      refl_ice  =0.
      refl_land =0.

      if(frac_water  >= 1.0E-7)   call fd_water_refl(tht=tht,sss=sss,sst=sst,wind=wind,phir=phir, refl=refl_water)

      if(frac_ice_adj>= 1.0E-7)   call fd_seaice_refl(tht,   refl_ice)

      if(frac_land   >=1.0E-7)  then

         ! using 1.01 rather than 1 accounts for possible roundoff error
         sm_adj=sm
         if(sm_adj.gt.1.01) sm_adj=0.35 !fill in missing values, looking at the sm images 0.35 is a typical value for islands
         if(sm_adj.gt.1.)   sm_adj=1.   !fixes any roundoff error

         if(lst.ge.273.15) then
            call fd_land_refl(1,tht,sm_adj, refl_land)  !1 denotes soil
         else
            call fd_land_refl(2,tht,sm_adj, refl_land)  !2 denotes snow
         endif

         call fd_tausq(idayjl,xlat,xlon,tht, tausq)

         refl_land = refl_land*0.74*tausq  !0.74=exp(-.3) 0.3 is rough hght,b param as in pellarin et al. 2003
         if(minval(refl_land).lt.0 .or. maxval(refl_land).gt.1) stop 'land refl oob, pgm stopped'

      endif ! calculate land reflectivity

      xrefl_tot=frac_water*refl_water + frac_ice_adj*refl_ice + frac_land*refl_land
      xtbsur=frac_water*(1-refl_water)*(sst+273.15) + frac_ice_adj*(1-refl_ice)*lst + frac_land*(1-refl_land)*lst

      if (present(refl_tot))  refl_tot=xrefl_tot
      if (present(tbsur))     tbsur=xtbsur

      return
   end subroutine find_refl_RTM




   subroutine fd_tausq(idayjl,xlat,xlon,tht, tausq)
      implicit none


      integer(4), intent(in)  :: idayjl
      real(4)   , intent(in)  :: xlat,xlon,tht
      real(4)   , intent(out) :: tausq

      integer(4),save ::  istart = 1
      integer(4)      ::  i1,i2,j1,j2,k1,k2,l1,l2
      real(4)         ::  a1,a2,b1,b2,c1,c2,d1,d2,brief,xmon


      call get_dir_paths(table_dir, data_dir)
      if(istart.eq.1) then
         istart=0
         open(iu,file=trim(table_dir) // '/' // tausq_file,access='stream',form='unformatted',status='old',action='read')
         read(iu) itau
         close(iu)
         if(maxval(itau).ne.10000) stop 'itau table was not read in, pgm stopped'
      endif

      if(idayjl.lt.1 .or. idayjl.gt.366.)  then
         write(*,*) idayjl
         stop 'idayjl oob in fd_vege_tran, pgm stopped'
      endif

      if(abs(xlat).gt.90.)                 then
         write(*,*) xlat
         stop 'xlat   oob in fd_vege_tran, pgm stopped'
      endif

      if(xlon.lt.0..or.xlon.gt.360.)       then
         write(*,*) xlon
         stop 'xlon   oob in fd_vege_tran, pgm stopped'
      endif
      if(tht.lt.0. .or.tht.gt.  90.)       then
         write(*,*) tht
         stop 'tht    oob in fd_vege_tran, pgm stopped'
      endif

      xmon=12.*(idayjl - 0.5)/365.25
      if(xmon.gt.11.9999) xmon=11.9999

      brief=xmon-0.5
      i1=1+brief
      i2=i1+1
      a1=i1-brief
      a2=1-a1
      if(i1.eq. 0) i1=12
      if(i2.eq.13) i2= 1

      brief=xlat+89.5
      j1=1+brief
      j2=j1+1
      b1=j1-brief
      b2=1.-b1
      if(j1.eq.  0) j1=  1
      if(j2.eq.181) j2=180

      brief=xlon-0.5
      k1=1+brief
      k2=k1+1
      c1=k1-brief
      c2=1-c1
      if(k1.eq.  0) k1=360
      if(k2.eq.361) k2=  1

      brief=tht
      if(brief.gt.89.999) brief=89.999
      l1=1+brief
      l2=l1+1
      d1=l1-brief
      d2=1-d1

      tausq=                                                          &
         d1*(a1*b1*(c1*itau(k1,j1,i1,l1)+c2*itau(k2,j1,i1,l1))+          &
         a1*b2*(c1*itau(k1,j2,i1,l1)+c2*itau(k2,j2,i1,l1))+          &
         a2*b1*(c1*itau(k1,j1,i2,l1)+c2*itau(k2,j1,i2,l1))+          &
         a2*b2*(c1*itau(k1,j2,i2,l1)+c2*itau(k2,j2,i2,l1)))+         &
         d2*(a1*b1*(c1*itau(k1,j1,i1,l2)+c2*itau(k2,j1,i1,l2))+      &
         a1*b2*(c1*itau(k1,j2,i1,l2)+c2*itau(k2,j2,i1,l2))+          &
         a2*b1*(c1*itau(k1,j1,i2,l2)+c2*itau(k2,j1,i2,l2))+          &
         a2*b2*(c1*itau(k1,j2,i2,l2)+c2*itau(k2,j2,i2,l2)))

      tausq=1.e-4*tausq

      return
   end subroutine fd_tausq



   subroutine fd_land_refl(itype,tht,sm, refl)
      implicit none

      integer(4), intent(in)  ::   itype  !1=spol above freezing    2=snow    3=sea ice
      real(4), intent(in)     ::   tht,sm

      real(4), dimension(2),intent(out)    ::   refl

      real(4)             ::  costht, sinsqtht
      complex(4)          ::  rv,rh, esqrt, permit, permit_alpha

      real(4), parameter  ::  acoef=1.850842
      !from fd_soil_coefs.f assuming freq= 1.413

      complex(4), parameter ::   bcoef=(16.06828,-0.7877967)
      ! from fd_soil_coefs.f assuming freq= 1.413

      real(4), parameter :: alpha =0.65
      real(4), parameter :: beta  =1.09

      if(itype.lt.1 .or. itype.gt.3) stop 'itype oob in fd_land_refl, pgm stopped'

      !   this is the ulaby,moore,fumg model, vol 3 page 2103 that sab used.
      if(itype.eq.1) then  !soil above freezing
         permit_alpha=acoef + bcoef*sm**beta
         permit = permit_alpha**(1./alpha)  !soil
      endif

      !    page 2026 and 2027 give the following for fresh-water ice at 1.4 ghz and t=-20:   permit_ice=(3.15,-0.0005)
      !    fig e.26 on page 2062 gives real part of snow to be about 1.73 for snow density of 0.4
      !    fig e.28 on page 2065 gives imag part of snow to be 0.23*(imag part of ice)=0.23*0.0005  for snow density of 0.4

      if(itype.eq.2)  permit=(1.73,-0.00012)  !snow

      if(itype.eq.3)  permit=(3.,-0.1)        !sea ice

      costht=cos((tht)*rad)
      sinsqtht=1.-costht*costht

      esqrt=csqrt(permit-sinsqtht)
      rh=(costht-esqrt)/(costht+esqrt)
      rv=(permit*costht-esqrt)/(permit*costht+esqrt)
      refl(1) =rv*conjg(rv)
      refl(2) =rh*conjg(rh)

      return
   end subroutine fd_land_refl


   subroutine fd_seaice_refl(tht, refl)
      implicit none

      real(4), intent(in)                 ::  tht
      real(4), dimension(2), intent(out)  :: refl

      real(4)     ::  costht, sinsqtht
      complex(4)  ::  rv,rh, esqrt, permit

      permit=(3.,-0.1)        !sea ice

      costht=cos((tht)*rad)
      sinsqtht=1.-costht*costht

      esqrt=csqrt(permit-sinsqtht)
      rh=(costht-esqrt)/(costht+esqrt)
      rv=(permit*costht-esqrt)/(permit*costht+esqrt)
      refl(1) =rv*conjg(rv)
      refl(2) =rh*conjg(rh)

      return
   end subroutine fd_seaice_refl



   subroutine fd_water_refl(tht,sss,sst,wind,phir, refl)
      implicit none

      real(4),  intent(in)                :: tht,sss,sst,wind
      real(4), intent(in), optional       :: phir
      real(4),  dimension(2), intent(out) :: refl


      real(4),dimension(2)  ::  dems,em0, dems1, demiss_w_res
      real(4),dimension(4)  ::  dems_phi


      call fdem0_meissner_wentz(freq_aq,tht,sst,sss, em0)
      refl=1-em0

      !     add wind induced reflective (minus emissivity)
      if (present(phir)) then
         call find_demiss_rough_wspd_SMAP_V3 (wind, sst, tht,   dems1)
         call find_demiss_rough_wdir_SMAP_V3 (wind, phir,    dems_phi)
         dems(1:2) = dems1(1:2) + dems_phi(1:2)
      else
         call find_demiss_rough_wspd_SMAP_V3 (wind, sst, tht,    dems)
      endif

      ! add small posthoc correction
      call find_demiss_wind_SMAP_RESIDUAL(wind, demiss_w_res)
      dems(1:2) = dems(1:2) + demiss_w_res(1:2)

      dems = dems/teff

      refl=refl - dems

      return
   end subroutine fd_water_refl




end module SMAP_ROUGHNESS_GMF_V3B_module


