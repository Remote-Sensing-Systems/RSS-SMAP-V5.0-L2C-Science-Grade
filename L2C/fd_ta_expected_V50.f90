! V5.0 processing code
! use HYCOM as reference salinity

! TM 01/28/2022
! SSS REF -> TB SUr0 -> TB SUR -> TB TOA -> TB TOI -> TA EARTH -> TA   exepcted
! apply sea-ice and land correction

!    No SSS retreival AND no land correction is performed if gland>0.1 or fland>0.1.
!    The TA expected calculation uses fland if no land correction is performed.
!    If a land correction is performed it applies the land correction and sets fland=0.
!    No SSS retreival AND no SIC correction is performed if icezone=5 or (icezone=6 and iceflag_amsr2(2)=1).
!    The TA expected calculation uses fice if no SIC is performed.
!    If a SIC is performed it applies the SIC and sets fice=0.

!
! iopt parameters
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref

subroutine fd_ta_expected
   use l2c_module_smap
   use SMAP_ROUGHNESS_GMF_V3B_module
   implicit none


   integer(4)              ::  idir, ilon, ilat
   real(4), dimension(4)   ::  tbsur0,  tbsur, tbtoa, tbtoi, taearth, ta
   real(4), dimension(2)   ::  refl0,refl,xbuf
   real(4)                 ::  wind, thtadj, tbdown, sst, phir, prot, xlat, xlon

   real(8)                 ::  time_2000,secyr,secdy
   integer(4)              ::  lyear,idayjl,imon,idaymo,isecdy

   real(4)                 ::  xfland, xfice


! initialize as missing
   tb_sur0_exp   = missing_val4
   tb_sur_exp    = missing_val4
   tb_toa_exp    = missing_val4
   tb_toi_exp    = missing_val4
   tb_toi_exp_2  = missing_val4
   ta_earth_exp  = missing_val4
   ta_ant_exp    = missing_val4

   fra_exp       = missing_val4
   pratot_exp    = missing_val4

   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            time_2000 = time(idir,ilon,ilat)
            call fd_date_2000(time_2000, secyr,lyear,idayjl,imon,idaymo,secdy)
            isecdy=nint(secdy)

            xlat=cellat(idir,ilon,ilat)
            xlon=cellon(idir,ilon,ilat)

            sst     = surtep(ilon,ilat) - 273.15
            wind    = winspd(ilon,ilat)
            phir    = eaa(idir,ilon,ilat)-windir(ilon,ilat) !relative wind direction angle, 0=upwind
            thtadj  = eia(idir,ilon,ilat)

            ! new in V5:
            ! avoid double counting land and sea-ice correction
            ! over land and sea-ice NO corrections are performed
            ! fland/fice are used
            ! However, if a correction is perfromed, then fland/fice need to be set to 0 to avoid double counting

            xfland = 0.0
            if (fland(idir,ilon,ilat)>0.1 .or. gland(idir,ilon,ilat)>0.1) xfland=fland(idir,ilon,ilat)

            xfice = 0.0
            if (icezone(ilon,ilat)==5) xfice = 1.0
            if (icezone(ilon,ilat)==6 .and. iceflag_amsr2(2,ilon,ilat)==1) xfice = 1.0


            ! TB sur0 (specular)
            call find_refl_RTM(idayjl=idayjl,xlat=xlat,xlon=xlon,tht=thtadj,&
               phir=phir,frac_land=xfland,frac_ice=xfice,sss=sss_ref(ilon,ilat),surtep=surtep(ilon,ilat),&
               wind=0.0,       lst=lst(ilon,ilat),sm=sm(ilon,ilat), refl_tot=refl0,tbsur=tbsur0(1:2))
            tbsur0(3:4) = 0.0

            ! SIC
            if (icezone(ilon,ilat)>=1 .and. icezone(ilon,ilat)<=4) then
               tbsur0(1:2) = tbsur0(1:2) + dtb_sea_ice_corr(1:2,ilon,ilat)
               ! dtb_sea_ice_corr set to 0 in fd_tb_sur_sic if icezone>=5
            endif

            tb_sur0_exp(1:4,idir,ilon,ilat) = tbsur0(1:4)

            ! TB sur (including roughness)
            ! Note: Value for refl is needed later in the TB TOA computation
            call find_refl_RTM(idayjl=idayjl,xlat=xlat,xlon=xlon,tht=thtadj,&
               phir=phir,frac_land=xfland,frac_ice=xfice,sss=sss_ref(ilon,ilat),surtep=surtep(ilon,ilat),&
               wind=wind,        lst=lst(ilon,ilat),sm=sm(ilon,ilat), refl_tot=refl,tbsur=tbsur(1:2))
            tbsur(3:4) = (1.0-xfland)*(1.0-xfice)*dtb_rough(3:4,idir,ilon,ilat)
            ! ensure to be 0 over land or sea ice

            ! SIC
            if (icezone(ilon,ilat)>=1 .and. icezone(ilon,ilat)<=4) then
               tbsur(1:2) = tbsur(1:2) + dtb_sea_ice_corr(1:2,ilon,ilat)
               ! dtb_sea_ice_corr set to 0 in fd_tb_sur_sic if icezone>=5
            endif

            tb_sur_exp(1:4,idir,ilon,ilat) = tbsur(1:4)

            ! TB TOA
            tbdown=tbdw(ilon,ilat)+ tran(ilon,ilat)*tb_cos
            ! Frank has subtracted tb_cos = 3K when deriving the reflected galactic tables
            ! atmospheric correction: tb_toa -> tb_sur
            tbtoa(1:2) = tbup(ilon,ilat) + tran(ilon,ilat)*(tbsur(1:2) +  refl(1:2)*tbdown)

            tbtoa(3:4) = tbsur(3:4)*(tran(ilon,ilat)**2)

            if (iopt(1)==1 .and. gland(idir,ilon,ilat)>0.0005) then
               tbtoa(1:2) = tbtoa(1:2) + ta_lnd(1:2,idir,ilon,ilat)
               ! ta_lnd set to 0 in fd_tbtoa_lc if fland>0.1 or gland>0.1
            endif
            tb_toa_exp(1:4,idir,ilon,ilat) = tbtoa(1:4)


            ! expected Faraday rotation and total PRA
            fra_exp(idir,ilon,ilat)    = tec(idir,ilon,ilat)*frdrot(idir,ilon,ilat)
            ! Frank had already included the xtec_scal=0.75 factor in the L2B processing fot TEC

            pratot_exp(idir,ilon,ilat) = fra_exp(idir,ilon,ilat) + pra(idir,ilon,ilat)

            ! TB TOI
            ! This uses the SMAP measured S3
            ! assuming that the S3 at the surface and TOA is zero.
            ! This is inaccurate at high wind speeds.
            xbuf(1:2) = tb_toa_exp(1:2,idir,ilon,ilat)
            call vh_to_stokes(xbuf)
            prot = pra_smap(idir,ilon,ilat)
            tbtoa(1:2) = xbuf(1:2)
            tbtoi(1) = tbtoa(1)
            tbtoi(2) = tbtoa(2)*cos((2.0*prot)*rad)
            tbtoi(3) = tbtoa(2)*sin((2.0*prot)*rad)
            tbtoi(4) = tbtoa(4)
            xbuf(1:2) = tbtoi(1:2)
            call stokes_to_vh(xbuf)
            tbtoi(1:2) = xbuf(1:2)
            tb_toi_exp(1:4,idir,ilon,ilat) = tbtoi(1:4)

            ! TB TOI is also calculated using the expected polarization rotation angle and the value is saved
            xbuf(1:2) = tb_toa_exp(1:2,idir,ilon,ilat)
            call vh_to_stokes(xbuf)
            prot = pratot_exp(idir,ilon,ilat)
            tbtoa(1:2) = xbuf(1:2)
            tbtoi(1) = tbtoa(1)
            tbtoi(2) = tbtoa(2)*cos((2.0*prot)*rad)
            tbtoi(3) = tbtoa(2)*sin((2.0*prot)*rad)
            tbtoi(4) = tbtoa(4)
            xbuf(1:2) = tbtoi(1:2)
            call stokes_to_vh(xbuf)
            tbtoi(1:2) = xbuf(1:2)
            tb_toi_exp_2(1:4,idir,ilon,ilat) = tbtoi(1:4)

            ! TA Earth
            ! from here on only the SMAP measured S3 is used
            xbuf(1:2) = tb_toi_exp(1:2,idir,ilon,ilat)
            call vh_to_stokes(xbuf)
            tbtoi(1:2) = xbuf(1:2)
            taearth = matmul(amat_inv_IQ,tbtoi)
            xbuf(1:2) = taearth(1:2)
            call stokes_to_vh(xbuf)
            taearth(1:2) = xbuf(1:2)
            ta_earth_exp(1:4,idir,ilon,ilat) = taearth(1:4)

            ! TA antenna
            xbuf(1:2) = ta_earth_exp(1:2,idir,ilon,ilat)
            call vh_to_stokes(xbuf)
            taearth(1:2) = xbuf(1:2)
            ta(1:3) = taearth(1:3) + &
               iopt(2)*ta_gal_dir(1:3,idir,ilon,ilat) + iopt(3)*ta_gal_ref(1:3,idir,ilon,ilat) + &
               iopt(4)*ta_sun_dir(1:3,idir,ilon,ilat) + iopt(5)*ta_sun_ref(1:3,idir,ilon,ilat)
            ta(4)   = taearth(4)
            xbuf(1:2) = ta(1:2)
            call stokes_to_vh(xbuf)
            ta(1:2) = xbuf(1:2)
            ta_ant_exp(1:4,idir,ilon,ilat) = ta(1:4)


         enddo ! idir
      enddo ! ilon
   enddo ! ilat

   return
end subroutine fd_ta_expected
