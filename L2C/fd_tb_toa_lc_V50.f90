! TM 01/29/2022
! adapted to V5.0
! No SSS retreival AND no land correction is performed if gland>0.1 or fland>0.1.


! TM 05/03/2019
! adapted to V4.0
! new land correction
! add DTB accounting for difference in land model TB and SMAP TB


! TM 04/25/2015
! TA -> TA EARTH -> TB TOI -> TB TOA + land correction
! 1. apply APC
! 2. remove Faraday rotation
! 3. apply land correction

!
! iopt parameters
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref
! 6: swh in roughness correction

subroutine fd_tb_toa_lc
   use l2c_module_smap
   use SMAP_ROUGHNESS_GMF_V3B_module

   implicit none

   integer(4)              ::  idir, ilon, ilat
   real(4), dimension(4)   ::  taearth, tbtoi, tbtoa
   real(4), dimension(2)   ::  xbuf
   real(4)                 ::  phix
   real(4), dimension(4)   ::  dems_phi

   real(4)                 ::  wind, phir, tbtoa3_wind, zval

   real(4)                 ::  dtb_2(2)
   real(8)                 ::  xtime
   real(8)                 ::  secyr,secdy
   integer(4)              ::  lyear,idayjl,imon,idaymo,isecdy
   real(4)                 ::  xhour, xlat, xlon




! initialize as missing
   tb_toi   = missing_val4
   tb_toa   = missing_val4
   tb_toa_lc= missing_val4
   pra_smap = missing_val4


   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            ! APC
            taearth(1:4) = ta_earth(1:4,idir,ilon,ilat)
            xbuf(1:2)  = taearth(1:2)
            call vh_to_stokes(xbuf)
            taearth(1:2) = xbuf(1:2) ! ta earth in IQ basis
            tbtoi = matmul(amat_IQ, taearth)
            xbuf(1:2) = tbtoi(1:2)
            call stokes_to_vh(xbuf)
            tbtoi(1:2) = xbuf(1:2)
            tb_toi(1:4,idir,ilon,ilat) = tbtoi(1:4) ! save array in VH basis

            ! Faraday Rotation Correction
            ! changed 11/15/2015: include estimate for surface 3rd Stokes from wind
            xbuf(1:2) = tbtoi(1:2)
            call vh_to_stokes(xbuf)
            tbtoi(1:2) = xbuf(1:2)

            wind    = winspd(ilon,ilat)
            phir    = eaa(idir,ilon,ilat)-windir(ilon,ilat) !relative wind direction angle, 0=upwind
            phir    = - phir ! TM 11/15/2015 sign change to account for S3/S4 convention

            call find_demiss_rough_wdir_SMAP_V3(wind, phir,   dems_phi) ! wdir emissivity signal at surface [*290K]
            tbtoa3_wind  = dems_phi(3)*(surtep(ilon,ilat)/teff)*(tran(ilon,ilat)**2) ! estimated 3rd Stokes from wind roughness

            ! This was changed on 04/17/2012: atand needs to be substituted by atan2d in order to work for negative Q
            ! which can occur during space maneuvers

            tbtoa(1)=tbtoi(1)

            zval = (tbtoi(2)**2 + tbtoi(3)**2 - tbtoa3_wind**2)
            if (zval<0.0001) zval=0.0001
            tbtoa(2)=sqrt(zval) !faraday rotation correction

            tbtoa(3)=tbtoa3_wind

            tbtoa(4)=tbtoi(4)

            ! this block is unchanged
            if (abs(tbtoi(3))>1.0E-8 .or. abs(tbtoi(2))>1.0E-8) then
               phix = atan2(tbtoi(3),tbtoi(2)) / rad
               pra_smap(idir,ilon,ilat) = phix/2.0 ! measured SMAP pol rot angle assuming that U(TOA)=0
            else ! pol rotation angle ill-defined
               pra_smap(idir,ilon,ilat) =  missing_val4
            endif

            xbuf(1:2) = tbtoa(1:2)
            call stokes_to_vh(xbuf)
            tbtoa(1:2) = xbuf(1:2)

            tb_toa(1:4,idir,ilon,ilat) = tbtoa(1:4)

            ! land correction
            if (iopt(1)==1 .and. gland(idir,ilon,ilat)>0.0005) then

               ! compute step 2 land correction part that accounts for difference between land model TB and land SMAP TB
               xlat = cellat(idir,ilon,ilat)
               xlon = cellon(idir,ilon,ilat)
               if (xlon<0.0001)   xlon=0.0001
               if (xlon>359.9999) xlon=359.9999
               xtime = time(idir,ilon,ilat)
               call fd_date_2000(xtime, secyr,lyear,idayjl,imon,idaymo,secdy)
               isecdy=nint(secdy)
               xhour = secdy/3600.

               call land_corr_step2(xlat, xlon, zang(idir,ilon,ilat), gland(idir,ilon,ilat), imon,  dtb_2)

               ta_lnd(1:2,idir,ilon,ilat) = ta_lnd(1:2,idir,ilon,ilat) + dtb_2(1:2)

               ! set land correction to 0 if fland>0.1 or gland 0.1
               ! we do not retreive SSS in that instance.
               if (fland(idir,ilon,ilat)>0.1) ta_lnd(1:2,idir,ilon,ilat)=0.0
               if (gland(idir,ilon,ilat)>0.1) ta_lnd(1:2,idir,ilon,ilat)=0.0

               tbtoa(1:2) = tbtoa(1:2) - ta_lnd(1:2,idir,ilon,ilat)

            endif ! do land correction

            tb_toa_lc(1:4,idir,ilon,ilat) = tbtoa(1:4)


         enddo ! idir
      enddo ! ilon
   enddo ! ilat



   return
end subroutine fd_tb_toa_lc
