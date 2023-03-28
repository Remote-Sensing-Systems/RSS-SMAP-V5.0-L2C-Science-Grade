! TM 01/28/2022
! adapted to V5
! include sea-ice correction
! sea-ice correction and gice_est scaled by 0.72 (28% reduction) from computed value based on Thomas Meissner's findings
! this was done in decode_l2b


! TM 04/25/2015
! TB TOA -> TB SUR -> TB SUR0
! 1. atmospheric correction
! 2. surface roughness correction

!
! iopt parameters
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref

subroutine fd_tb_sur_sic
   use l2c_module_smap
   use SMAP_ROUGHNESS_GMF_V3B_module
   implicit none


   integer(4)              ::  idir, ilon, ilat
   real(4), dimension(4)   ::  tbtoa, tbsur
   real(4)                 ::  wind, thtadj, tbdown, sst, phir
   real(4), dimension(2)   ::  xval, emiss, demiss_rough, demiss_w ! [*290K]
   real(4), dimension(2)   ::  demiss_w_res ! [*290K]
   real(4), dimension(4)   ::  demiss_phi !1=V 2=H 3=S3 4=S4

! initialize as missing
   tb_sur   = missing_val4
   tb_sur0  = missing_val4
   dtb_rough  = missing_val4

   tb_sur0_sic = missing_val4

! roughness correction
   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            sst     = surtep(ilon,ilat) - 273.15

            wind    = winspd(ilon,ilat)

            phir    = eaa(idir,ilon,ilat)-windir(ilon,ilat) !relative wind direction angle, 0=upwind
            phir    = - phir ! TM 11/15/2015 sign change to account for S3/S4 convention

            thtadj  = eia(idir,ilon,ilat)

            tbdown=tbdw(ilon,ilat) + tran(ilon,ilat)*tb_cos
            ! Frank has subtracted tb_cos = 3K when deriving the reflected galactic tables

            ! atmospheric correction: tb_toa -> tb_sur
            tbtoa(1:4) = tb_toa_lc(1:4,idir,ilon,ilat)
            xval(1:2) =(tbtoa(1:2) - tbup(ilon,ilat))/tran(ilon,ilat)
            emiss(1:2)=(xval(1:2) - tbdown)/(surtep(ilon,ilat) -tbdown)
            tbsur(1:2)=emiss(1:2)*surtep(ilon,ilat)
            tbsur(3:4)=tbtoa(3:4)/(tran(ilon,ilat)**2)
            tb_sur(1:4,idir,ilon,ilat) =  tbsur(1:4)

            ! skip roughness correction if more than 50% land or within sea-ice zone 5
            tb_sur0(1:2,idir,ilon,ilat) = tb_sur(1:2,idir,ilon,ilat)
            tb_sur0(3:4,idir,ilon,ilat) = 0.0

            if (gland(idir,ilon,ilat) > 0.5) cycle
            ! changed in V5.0
            if (icezone(ilon,ilat) == 5) cycle
            if (icezone(ilon,ilat) == 6 .and. iceflag_amsr2(2,ilon,ilat) /=0) cycle !AMSR2 land contaminated. Frank's sea-ice flag set.
            if (icezone(ilon,ilat) == 7 .and. iceflag_amsr2(1,ilon,ilat) /=0) cycle !AMSR2 land contaminated. Clim. sea-ice flag set.


            ! surface roughness correction
            call find_demiss_rough_wspd_SMAP_V3 (wind, sst, thtadj,     demiss_w)
            call find_demiss_rough_wdir_SMAP_V3 (wind, phir,          demiss_phi)

            ! small post-hoc residuals
            call find_demiss_wind_SMAP_RESIDUAL(wind, demiss_w_res)

            demiss_rough(1:2) = demiss_w(1:2) + demiss_phi(1:2) + demiss_w_res(1:2)

            dtb_rough(1:2,idir,ilon,ilat) =  demiss_rough(1:2)*surtep(ilon,ilat)/teff
            dtb_rough(3:4,idir,ilon,ilat) =  demiss_phi(3:4)*surtep(ilon,ilat)/teff

            tb_sur0(1:2,idir,ilon,ilat) = tb_sur(1:2,idir,ilon,ilat) - dtb_rough(1:2,idir,ilon,ilat)
            tb_sur0(3:4,idir,ilon,ilat) = 0.0


         enddo ! idir
      enddo ! ilon
   enddo ! ilat

! sea-ice correction
   tb_sur0_sic = tb_sur0

   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            ! set sea-ice correction to 0 in zones 0,5,6,7
            if (icezone(ilon,ilat) == 0) dtb_sea_ice_corr(1:2,ilon,ilat)=0.0
            if (icezone(ilon,ilat) >= 5) dtb_sea_ice_corr(1:2,ilon,ilat)=0.0

            tb_sur0_sic(1:2,idir,ilon,ilat) = tb_sur0(1:2,idir,ilon,ilat) - dtb_sea_ice_corr(1:2,ilon,ilat)

         enddo ! idir
      enddo ! ilon
   enddo ! ilat

   return
end subroutine fd_tb_sur_sic
