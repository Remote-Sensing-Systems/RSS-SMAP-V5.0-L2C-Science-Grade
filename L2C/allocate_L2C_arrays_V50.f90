subroutine allocate_L2C_arrays
   use l2c_module_smap
   implicit none

   allocate(xobs_sat_ccmp(nlon,nlat))
   allocate(sst_err(nlon,nlat))

   allocate(surtep_save(nlon,nlat))
   allocate(winspd_save(nlon,nlat))
   allocate(windir_save(nlon,nlat))

   allocate(tb_consistency_save(2,nlon,nlat))


   allocate(sss_0(2,nlon,nlat))
   allocate(sss_P(2,nlon,nlat))
   allocate(sss_M(2,nlon,nlat))
   allocate(dsss_tot_40km(2,nlon,nlat))
   allocate(dsss_tot_70km(2,nlon,nlat))
   allocate(dsss_40km(2,nlon,nlat,nuncertainties))
   allocate(dsss_70km(2,nlon,nlat,nuncertainties))
   allocate(tb_sur0_save(4,2,nlon,nlat))
   allocate(dsss_v(2,nlon,nlat))
   allocate(dsss_h(2,nlon,nlat))

   allocate(gice_est(nlon,nlat))
   allocate(icezone(nlon,nlat))
   allocate(iceflag_amsr2(3,nlon,nlat))
   allocate(dtb_sea_ice_corr(2,nlon,nlat))
   allocate(dtb_sea_ice_error_est(2,nlon,nlat))

   return
end subroutine allocate_L2C_arrays
