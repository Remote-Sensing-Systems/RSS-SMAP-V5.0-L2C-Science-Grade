module  external_files_l2c_module
   implicit none

! APC from Earth via integration over AP
! TM March 2022
   character(len=*), parameter   ::   apc_file = 'SMAP_AMATRIX_V50_TM_MARCH_2022_EARTH_FOV.txt'

! correction for reflector temperature
   character(len=*), parameter   ::   temp_tab_file = 'delta_refl_temp_V3A.dat'

! adjsutment to galactic map based on Frank's fore-aft analysis
   character(len=*), parameter   ::   gal_adj_map_file = 'smooth_gal_adjustment_maps.dat'

! vegetation tau square table
! Aquarius ATBD
   character(len=*), parameter   ::   tausq_file='mk_vege_tausq_tab.dat'

! surface roughness emissivity table
   character(len=*), parameter   ::   emiss_coeff_harm_file = 'dew_phi_VH34_harmonic_tab.dat'
   character(len=*), parameter   ::   emiss_coeff_harm_file_AQ = 'deW_harm_coeffs_V9A_MI.dat'
   character(len=*), parameter   ::   delta_file = 'delta_EW_V5_B.dat'
   character(len=*), parameter   ::   demiss_res_file = 'dw_res_tab_spline.dat'

! adjusted land correction
   character(len=*), parameter   ::   land_delta_tb_file = 'mk_smap_minus_model_tbmap.dat'

! fland and gland gradient map
   character(len=*), parameter   ::   land_gradient_file = 'mk_land_gradient_map.dat'

! tb meas - exp table direct access file
   character(len=*), parameter   ::   dtb_statfile = 'smap_dtb_stats_ocean_V50.dat'

! wind uncertainty estimates
! Carl's file with CCMP - buoy statistics for CCMP RT
   character(len=*), parameter   ::   wind_unc_file = 'CCMP_RT_minus_Buoy_W_sat_no_sat.txt'

! climatological ice mask
   character(len=*), parameter   ::  ICEFILE = 'ice4.dat'

end module external_files_l2c_module


