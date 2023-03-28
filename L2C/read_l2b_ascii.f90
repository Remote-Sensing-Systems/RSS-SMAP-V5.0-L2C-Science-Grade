
subroutine read_l2b_ascii(filename_l2b, ilon, ilat, idir)
   use l2c_module_smap

   implicit none

   character(len=250), intent(in)      ::   filename_l2b
   integer(4), intent(out)             ::   ilon,ilat,idir



   OPEN(10,FILE=filename_l2b,STATUS='old')
   WRITE(*,*) filename_l2b
   read(10,'(2I5,I2,1x,I1,1x,F7.3,I5.3,1x,F7.3,1x,F16.6,1x,78F10.5,I4,3I4)') ilon, ilat, idir, ifill(idir,ilon,ilat),&
   &alpha(idir,ilon,ilat),iscan(idir,ilon,ilat),dist(idir,ilon,ilat), time(idir,ilon,ilat), zang(idir,ilon,ilat),&
   &cellat(idir,ilon,ilat),cellon(idir,ilon,ilat),eia(idir,ilon,ilat),&
   &eaa(idir,ilon,ilat), pra(idir,ilon,ilat), sunglt(idir,ilon,ilat),&
   &monglt(idir,ilon,ilat), gallat(idir,ilon,ilat),gallon(idir,ilon,ilat),frdrot(idir,ilon,ilat), teclat(idir,ilon,ilat),&
   &teclon(idir,ilon,ilat), loss_ant(:,idir,ilon,ilat),temp_ant(:,idir,ilon,ilat),sun_alpha(idir,ilon,ilat),&
   &sun_beta(idir,ilon,ilat),ta_ant(:,idir,ilon,ilat),ta_ant_filtered(:,idir,ilon,ilat),ta_sun_dir(:,idir,ilon,ilat),&
   &ta_sun_ref(:,idir,ilon,ilat),ta_gal_dir(:,idir,ilon,ilat),ta_gal_ref_tab(:,idir,ilon,ilat,1),&
   &ta_gal_ref_tab(:,idir,ilon,ilat,2),ta_gal_ref_tab(:,idir,ilon,ilat,3),ta_gal_ref_tab(:,idir,ilon,ilat,4),&
   &ta_gal_ref_tab(:,idir,ilon,ilat,5),ta_lnd(:,idir,ilon,ilat),fland(idir,ilon,ilat),&
   &gland(idir,ilon,ilat), wt_sum(idir,ilon,ilat),tec(idir,ilon,ilat),sss_ref(ilon,ilat),surtep(ilon,ilat),&
   &winspd_anc(ilon,ilat), windir(ilon,ilat),absp_oxy(ilon,ilat),absp_vap(ilon,ilat),&
   &absp_liq(ilon,ilat),tran(ilon,ilat),tbup(ilon,ilat),tbdw(ilon,ilat),lst(ilon,ilat),sm(ilon,ilat),&
   &solar_flux(ilon,ilat),rain(ilon,ilat),sst_err(ilon,ilat),xobs_sat_ccmp(ilon,ilat),&
   &gice_est(ilon,ilat),dtb_sea_ice_corr(1:2,ilon,ilat),dtb_sea_ice_error_est(1:2,ilon,ilat),&
   &icezone(ilon,ilat),iceflag_amsr2(1:3, ilon,ilat)

   close(10)

end subroutine read_l2b_ascii
