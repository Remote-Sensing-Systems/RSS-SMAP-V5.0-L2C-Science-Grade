! TM 07/27/2015
! TA -> TA cal drift corrected
! V3.0/V4.0/V5.0


subroutine correct_cal_drift_V30
   use l2c_module_smap
   implicit none


   integer(4)              ::  idir, ilon, ilat
   real(4), dimension(4)   ::  xta, xta_corr
   real(4), parameter      ::  T_DICKE=293.0
   real(4), dimension(2)   ::  afac, bfac


   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            xta(1:4) = ta_ant_calibrated(1:4,idir,ilon,ilat)

            if (abs(xta(1)-missing_val4)<0.1) cycle
            if (abs(xta(2)-missing_val4)<0.1) cycle
            if (abs(xta(3)-missing_val4)<0.1) cycle
            if (abs(xta(4)-missing_val4)<0.1) cycle

            xta_corr(3:4) = xta(3:4) - dtb_bias_orbit(3:4)

            afac(1:2) =dtb_bias_orbit(1:2)/(tf_ave_orbit(1:2)-T_DICKE)
            ! this is a placeholder for changing TND
            ! using this correction makes sure that the AMAZON does not change much when doing the drift correction

            bfac(1:2) = afac(1:2)*(xta(1:2)-T_DICKE)

            xta_corr(1:2) = xta(1:2) - bfac(1:2)

            ta_ant_calibrated(1:4,idir,ilon,ilat) = xta_corr(1:4)

         enddo ! idir
      enddo ! ilon
   enddo ! ilat



   return
end subroutine correct_cal_drift_V30
