! V3.0
! routines for correcting emissive reflector



subroutine correct_reflector_emissivity_V30
! done during L2C processing
! using intepolation for physical reflector temperature
! use 1% for reflector emissivity
   use l2c_module_smap
   implicit none

   integer(4)              :: idir, ilon, ilat
   real(4), dimension(4)   :: xtf

   integer(4)              :: lyear,idayjl,imon,idaymo !year, julian day, month (1-12), day of month (1=31))
   real(8)                 :: secyr,secdy              !seconds in year, second of day

   real(4)                 :: xang
   real(8), dimension(2)   :: delta_temp_refl, temp_refl



   ta_ant_calibrated=ta_ant_filtered
   dtemp_ant=0.0

   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            call fd_date_2000(time(idir,ilon,ilat),   secyr,lyear,idayjl,imon,idaymo,secdy)

            xtf(1:4) = ta_ant_filtered(1:4,idir,ilon,ilat)

            if (abs(xtf(1)-missing_val4)<0.1) cycle
            if (abs(xtf(2)-missing_val4)<0.1) cycle
            if (abs(xtf(3)-missing_val4)<0.1) cycle
            if (abs(xtf(4)-missing_val4)<0.1) cycle

            xang = zang(idir,ilon,ilat)
            if (xang<0.001) xang=0.001
            if (xang>359.999) xang=359.999
            call fd_delta_temp_refl_V3(secyr, xang,  delta_temp_refl)

            temp_refl(1:2) = temp_ant(1:2,idir,ilon,ilat)
            if (abs(temp_refl(1)-missing_val4)<0.1) cycle
            if (abs(temp_refl(2)-missing_val4)<0.1) cycle

            temp_refl(1:2) = temp_refl(1:2) +  delta_temp_refl(1:2)
            dtemp_ant(1:2,idir,ilon,ilat) =  delta_temp_refl(1:2)

            xtf(1)   = (xtf(1) - emiss_refl_adj(1)*temp_refl(1)) / (1.0-emiss_refl_adj(1))
            xtf(2)   = (xtf(2) - emiss_refl_adj(2)*temp_refl(2)) / (1.0-emiss_refl_adj(2))
            xtf(3)   =  xtf(3)/sqrt((1.0-emiss_refl_adj(1))*(1.0-emiss_refl_adj(2)))
            xtf(4)   =  xtf(4)/sqrt((1.0-emiss_refl_adj(1))*(1.0-emiss_refl_adj(2)))

            ta_ant_calibrated(1:4,idir,ilon,ilat)  = xtf(1:4)

         enddo ! idir
      enddo ! ilon
   enddo ! ilat

   return
end subroutine correct_reflector_emissivity_V30
