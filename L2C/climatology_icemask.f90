! TM 03/09/2022

! climatology icemask (ice4) put onto SMAP L2C grid
! saved to iceflag_amsr2(1,:,:)

subroutine climatology_icemask
   use l2c_module_smap
   implicit none


   integer(4)              ::  idir, ilon, ilat, iice
   integer(4)              ::  lyear,idayjl,imon,idaymo !year, julian day, month (1-12), day of month (1=31))
   real(8)                 ::  secyr,secdy              !seconds in year, second of day
   real(4)                 ::  xlat, xlon

! default
   iceflag_amsr2(1,:,:) = 0
   iceflag_amsr2(3,:,:) = 0


   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            call fd_date_2000(time(idir,ilon,ilat),   secyr,lyear,idayjl,imon,idaymo,secdy)

            xlat= cellat(idir,ilon,ilat)
            xlon= cellon(idir,ilon,ilat)

            call FDICE4(IMON,XLAT,XLON, IICE)
            if (iice == 1) iceflag_amsr2(1,ilon,ilat) = 1

         enddo ! idir
      enddo ! ilon
   enddo ! ilat

   return
end subroutine climatology_icemask
