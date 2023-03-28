! adapted for 5.0
! TM 01/31/2022

! TM 05/03/2015
! TB SUR0 -> SSS
!

subroutine fd_sss
   use l2c_module_smap
   use sss_module
   implicit none


   integer(4)              ::  idir, ilon, ilat, iflag_sss
   real(4), dimension(2)   ::  tbsur0
   real(4)                 ::  thtadj, sst, xsss, xchisq_sss

! initialize as missing
   sss_smap_40km       = missing_val4
   iflag_sss_conv =-1 ! default for no data
   tb_consistency = missing_val4

   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            sst         = surtep(ilon,ilat) - 273.15
            thtadj      = eia(idir,ilon,ilat)
            tbsur0(1:2) = tb_sur0_sic(1:2,idir,ilon,ilat) ! V5: sea-ice corrected tb_sur0


            ! no salinity retrieval if more than 10% fland or gland
            if (gland(idir,ilon,ilat) > 0.1) cycle
            if (fland(idir,ilon,ilat) > 0.1) cycle ! added in V5

            ! no salinity retrieval in sea ice zone 5
            ! no salinity retrieval in sea ice zone 6 (AMSR2 VL resolution close to land) if ASMR V8.1 iceflag is set
            ! no salinity retrieval in sea ice zone 7 (AMSR2  L resolution close to land) if climatological iceflag is set

            if (icezone (ilon,ilat) == 5) cycle
            if (icezone (ilon,ilat) == 6 .and. iceflag_amsr2(2,ilon,ilat)==1) cycle
            if (icezone (ilon,ilat) == 7 .and. iceflag_amsr2(1,ilon,ilat)==1) cycle

            call sss_MLE(thtadj,sst,tbsur0,   xsss, xchisq_sss, iflag_sss)
            iflag_sss_conv(idir,ilon,ilat) = iflag_sss
            if (xchisq_sss<1.0E-8) xchisq_sss=1.0E-8
            tb_consistency(idir,ilon,ilat) = sqrt(xchisq_sss)

            if (iflag_sss /=0) then
               write(*,*) ilat,ilon,idir
               write(*,*) gland(idir,ilon,ilat), icezone(ilon,ilat)
               write(*,*) tbsur0(1:2)
               write(*,*) thtadj,sst
               write(*,*) xchisq_sss, iflag_sss
               write(*,*) ' sss MLE has not converged'
               write(*,*)
               cycle
            endif

            sss_smap_40km(idir,ilon,ilat) = xsss

         enddo ! idir
      enddo ! ilon
   enddo ! ilat

   return
end subroutine fd_sss
