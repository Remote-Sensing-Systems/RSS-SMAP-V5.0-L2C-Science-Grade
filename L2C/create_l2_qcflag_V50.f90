! adapted to V5.0

subroutine create_l2_qcflag
   use l2c_module_smap
   implicit none

   integer(4)              ::  idir, ilon, ilat, iflag_sun
   real(4)                 ::  xsst

! initialize
   iqc_flag = 0

   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            ! bit 0: no or invalid radiometer observation in cell
            if(ifill(idir,ilon,ilat) == 0) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),0)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif
            if (abs(ta_ant_filtered(1,idir,ilon,ilat)-missing_val4)<0.1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),0)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif
            if (abs(ta_ant_filtered(2,idir,ilon,ilat)-missing_val4)<0.1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),0)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif
            if (abs(ta_ant_filtered(3,idir,ilon,ilat)-missing_val4)<0.1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),0)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif

            ! bit1: problem with OI
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),1)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif

            ! bit2: strong land contamination
            if (gland(idir,ilon,ilat) > 0.1  .or. fland(idir,ilon,ilat) > 0.1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),2)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif

            ! bit3: strong sea ice contamination
            ! use new sea-ice flag. Changed in V5
            ! ice zone 5
            ! ice zone 6 (AMSR2 VL resolution close to land) and standard AMSR2 V8.2 flag set
            if (icezone(ilon,ilat)==5) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),3)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif
            if (icezone(ilon,ilat)==6 .and. iceflag_amsr2(2,ilon,ilat)==1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),3)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif

            ! bit4: MLE in SSS retrieval algo has not converged
            ! changed in V5. Also flagged if sss is < 1.5 or sss > 43.5
            if (iflag_sss_conv(idir,ilon,ilat) /=0) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),4)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
               cycle
            endif
            if (sss_smap_40km(idir,ilon,ilat) <= sss_val_min  .or. sss_smap_40km(idir,ilon,ilat) >= sss_val_max) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),4)
               sss_smap_40km(idir,ilon,ilat) = missing_val4
               cycle
            endif

            ! bit5: sunglint
            ! changed in V4.0
            call compute_sun_qc(sunglt(idir,ilon,ilat), winspd(ilon,ilat),   iflag_sun)
            if (iflag_sun==1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),5)
            endif

            ! bit6: moonglint
            if (abs(monglt(idir,ilon,ilat))<15.) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),6)
            endif

            ! bit7: high reflected galaxy
            if (ta_gal_ref(1,idir,ilon,ilat)/2.0 > 2.0) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),7)
            endif

            ! bit8: moderate land contamination
            if (gland(idir,ilon,ilat) > 0.04  .or. fland(idir,ilon,ilat) > 0.005) then
               ! threshold changed in V4.0
               ! use 0.005 for fland consistently
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),8)
            endif

            ! bit9: moderate sea ice contamination
            ! V5: sea-ice zones  3 and 4
            if (icezone(ilon,ilat)==3 .or. icezone(ilon,ilat)==4) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),9)
            endif

            ! bit10: MLE in SSS retrieval algo. has poorly converged (tb_consistency check)
            if (tb_consistency(idir,ilon,ilat)> 1.0) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),10)
            endif

            ! bit11: low SST
            xsst = surtep(ilon,ilat)-273.15
            if (xsst < 5.0) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),11)
            endif

            ! bit12: high wind
            if (winspd(ilon,ilat) > 15.0) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),12)
            endif

            ! bit13: light land contamination
            if (gland(idir,ilon,ilat)> 0.001) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),13)
            endif

            ! bit14: light sea ice contamination
            ! V5: sea-ice zones 1 and 2
            if (icezone(ilon,ilat)==1 .or. icezone(ilon,ilat)==2) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),14)
            endif

            ! bit15: rain flag (IMERG)
            if (rain(ilon,ilat)>0.1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),15)
            endif

            ! bit16: no sea-ice check or correction possible because AMSR2 TB missing
            if (icezone(ilon,ilat)==7 .and. iceflag_amsr2(1,ilon,ilat)==1) then
               iqc_flag(idir,ilon,ilat) = ibset(iqc_flag(idir,ilon,ilat),16)
               sss_smap_40km(idir,ilon,ilat)=missing_val4
            endif


         enddo ! idir
      enddo ! ilon
   enddo ! ilat


   return
end subroutine create_l2_qcflag

