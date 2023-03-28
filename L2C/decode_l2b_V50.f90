!     Version 5
!     added AMSR2 sea-ice
!     took out SWH

!     Version 3
!     added IMERG rain rate as last record

!     same as Frank's routine
!     last edited on 09/24/2015: added sun angles


!     I am reducing the recorded vlaues for the sea-ice correction and gice_est by 28% in zones 1,2,3
!     When running it with the computed values, I found a slight over-correction

!
subroutine decode_l2b(iorbit,filename_l2b, start_time, ierr)
   use l2c_module_smap
   implicit none

   integer(4), intent(in)      ::   iorbit
   character(120),intent(in)   ::   filename_l2b
   real(8) , intent(out)       ::   start_time
   integer(4), intent(out)     ::   ierr

   integer(4)                  ::   iorbitx
   real(4)                     ::   buf(84) ! changed from 80 in V5
   integer(4)                  ::   jbuf(4) ! added in V5

   integer(4)                  ::   jdir,idir,icel,nout
   integer(4)                  ::   ilat,ilon
   integer(4)                  ::   itrg_tot,itrg,kscan

   ierr=0

   open(unit=2,file=filename_l2b,status='old',action='read',form='binary')

   read(2) vernum,iorbitx,start_time,xlon_node1,xlon_node2,klon_end
   if(iorbitx.ne.iorbit) then
      write(*,*) iorbit,iorbitx
      write(*,*)' iorbit out of sync in smap_l2b_processor_rt, pgm error'
      ierr=1
      return
   endif

   ifill=0

   do jdir=1,2    !1 is fore look, 2 is aft look

      read(2) idir,nout
      if( idir.ne. jdir) then
         write(*,*) idir, jdir
         write(*,*) ' idir out of sync, pgm error'
         ierr=1
         return
      endif

      do icel=1,nout
         read(2) buf
         read(2) jbuf
         ilat=nint(buf(1))
         ilon=nint(buf(2))
         if(ilat.lt.1 .or. ilat.gt.nlat) then
            write(*,*) ilat
            write(*,*) ' ilat oob in decode_l2b, pgm error'
            ierr=1
            return
         endif
         if(ilon.lt.1 .or. ilon.gt.nlon) then
            write(*,*) ilon
            write(*,*) ' ilon oob in decode_l2b, pgm error'
            ierr=1
            return
         endif

         ifill( idir,ilon,ilat)=1

         itrg_tot=nint(buf(3))
         itrg=1 + int((itrg_tot-1)/7)
         kscan=itrg_tot - 7*(itrg-1) - 4
         alpha( idir,ilon,ilat)=0.5*(itrg-1)

         iscan( idir,ilon,ilat)              =buf( 4)
         dist(  idir,ilon,ilat)              =buf( 5)
         time(  idir,ilon,ilat)              =buf( 6) + start_time
         zang(  idir,ilon,ilat)              =buf( 7)
         cellat(idir,ilon,ilat)              =buf( 8)
         cellon(idir,ilon,ilat)              =buf( 9)
         eia(   idir,ilon,ilat)              =buf(10)
         eaa(   idir,ilon,ilat)              =buf(11)
         pra(   idir,ilon,ilat)              =buf(12)
         sunglt(idir,ilon,ilat)              =buf(13)
         monglt(idir,ilon,ilat)              =buf(14)
         gallat(idir,ilon,ilat)              =buf(15)
         gallon(idir,ilon,ilat)              =buf(16)
         frdrot(idir,ilon,ilat)              =buf(17)  !faraday rotation for tec=1
         teclat(idir,ilon,ilat)              =buf(18)
         teclon(idir,ilon,ilat)              =buf(19)

         loss_ant(   :,idir,ilon,ilat)       =buf(20:21)
         temp_ant(   :,idir,ilon,ilat)       =buf(22:23)

         ! added sun angles
         sun_alpha(idir,ilon,ilat)           =buf(24)
         sun_beta(idir,ilon,ilat)            =buf(25)

         ta_ant(    :,idir,ilon,ilat)            =buf(26:29) ! this was called ta_ert in Frank's routine.
         ta_ant_filtered(    :,idir,ilon,ilat)   =buf(30:33) ! added ta filtered.

         ta_sun_dir(:,idir,ilon,ilat)        =buf(34:36)
         ta_sun_ref(:,idir,ilon,ilat)        =buf(37:39)

         ta_gal_dir(:,idir,ilon,ilat)        =buf(40:42)
         ta_gal_ref_tab(:,idir,ilon,ilat,1)  =buf(43:45) !W=0
         ta_gal_ref_tab(:,idir,ilon,ilat,2)  =buf(46:48) !W=5
         ta_gal_ref_tab(:,idir,ilon,ilat,3)  =buf(49:51) !W=10
         ta_gal_ref_tab(:,idir,ilon,ilat,4)  =buf(52:54) !W=15
         ta_gal_ref_tab(:,idir,ilon,ilat,5)  =buf(55:57) !W=20

         ta_lnd(    :,idir,ilon,ilat)        =buf(58:59)
         fland(       idir,ilon,ilat)        =buf(60)
         gland(       idir,ilon,ilat)        =buf(61)
         wt_sum(      idir,ilon,ilat)        =buf(62)

         tec(         idir,ilon,ilat)        =buf(63)  !is different for fore vs aft

         !  remaining parameters are the same for fore and aft
         sss_ref(    ilon,ilat)               =buf(64)

         surtep(     ilon,ilat)               =buf(65)

         winspd_anc(ilon,ilat)                =buf(66)
         windir(     ilon,ilat)               =buf(67)

         absp_oxy(   ilon,ilat)               =buf(68)  !oxy absorption
         absp_vap(   ilon,ilat)               =buf(69)  !vap absorption
         absp_liq(   ilon,ilat)               =buf(70)  !liq absorption
         tran(       ilon,ilat)               =buf(71)
         tbup(       ilon,ilat)               =buf(72)
         tbdw(       ilon,ilat)               =buf(73)

         lst(        ilon,ilat)               =buf(74)
         sm(         ilon,ilat)               =buf(75)

         solar_flux( ilon,ilat)               =buf(76)

         ! added in V3
         rain(       ilon,ilat)               =buf(77)

         ! added in V5
         sst_err(    ilon,ilat)               =buf(78)
         xobs_sat_ccmp(ilon,ilat)             =buf(79)

         ! V5 AMSR2 sea-ice
         gice_est(   ilon,ilat)               =buf(80)

         dtb_sea_ice_corr(     1:2,ilon,ilat)        =buf(81:82)
         dtb_sea_ice_error_est(1:2,ilon,ilat)        =buf(83:84)

         icezone(           ilon,ilat)        =jbuf(1)
         iceflag_amsr2(1:3, ilon,ilat)        =jbuf(2:4)


      enddo  !icel
   enddo  !jdir
   close(2)

! V5 invalid time
   where(time < start_time)
      time = 0.d0
      ifill = 0
      cellat = missing_val4
      cellon = missing_val4
   endwhere

! adjust sea-ice correction in zones 1,2,3
! lower sea-ice error estimate in zones 3 and 4
   do ilat=1,nlat
      do ilon=1,nlon

         if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
         if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

         if (icezone(ilon,ilat)>=1 .and. icezone(ilon,ilat)<4) then
            gice_est(ilon,ilat)             = gice_est(ilon,ilat) * sea_ice_adj_1
            dtb_sea_ice_corr(1:2,ilon,ilat) = dtb_sea_ice_corr(1:2,ilon,ilat) * sea_ice_adj_1
         endif

         if (icezone(ilon,ilat)==3) dtb_sea_ice_error_est(1:2,ilon,ilat) = dtb_sea_ice_error_est(1:2,ilon,ilat)*sea_ice_adj_2
         if (icezone(ilon,ilat)==4) dtb_sea_ice_error_est(1:2,ilon,ilat) = dtb_sea_ice_error_est(1:2,ilon,ilat)*sea_ice_adj_3



      enddo
   enddo

   return
end subroutine decode_l2b
