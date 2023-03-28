! TM 04/25/2015
! TA antenna -> TA EARTH
! removes galaxy, sun, moon, CS

!
! iopt parameters
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref

subroutine fd_ta_earth
   use l2c_module_smap
   use SMAP_ROUGHNESS_GMF_V3B_module
   implicit none


   integer(4)              ::  idir, ilon, ilat
   integer(4)              ::  iwin1, iwin2
   real(4), dimension(4)   ::  tamea, taearth
   real(4)                 ::  wind, brief, c1, c2

   real(4), parameter      ::  wsss0=35.0, wsst0=20.0, transq0=1.0  !nominal values use when generating reflected galaxy tables
   real(4), dimension(2)   ::  xbuf, refl0, refl
   real(4)                 ::  tht, xlat, xlon, phir, transq
   real(4), dimension(3)   ::  ta_ref_adj

   integer(4)              ::  lyear,idayjl,imon,idaymo !year, julian day, month (1-12), day of month (1-31)
   real(8)                 ::  secyr,secdy              !seconds in year, second of day

   real(4)                 ::  xta_gal_refl_tab(3,5), xgallat, xgallon, xta_gal_refl(2)

   real(4)                 ::  zice


   ta_earth = missing_val4  ! initialize as missing


   do ilat=1,nlat
      do ilon=1,nlon
         do idir=1,2

            if(ifill(idir,ilon,ilat) == 0) cycle ! missing observation
            if (abs(wt_sum(idir,ilon,ilat)-1.0)>0.01) cycle

            call fd_date_2000(time(idir,ilon,ilat),   secyr,lyear,idayjl,imon,idaymo,secdy)


            ! reflected galaxy
            ! need to interpolate to actual wind speed
            wind = winspd(ilon,ilat)
            if (wind < 0.0) wind=0.0

            if (gland(idir,ilon,ilat)>0.1) wind = 0.0  ! over land use specular value
            if (fland(idir,ilon,ilat)>0.1) wind = 0.0  ! over land use specular value

            ! changed in V5
            ! over sea ice use specular value
            if (icezone(ilon,ilat)==5) wind=0.0
            if (icezone(ilon,ilat)==6 .and. iceflag_amsr2(2,ilon,ilat)/=0) wind=0.0


            if (igal_wspd ==1) then ! use GO galaxy and NCEP wind speed
               brief=wind/5.
               if(brief.gt.3.999) brief=3.999
               iwin1=1+floor(brief)
               if (iwin1<1) iwin1=1
               iwin2=iwin1+1
               c1=iwin1-brief
               c2=1.0-c1
               ta_gal_ref(:,idir,ilon,ilat) = c1*ta_gal_ref_tab(:,idir,ilon,ilat,iwin1) + c2*ta_gal_ref_tab(:,idir,ilon,ilat,iwin2)
               ! interpolate
            else if (igal_wspd ==2) then ! use Frank's adjustement that he derived from SMAP
               xta_gal_refl_tab(1:3,1:5) = ta_gal_ref_tab(1:3,idir,ilon,ilat,1:5) !IQ Basis
               xgallat = gallat(idir,ilon,ilat)
               xgallon = gallon(idir,ilon,ilat)
               if (xgallon > 359.999) xgallon=359.999
               call find_ta_gal_refl(xta_gal_refl_tab,xgallat,xgallon,wind, xta_gal_refl)
               ! convert back from VH to IQ
               ta_gal_ref(1,idir,ilon,ilat) = xta_gal_refl(1)+xta_gal_refl(2)
               ta_gal_ref(2,idir,ilon,ilat) = xta_gal_refl(1)-xta_gal_refl(2)
               ta_gal_ref(3,idir,ilon,ilat) = 0.0
            endif




            ! sun. need to multiply be actual value for solar flux
            ta_sun_dir(:,idir,ilon,ilat) = ta_sun_dir(:,idir,ilon,ilat)*solar_flux(ilon,ilat) !convert to actual flux for current day
            ta_sun_ref(:,idir,ilon,ilat) = ta_sun_ref(:,idir,ilon,ilat)*solar_flux(ilon,ilat) !convert to actual flux for current day

            tamea(1:4) = ta_ant_calibrated(1:4,idir,ilon,ilat) ! chaged from ta_ant_filtered(1:4,idir,ilon,ilat)
            xbuf(1:2)  = tamea(1:2)
            call vh_to_stokes(xbuf)
            tamea(1:2) = xbuf(1:2)

            taearth(1:3) = tamea(1:3) - iopt(2)*ta_gal_dir(1:3,idir,ilon,ilat) - iopt(4)*ta_sun_dir(1:3,idir,ilon,ilat)
            taearth(4)   = tamea(4)
            ! the tagal_dir does contain the isotropic cosmic background cold space
            ! when computing the tagal_ref tables Frank has subtracted out the isotropic cosmic background cold space
            ! assuming an arbitrary value of 3 K.
            ! this is part of the downwelling reflected radiation

            tht = eia(idir,ilon,ilat)
            xlat= cellat(idir,ilon,ilat)
            xlon= cellon(idir,ilon,ilat)
            phir=eaa(idir,ilon,ilat)-windir(ilon,ilat) !relative wind direction angle, 0=upwind
            ! compute actual reflectivity
            zice=0.0
            if (icezone(ilon,ilat)==5) zice=1.0
            if (icezone(ilon,ilat)==6 .and. iceflag_amsr2(2,ilon,ilat)/=0) zice=1.0
            call find_refl_RTM(idayjl=idayjl,xlat=xlat,xlon=xlon,tht=tht,phir=phir,&
               frac_land=fland(idir,ilon,ilat),frac_ice=zice,sss=sss_ref(ilon,ilat),&
               surtep=surtep(ilon,ilat),wind=winspd(ilon,ilat),lst=lst(ilon,ilat),sm=sm(ilon,ilat), &
               refl_tot=refl)

            ! nominal reflectivity and tran**2 used when computing galaxy refl. table
            call fd_water_refl(tht=tht,sss=wsss0,sst=wsst0,wind=wind, refl=refl0)  !refl. coef. assume for galaxy table
            transq = tran(ilon,ilat)**2 ! actual tran**2
            call adjust_tagal_ref(refl0,refl,transq0,transq,taearth,ta_gal_ref(1:3,idir,ilon,ilat), ta_ref_adj)
            ta_gal_ref(1:3,idir,ilon,ilat)=ta_ref_adj(1:3)

            ! nominal reflectivity and tran**2 used when computing sun refl. table
            call fd_water_refl(tht=tht,sss=wsss0,sst=wsst0,wind=0.0, refl=refl0)  !refl. coef. assume for galaxy table
            transq = tran(ilon,ilat)**2 ! actual tran**2
            call adjust_tagal_ref(refl0,refl,transq0,transq,taearth,ta_sun_ref(1:3,idir,ilon,ilat), ta_ref_adj)
            ta_sun_ref(1:3,idir,ilon,ilat)=ta_ref_adj(1:3)

            taearth(1:3) = taearth(1:3) - iopt(3)*ta_gal_ref(1:3,idir,ilon,ilat) - iopt(5)*ta_sun_ref(1:3,idir,ilon,ilat)
            taearth(4)   = tamea(4)

            ! store as V/H
            xbuf(1:2)  = taearth(1:2)
            call stokes_to_vh(xbuf)
            taearth(1:2) = xbuf
            ta_earth(1:4,idir,ilon,ilat)=taearth(1:4)

         enddo ! idir
      enddo ! ilon
   enddo ! ilat



   return
end subroutine fd_ta_earth
