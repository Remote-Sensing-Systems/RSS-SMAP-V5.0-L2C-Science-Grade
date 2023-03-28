
!     surface temperature (K)
!     dtb to be added to tb land correction, which is then subtracted for tb
!     dsss to be added to salinity retrieval

      subroutine land_corr_step2(cellat,cellon,zang,gland,imonth,  dtb)
         use external_files_L2C_module
         use dir_paths_module
         implicit none

         integer(4), parameter   ::  nres=4
         integer(4), parameter   ::  mlon=360*nres, mlat=180*nres

         integer(4) imonth
         real(4) cellat,cellon,zang,gland, dtb(2)

         integer(4) istart
         integer(4) ilat,ilon,iasc,jlat,jlon

!     tables that are read in
         real(4) smap_minus_mod_tb(2,mlon,mlat,2,12)
         real(4) grad(2,2880,1441)

         data istart/1/

         if(istart.eq.1) then
            istart=0
            call get_dir_paths(table_dir, data_dir)

!     =================================================================================================
!     ============================ read in smap minus model tb diffs  =================================
!     =================================================================================================

            call openbig(3,trim(table_dir) // '/' // land_delta_tb_file,'old')
            read(3) smap_minus_mod_tb
            close(3)

!     =================================================================================================
!     ============================ read in fland,gland gradient map   =================================
!     =================================================================================================
            call openbig(3,trim(table_dir) // '/' // land_gradient_file,'old')
            read(3) grad
            close(3)

         endif

         ilat=1+int(4*(cellat+90.))
         ilon=1+int(4*cellon)
         if(ilat.lt.1 .or. ilat.gt.mlat) then
            write(*,*) '1st ', ilat, cellat
            stop 'ilat oob, pgm stopped'
         endif
         if(ilon.lt.1 .or. ilon.gt.mlon) then
            write(*,*) '1st ', ilon, cellon
            stop 'ilon oob, pgm stopped'
         endif

         jlat=nint(8.*cellat)+721
         jlon=nint(8.*cellon)+1
         if (jlon .eq. 2881) jlon=1
         if(jlat.lt.1 .or. jlat.gt.1441) then
            write(*,*) '2nd', ilat, cellat
            stop 'ilat oob, pgm stopped'
         endif
         if(jlon.lt.1 .or. jlon.gt.2880) then
            write(*,*) '2nd ', ilon, cellon
            stop 'ilon oob, pgm stopped'
         endif

         iasc=1
         if(zang.gt.180) iasc=2

         ! the second term accounts for the difference between land model TB and land SMAP TB
         ! the first term is an adjsutment that accounts for the residual gradient in the 1/8 deg land field
         ! see 'O:/smap/land_correction/mk_land_gradient_map.f'
         ! 1st index in grad = gradient of fland   2nd index in grad = gradient of gland
         ! typical TB land - TB ocean: V-pol 150K   h-pol 160K

         dtb(1)= 150.*grad(2,jlon,jlat) + gland*smap_minus_mod_tb(1,ilon,ilat,iasc,imonth)
         dtb(2)= 160.*grad(2,jlon,jlat) + gland*smap_minus_mod_tb(2,ilon,ilat,iasc,imonth)


         ! adjustment at high gland
         if(gland.gt.0.02) then
            dtb(1) = dtb(1) + 78*(gland-0.02)
            dtb(2) = dtb(2) + 83*(gland-0.02)
         endif

         return
      end subroutine land_corr_step2

