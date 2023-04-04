program MAKE_SMAP_L2C_V50_ascii
   use l2c_module_smap
   use dir_paths_module

   implicit none

   character(len=250)      ::  filename_l2b

   integer(4)              ::  iorbit1,iorbit2, ibad
   character(len=5)        ::  corbit1, corbit2
   character(len=2)        ::  cpixel

   integer(4)              ::  iorbit,ires_opt,ires

   logical(4)              ::  lexist

   real(8)                 ::  start_time

   real(8)                 ::  secyr,secdy
   integer(4)              ::  lyear,idayjl,imon,idaymo,isecdy
   real(4)                 ::  xhour

   integer(4)              ::  iflag
   real(8), dimension(4)   ::  dtf_bias
   real(8), dimension(2)   ::  tf_ave


   character(len=10),dimension(3)          ::  sbuf
   integer(4), dimension(8)                ::  date_time

   integer(4)              ::  ipixel
   integer(4)              ::  cnt
   integer(4)              ::  ilon, ilat, idir


   write(*,*) ' V5.0 SSS L2C processing'

! begin and end orbit are passed through command line
   call get_command_argument(1, corbit1)
   call get_command_argument(2, corbit2)
   call get_command_argument(3, cpixel)

   cnt = command_argument_count()

   if (cnt .eq. 3) then
      read(corbit1,*) iorbit1
      read(corbit2,*) iorbit2
      read(cpixel, *) ipixel
   else
      ! Incorrect number of command line arguments
      error stop "incorrect number of command line arguments. Must be: start orbit, end orbit, pixel number"
   end if


   write(*,'(a18,1x,i5.5)') ' begin orbit ',iorbit1
   write(*,'(a18,1x,i5.5)') ' end orbit '  ,iorbit2

   write(*,*)


! start processing
   iopt=(/1,1,1,1,1,  0,0,0,0,0/)
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref

   igal_wspd = 2
! =1 use provided wind speed for galaxy roughness
! =2 Frank's new galaxy model

   ires_opt=1
   ires=40
   write(*,*) ' start SMAP FINAL SSS L2C PROCESSING V5.0'
   write(*,*)

   call allocate_L2C_arrays

   call initialize_APC_V5


   do iorbit=iorbit1,iorbit2

      call check_orbit(iorbit, ibad)
      if(ibad.eq.1)     then
         write(*,8866) iorbit,' in bad orbit list. skip orbit'
8866     format(1x,i5.5,a75)
         cycle
      endif
      if(ibad.eq.2)     then
         write(*,8867) iorbit,' in bad ancillary orbit list. skip orbit'
8867     format(1x,i5.5,a75)
         cycle
      endif

      call get_filename_l2B_ascii(iorbit, ipixel, filename_l2b)
      inquire(file=filename_l2B,exist=lexist)
      if(.not.(lexist))  then
         write(*,*) filename_l2b
         write(*,*) ' L2B file does not exist. no L2C processing.'
         cycle
      endif

      write(*,*)
      write(*,*) ' processing orbit ',iorbit

      CALL DATE_AND_TIME (sbuf(1),sbuf(2),sbuf(3),date_time)
      write(*,8857) iorbit,' SSS V5.0 L2C processing start     ',&
         date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
8857  format(1x,i5.5,a35,i4.4,'-',i2.2,'-',i2.2,' T',i2.2,':',i2.2,':',i2.2)

      call read_l2b_ascii(filename_l2b, ilon, ilat, idir)

      write(*,*) '     calibration and reflector emissivity'
      call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
      isecdy=nint(secdy)
      xhour = real(secdy/3600., 4)

      ! correct reflector temperature
      call correct_reflector_emissivity_V30

      ! ocean target calibration and caibration drift correction
      call find_dtb_bias(iorbit,    dtf_bias, tf_ave, iflag)
      if (iflag /=0 ) then
         write(*,*) iorbit,' has no dtb_bias records. skip'
         cycle
      endif

      dtb_bias_orbit(1:4) = real(dtf_bias(1:4), 4)  !save
      tf_ave_orbit(1:2)   = real(tf_ave(1:2), 4)

      call correct_cal_drift_V30

      call climatology_icemask ! remap ice4 climatology on SMAP swath

      write(*,*) '     SSS retrieval'

      winspd=winspd_anc ! this is the CCMP wind speed (for SMAP SSS V5 files)

      call fd_ta_earth
      call fd_tb_toa_lc
      call fd_tb_sur_sic

      call fd_ta_expected

      call fd_sss

      call create_l2_qcflag

      call write_l2c_ascii(iorbit, ipixel, ilon, ilat, idir)

   enddo ! iorbit
   write(*,*)
   stop ' normal end MAKE_SMAP_L2C_V50_ascii '
end program MAKE_SMAP_L2C_V50_ascii