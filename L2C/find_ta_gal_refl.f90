subroutine find_ta_gal_refl(ta_gal_refl_tab,gallat,gallon,wind, ta_gal_refl)
   use external_files_l2c_module
   use dir_paths_module

   implicit none
   ! this is Frank's routine that outputs the adjusted glalxy based on SMAP for-aft
   ! the outpout of Frank's routine is in V/H


   real(4)         ::  ta_gal_refl_tab(3,5),gallat,gallon,wind,  ta_gal_refl(2)
   real(4), save   ::  dgalta_smooth(1440,720,2,4)

   integer(4),save ::  istart=1
   integer(4)      ::  iwin1,iwin2,jwin1,jwin2,ilatg,ilong,ipol
   real(4)         ::  brief,a1,a2,b1,b2,ta1,ta2,dta1,dta2


   if(istart.eq.1) then
      istart=0
      call get_dir_paths(table_dir, data_dir)
      call openbig(3,trim(table_dir) // '/' // gal_adj_map_file,'old')
      read(3) dgalta_smooth
      close(3)
   endif

   brief=(wind+2)/5.  !add 2ms to wind
   if(brief.gt.3.999) brief=3.999
   iwin1=int(1+brief, 4)
   iwin2=iwin1+1
   a1=iwin1-brief
   a2=1-a1


   if(wind.lt.5) then
      brief=(wind-1.5)/3.5
      if(brief.lt.0) brief=0
      jwin1=int(1+brief, 4)
      jwin2=jwin1+1
      b1=jwin1-brief
      b2=1-b1
   else
      brief=wind/5.
      if(brief.gt.2.999) brief=2.999
      jwin1=int(1+brief, 4)
      jwin2=jwin1+1
      b1=jwin1-brief
      b2=1-b1
   endif

   ilatg=1+int(4*(gallat+90.))
   ilong=1+int(4* gallon)
   if(ilatg.lt.1 .or. ilatg.gt. 720) then
      write(*,*) gallat,ilatg
      stop 'error in ilatg in '
   endif

   if(ilong.lt.1 .or. ilong.gt.1440) then
      write(*,*) gallon,ilong
      stop 'error in ilong'
   endif

   do ipol=1,2

      if(ipol.eq.1) then
         ta1= 0.5*(ta_gal_refl_tab(1,iwin1) + ta_gal_refl_tab(2,iwin1))
         ta2= 0.5*(ta_gal_refl_tab(1,iwin2) + ta_gal_refl_tab(2,iwin2))
      else
         ta1= 0.5*(ta_gal_refl_tab(1,iwin1) - ta_gal_refl_tab(2,iwin1))
         ta2= 0.5*(ta_gal_refl_tab(1,iwin2) - ta_gal_refl_tab(2,iwin2))
      endif

      dta1=dgalta_smooth(ilong,ilatg,ipol,jwin1)
      dta2=dgalta_smooth(ilong,ilatg,ipol,jwin2)

      ta_gal_refl(ipol)=a1*ta1 + a2*ta2 + b1*dta1 + b2*dta2
   enddo  !ipol

   return
end subroutine find_ta_gal_refl

