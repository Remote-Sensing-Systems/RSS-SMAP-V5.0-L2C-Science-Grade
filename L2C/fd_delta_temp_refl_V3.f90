! TM 03/22/2018
! adjustment to reflector temperature for SMAP V3/V4/V5

subroutine fd_delta_temp_refl_V3(secyr, xang,  delta_temp_refl)
   use external_files_l2c_module
   use dir_paths_module
   use l2c_module_smap

   implicit none
! interpolate physical reflector temperature from tabulated values
   real(8), intent(in)         ::  secyr
   real(4), intent(in)         ::  xang
   real(8), dimension(2), intent(out)  :: delta_temp_refl

   integer(4), save            :: istart=1


! tabulated reflector temperature
! time dimensions changed from V2.
! in V3, a time wrap is performed at both ends
   real(8), dimension(2,0:366,-1:360), save  :: temp_refl_tab

   real(4)                     :: yang
   integer(4)                  :: k0, k1, i0, i1
   real(8)                     :: brief_k, brief_i
   real(8), dimension(2)       :: t0, t1, t

   real(8), parameter          :: secyr_tot = 86400.d0*365.d0
   real(8)                     :: ysecyr, xday


   if (istart==1) then
      istart=0
      ! read in tabulated reflector temperature
      call get_dir_paths(table_dir, data_dir)
      open(unit=3,file=trim(table_dir) // '/' // temp_tab_file,form='unformatted',access='stream',action='read',status='old')
      read(3) temp_refl_tab(:,:,0:360)
      temp_refl_tab(:,:,-1) = temp_refl_tab(:,:,359) ! full wrap
      close(3)
   endif


   yang=xang
   if (yang < 0.001)   yang=0.001
   if (yang > 359.999) yang=359.999
   i0 = floor(yang-0.5) ! centered on half angle
   i1 = i0 + 1
   brief_i = yang-i0-0.5


   ysecyr=secyr
   if (ysecyr>secyr_tot-0.1) ysecyr=secyr_tot-0.1 ! leap year: set to last day, last second

   xday = ysecyr/86400.d0 + 1
   k0 = floor(xday-0.5) ! centered on half day
   if (k0==366) k0=365  ! leap year
   if (k0<0 .or. k0>365) then
      write(*,*) secyr,xday,k0
      stop 'index k0 invalid in subroutine fd_delta_temp_refl'
   endif
   k1 = k0 + 1
   brief_k = xday-k0-0.5


   t0(:) = temp_refl_tab(:,k0,i0)*(1.0-brief_k) + temp_refl_tab(:,k1,i0)*brief_k
   t1(:) = temp_refl_tab(:,k0,i1)*(1.0-brief_k) + temp_refl_tab(:,k1,i1)*brief_k

   t(:) = t0(:)*(1.0-brief_i) + t1(:)*brief_i

   delta_temp_refl(:) = t(:)

   return
end subroutine fd_delta_temp_refl_V3
