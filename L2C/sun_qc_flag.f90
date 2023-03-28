 ! Richard Lindsley
 ! sun glint flag for SMAP V4.0 and V5.0

 ! Sunglint QC flag
 !
 ! sunglt: the sun glint angle in degrees
 ! winspd: the wind speed in m/s
 !
 ! output = 1 if the QC flag should be set   = 0 otherwise


subroutine compute_sun_qc(sunglt, winspd,   iflag_sun)
   implicit none

   real(4), intent(in)         :: sunglt, winspd
   integer(4), intent(out)     :: iflag_sun

   real(4)                     :: thresh_wind

   ! Create a customized sunglint QC flag: it's a polynomial
   ! function of sun glint angle between 30 and 50 degrees, where
   ! if the wind speed threshold is larger than the function, the
   ! mask applies. Greater than 50 degrees, it's never masked, and
   ! between 0 and 30 degrees it's always masked.

   if (sunglt < 30.0) then
      iflag_sun = 1
   else if (sunglt >= 50.0) then
      iflag_sun = 0
   else
      ! When the sun glint is between 30 and 50 degrees, compare the
      ! wind speed against the threshold level
      thresh_wind = 1.0 / 8000.0 * (sunglt - 30.0)**4

      if (winspd > thresh_wind) then
         iflag_sun = 1
      else
         iflag_sun = 0
      endif

   endif

   return
end subroutine compute_sun_qc

