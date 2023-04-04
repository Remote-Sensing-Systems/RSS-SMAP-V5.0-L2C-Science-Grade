!     this routine makes a small adjustment for the reflected galaxy and sun radiation
!     the refelcted galaxy and sun tables were calculated for a nominal ocean reflectivity and atmosphere and no Faraday rotation
!     this routine adjusts for the actual reflectivity
!     adapted for SMAP
!     TM 04/27/2015


subroutine adjust_tagal_ref(refl0,refl,transq0,transq,taearth,tagal_ref, tagal_ref_adj)
   use l2c_module_smap
   implicit none

   real(4),    intent(in)      :: refl0(2)           ! nominal reflectivity
   real(4),    intent(in)      :: refl(2)            ! actual  reflectivity
   real(4),    intent(in)      :: transq0            ! nominal tran**2
   real(4),    intent(in)      :: transq             ! actual  tran**2
   real(4),    intent(in)      :: taearth(4)         ! TA of Earth radiation BEFORE the galaxy adjsutment
                                                     ! this is needed to get an approximate value for Faraday rotation
   real(4),    intent(in)      :: tagal_ref(3)       ! value of TA_gal_ref or TA_sun_ref from the tables

   real(4),    intent(out)     :: tagal_ref_adj(3)   ! actual vlaue of TA_gal_ref or TA_sun_ref



   integer(4)                  :: istokes
   real(4)                     :: tbtoi(4),tb_gal_toa(3),tb_gal_toi(3), faraday_deg
   real(4)                     :: tbv,tbh

!     ===============================================================================================
!     ==================== find estimates of earth tb toi and faraday rotation angle ================
!     ===============================================================================================

   tbtoi = real(matmul(amat_IQ, taearth), 4)
   faraday_deg=real(0.5*atan(tbtoi(3)/tbtoi(2))/rad, 4)

!     =================================================================================================================
!     ========== convert normalized tagal to tbgal at toa (tagal_ref is found assuming no faraday rotation) ===========
!     =================================================================================================================
   do istokes=1,2
      tb_gal_toa(istokes)=real(dot_product(amat_IQ(istokes,1:3),tagal_ref(1:3)), 4)
   enddo
   tb_gal_toa(3)=0.0


!     =================================================================================================================
!     ==================== convert normalized tbgal to tbgal corresponding to actual env conditions ===================
!     =================================================================================================================

   tbv=0.5*(tb_gal_toa(1)+tb_gal_toa(2))
   tbh=0.5*(tb_gal_toa(1)-tb_gal_toa(2))

   tbv=tbv*refl(1)/refl0(1)
   tbh=tbh*refl(2)/refl0(2)

   tb_gal_toa(1)=tbv+tbh
   tb_gal_toa(2)=tbv-tbh

   tb_gal_toa(1:2)=tb_gal_toa(1:2)*transq/transq0

!     =================================================================================================================
!     =============================== apply faraday rotation to tbgal =================================================
!     =================================================================================================================

   tb_gal_toi(1)=tb_gal_toa(1)
   tb_gal_toi(2)=real(tb_gal_toa(2)*cos((2*faraday_deg)*rad), 4)
   tb_gal_toi(3)=real(tb_gal_toa(2)*sin((2*faraday_deg)*rad), 4)

!     =================================================================================================================
!     =============================== convert tb toi galatic radiation to antenna temperature =========================
!     =================================================================================================================

   do istokes=1,3
      tagal_ref_adj(istokes)=real(dot_product(amat_inv_IQ(istokes,1:3),tb_gal_toi(1:3)), 4)
   enddo

   return
end subroutine adjust_tagal_ref
