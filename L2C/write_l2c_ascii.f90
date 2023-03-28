subroutine write_l2c_ascii(iorbit, ipixel, ilon, ilat, idir)
   use l2c_module_smap
   use dir_paths_module

   implicit none

   integer(4), intent(in) :: iorbit
   integer(4), intent(in) :: ipixel
   integer(4), intent(in) :: ilon, ilat, idir
   character(len=250)     :: filename_l2c
   integer(4)             :: i

   call get_dir_paths(table_dir, data_dir)
   write(filename_l2c, 2002) trim(data_dir), iorbit, ipixel
2002 format(A,'/sample_l2c_',i5,'_',i2.2,'_output.txt')
   OPEN(10,FILE=filename_l2c,STATUS='REPLACE')
   WRITE(*,*) 'writing', filename_l2c


   write(10,'(4F7.3,/,4F7.3,/,4F7.3,/,4F7.3,/,F7.3)') &
   &tb_toi(1,idir,ilon,ilat),tb_toi(2,idir,ilon,ilat),tb_toi(3,idir,ilon,ilat),tb_toi(4,idir,ilon,ilat),&
   &tb_toa(1,idir,ilon,ilat), tb_toa(2,idir,ilon,ilat), tb_toa(3,idir,ilon,ilat),tb_toa(4,idir,ilon,ilat),&
   &tb_sur(1,idir,ilon,ilat),tb_sur(2,idir,ilon,ilat),tb_sur(3,idir,ilon,ilat),tb_sur(4,idir,ilon,ilat),&
   &tb_sur0(1,idir,ilon,ilat),tb_sur0(2,idir,ilon,ilat),tb_sur0(3,idir,ilon,ilat),tb_sur0(4,idir,ilon,ilat),&
   &sss_smap_40km(idir,ilon,ilat)

   do i=0,16
      if (btest(iqc_flag(idir,ilon,ilat),i)) then
         write(10,'(I4)',advance='no') i
      endif
   enddo

   write(10,'(/2F9.3)') cellat(idir,ilon,ilat), cellon(idir,ilon,ilat)

   close(10)

   write(*,*)
end subroutine write_l2c_ascii
