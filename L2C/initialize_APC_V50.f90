! read in external APC table
! Increased A_II from 1.0929 (V3.0, V4.0, pre-launch AP) to 1.1046
! that is about half way back from V2.0
! this is because the V4 claibration + Piepmeiers V5 L1B in the Amazon comes out about 3 K too low.

subroutine initialize_APC_V5
   use external_files_l2c_module
   use l2c_module_smap
   use dir_paths_module
   use MATRIX_ROUTINES
   implicit none

   integer(4), save        ::  istart=1, iu=3
   integer(4)              ::  i, j, ising

   if (istart==1) then

      istart=0
      ! read in APC
      call get_dir_paths(table_dir, data_dir)
      open(unit=iu,file=trim(table_dir) // '/' // apc_file,form='formatted',action='read',status='old')
      read(iu,*)
      do i=1,4
         read(iu,*) j, amat_IQ(i,1:4)
         if (i/=j) stop ' error reading APC matrix'
      enddo
      close(iu)

      call INVERT_MAT(4, AMAT_IQ,  AMAT_INV_IQ,ISING)
      if (ising /=0 ) stop ' singular APC matrix'

701   format(1x,i1,5x,4(f8.4,1x))


      write(6,*) ' A Earth Integration IQ basis'
      do i=1,4
         write(6,701) i, amat_IQ(i,1:4)
      enddo
      write(6,*)

      write(6,*) ' A inv Earth Integration IQ basis'
      do i=1,4
         write(6,701) i, amat_inv_IQ(i,1:4)
      enddo
      write(6,*)


   endif ! istart


   return
end subroutine initialize_APC_V5
