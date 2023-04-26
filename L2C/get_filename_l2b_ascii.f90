subroutine get_filename_l2B_ascii(iorbit, ipixel, filename_l2b)
use dir_paths_module
implicit none

integer(4), intent(in)              ::  iorbit
integer(4), intent(in)              ::  ipixel
character(len=250), intent(out)     ::  filename_l2b

character(len=220)                  ::  pathname
character(len=100)                  ::  str


call get_dir_paths(table_dir, data_dir)
pathname = data_dir

write(str,1001) iorbit, ipixel
1001 format('sample_l2b_',i5.5,'_',i2.2,'.txt')
 
filename_L2B = trim(pathname) // '/' // trim(str)

return
end subroutine get_filename_l2B_ascii
