module dir_paths_module
   implicit none
   character(len=220) ::  table_dir
   character(len=220) ::  data_dir


contains
   subroutine get_dir_paths(table_dir, data_dir)
      implicit none

      character(len=220), intent(out)      ::  table_dir
      character(len=220), intent(out)      ::  data_dir

      integer(4)     ::  env_stat

      call get_environment_variable("SMAP_TABLE_DIR", table_dir, status=env_stat)
      if (env_stat == 1) then
         ! Variable does not exist so use the fallback value
         table_dir = "../tables_L2C"
      else if (env_stat == 0) then
         ! Environment variable sucessfully read
         continue
      else
         ! No environment variable support or the string is too short
         error stop "did not read 'SMAP_TABLE_DIR'"
      end if

      call get_environment_variable("SMAP_DATA_DIR", data_dir, status=env_stat)
      if (env_stat == 1) then
         ! Variable does not exist so use the fallback value
         data_dir = "../sample_data"
      else if (env_stat == 0) then
         ! Environment variable sucessfully read
         continue
      else
         ! No environment variable support or the string is too short
         error stop "did not read 'SMAP_DATA_DIR'"
      end if


   end subroutine



endmodule



