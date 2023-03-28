project('smap-l2c-science-grade',
        'fortran',
        version : '0.1',
        default_options : [
          'warning_level=0',
          # TODO: enable these two eventually
          # 'warning_level=3',
          # 'fortran_std=f2008',
        ],
        meson_version : '>=0.56')

if meson.get_compiler('fortran').get_argument_syntax() == 'gcc'
  add_project_arguments(
    # This is just a safety check
    '-fimplicit-none',
    language : 'fortran')

  if get_option('debug')
    add_project_arguments('-fcheck=all', language : 'fortran')
  endif
else
  warning('Not adding any compiler options')
endif

executable('MAKE_SMAP_L2C_V50_ascii',
           sources : [
             'L2C/adjust_tagal_ref.f90',
             'L2C/allocate_L2C_arrays_V50.f90',
             'L2C/check_orbit.f90',
             'L2C/climatology_icemask.f90',
             'L2C/correct_cal_drift_V30.f90',
             'L2C/correct_reflector_emissivity_V30.f90',
             'L2C/create_l2_qcflag_V50.f90',
             'L2C/decode_l2b_V50.f90',
             'L2C/dir_paths_module.f90',
             'L2C/external_files_L2C_module.f90',
             'L2C/fd_date_2000.f90',
             'L2C/fd_delta_temp_refl_V3.f90',
             'L2C/FDICE4.f90',
             'L2C/fd_sss_V50.f90',
             'L2C/fd_ta_earth_V50.f90',
             'L2C/fd_ta_expected_V50.f90',
             'L2C/fd_tb_sur_sic_V50.f90',
             'L2C/fd_tb_toa_lc_V50.f90',
             'L2C/find_dtb_bias_V50.f90',
             'L2C/find_ta_gal_refl.f90',
             'L2C/get_filename_l2b_ascii.f90',
             'L2C/initialize_APC_V50.f90',
             'L2C/l2c_module_smap_V50.f90',
             'L2C/land_corr_step2.f90',
             'L2C/MAKE_SMAP_L2C_V50_ascii.f90',
             'L2C/MATRIX.f90',
             'L2C/meissner_wentz_dielectric.f90',
             'L2C/openbig.f90',
             'L2C/read_l2b_ascii.f90',
             'L2C/SMAP_ROUGHNESS_GMF_V3B_module.f90',
             'L2C/sss_module.f90',
             'L2C/stokes_converters.f90',
             'L2C/sun_qc_flag.f90',
             'L2C/write_l2c_ascii.f90',
           ])
