# Remote Sensing Systems Science Grade SMAP V5 L2C Processing Algorithm

This document was written to illustrate how users can run the Remote Sensing Systems Soil Moisture Active Passive (SMAP) Version 5 (V5) Level 2C (L2C) science-grade processing algorithm.  The science-grade processing algorithm is meant to be a user-friendly simplification of the full RSS SMAP V5 L2C processing code.

This code is written in Fortran and is set up to be complied using the `gfortran` compiler.  It inputs L2B ASCII files located in the `sample_data` subfolder and will output ASCII files containing relevant variables from the L2C processing.  Information on the contents of the input and output ASCII files can be found in `l2b_ascii_file_variable_list.md` and `l2c_ascii_variable_list.md`, respectively.

Also contained in the `sample_data` subfolder are example output L2C ASCII files (these end in `output_rss.txt`).  Users should first check to see that they can replicate these ASCII output files by following the instructions below.  Once done, users may edit and adapt the science-grade L2C processing code to their own needs.  

## Compiling the code

All source files are in the `L2C` folder, where `MAKE_SMAP_L2C_V50_ascii.f90` is
the main program and it depends on the remaining `.f90` files included.

[Meson](https://mesonbuild.com/) is used to build the code. Using a build
directory named `build`:

```bash
meson setup build .
meson compile -C build
```

The output executable is `build/MAKE_SMAP_L2C_V50_ascii`.

## Setting the directories
At present, there are two folders containing external files that the routines and subroutines of the science-grade L2C code require in order to run. These are the `tables_L2C` and `sample_data` subfolders.  If users decide to keep the directory structure the same as it is on GitHub, then they may skip this step and move on to running the code.  However, if users wish to place the data located in `tables_L2C` and `sample_data` in different locations, they must explicitly tell the code these locations before it can be run.  This is done via the following terminal commands:

`export SMAP_TABLE_DIR="TABLE_DIRNAME"`
`export SMAP_DATA_DIR="DATA_DIRNAME"`

Where `TABLE_DIRNAME` points to the location of the contents in the `tables_L2C` folder and `DATA_DIRNAME` points to the location of the contents of the `sample_data` folder.  Note that, if the users decide to change the directory/folder structure, all data from the `tables_L2C` folder must remain in a folder together and all data from the `sample_data` folder must remain in a folder together. 

## Running the code

Once the code has been compiled and the directories have been defined, the user will be able to run the code with the following command in their terminal:

`./MAKE_SMAP_L2C_V50_ascii 38362 38362 ##`

In this command, `38362` represents the start and end SMAP orbit.  This should remain unchanged for the purposes of attempting to replicate the output ASCII files located in the `sample_data` folder.  `##` represents the pixel number that the user wishes to run the code for.  There are 10 different 0.25x0.25 degree Earth grid pixels that can be used in this demo, each assigned a value between `1` and `10`.  Therefore, `##` should be an integer between `1` and `10`.  Information on these 10 pixels is given in `pixel_descriptions.md`.

After running this code, the output L2C ASCII files will be written to the `sample_data` folder and have the file name `sample_l2c_38362_##_output.txt` where `##` is, again, the pixel number.  The output file created by the user should be compared to the `.txt` file located in the `sample_data` folder with the name `sample_l2c_38362_##_output_rss.txt`.

