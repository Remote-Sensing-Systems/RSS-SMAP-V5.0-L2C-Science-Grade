This file lists the variables in the SMAP L2C output ASCII files (in order) as well as what they represent.
The SMAP L2C output ASCII files have the naming convention `sample_l2c_38362_##_output_rss.txt` where `38362` is the SMAP
orbit number and `##` represents the pixel number (an integer between 1 and 10).
Each of these values contained in a given file is for the pixel number in the L2C input ASCII file name.

| Variable | Description | Units |
| --- | --- | --- |
| tb_toi | Brightness temperature at the top of the ionosphere.  There are 4 values for V-pol, H-pol, 3rd Stokes, and 4th Stokes, respectively, in the first row | K |
| tb_toa | Brightness temperature at top of atmosphere BEFORE applying land correction.  There are 4 values for V-pol, H-pol, 3rd Stokes, and 4th Stokes, respectively, in the second row | K |
| tb_sur | Brightness temperature of rough ocean surface BEFORE applying roughness correction.  There are 4 values for V-pol, H-pol, 3rd Stokes, and 4th Stokes, respectively, in the third row | K |
| tb_sur0 | Brightness temperature of flat ocean surface AFTER applying roughness correction. There are 4 values for V-pol, H-pol, 3rd Stokes, and 4th Stokes, respectively, in the 4th row | K |
| sss_smap_40km | SMAP sea surface salinity at original 40km resolution | psu |
| iqc_flag | quality flags associated with the pixel measurement.  All are located in the 6th row.  For quality flag definitions please see the [RSS SMAP V5 ATBD](https://data.remss.com/smap/SSS/V05.0/documents/SMAP_NASA_RSS_Salinity_Release_V5.0.pdf) | value between 0 and 16 |
| cellat | latitude of the pixel (first column, 7th row) | degrees N |
| cellon | longitude of the pixel (second column, 7th row) | degrees E |
