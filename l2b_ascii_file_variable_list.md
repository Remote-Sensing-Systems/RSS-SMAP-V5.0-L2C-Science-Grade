This file lists the variables in the SMAP L2B input ASCII files (in order) as well as what they represent.
The SMAP L2B input ASCII files have the naming convention `sample_l2b_38362_##.txt` where `38362` is the SMAP
orbit number and `##` represents the pixel number (an integer between 1 and 10).
Each of these values contained in a given file is for the pixel number in the L2B input ASCII file name.
Note that all listed variables from `tec` onward are ancillary data needed for the RSS SMAP V5 L2C processing.

| Variable | Description | Units |
| --- | --- | --- |
| ilon | the index of the pixel's x-grid | integer between 1 and 1560 |
| ilat | the index of the pixel's y-grid | integer between 1 and 720 |
| idir | the index of the pixel's look direction | 1=for, 2=aft |
| ifill | indicates if the Earth cell contains any resampled observations | 1=contains observations or 0=no observations |
| alpha | SMAP scan angle | degrees |
| iscan | center scan index of the target cell. Dummy argument not used in L2C processing. | N/A |
| dist | distance from Earth to target cell.  Dummy argument not used in L2C processing. | km |
| time | time of observation | seconds since 2000-01-01 00:00:00 |
| zang | orbital position of spacecraft | degrees |
| cellat | latitude of pixel | degrees N |
| cellon | longitude of pixel | degrees E |
| eia | Earth incidence angle | degrees |
| eaa | azimuthal look angle relative to North | degrees |
| pra | polarization basis rotation angle | degrees |
| sunglt | sunglint angle | degrees |
| monglt | moonglint angle | degrees |
| gallat | latitudinal (polar) angle of specular reflection ray in ECI2000 | degrees |
| gallon | longitudinal (azimuthal) angle of specular reflection ray in ECI2000 | degrees |
| frdrot | Faraday rotation | degrees |
| teclat | lat of the intersection of the ray with the ionospheric shell (450 km). Dummy argument not used in L2C processing. | degrees N |
| teclon | lon of the intersection of the ray with the ionospheric shell (450 km). Dummy argument not used in L2C processing. | degrees E |
| loss_ant | antenna loss | K |
| temp_ant | physical temperature of reflector from L1B files | K |
| sun_alpha | sun azimuth angle in spacecraft coordinate system | degrees |
| sun_beta | sun zenith angle in spacecraft coordinate system | degrees |
| ta_ant | unfiltered antenna temperature | K |
| ta_ant_filtered | RFI filtered antenna temperature | K |
| ta_sun_dir | TA of direct sun intrusion | K |
| ta_sun_ref | TA of reflected sun intrusion | K |
| ta_gal_dir | TA of direct galaxy intrusion | K |
| ta_gal_ref_tab | TA of reflected galaxy from a pre-computed table. 5 different values, 1 for each of 5 different wind speeds | K |
| ta_lnd | TA of land | K |
| fland | land fraction within footprint | value between 0 and 1 |
| gland | land fraction weighted by antenna gain | value between 0 and 1 |
| wt_sum | sum of the Backus Gilbert Optimum Interpolation weights | N/A |
| tec | vertically integrated electron content between surface and spacecraft | 1e16 m-2 |
| sss_ref | reference sea surface salinity from HYCOM | psu |
| surtep | ancillary sea surface temperature (from CMC) | K |
| winspd_anc | ancillary sea surface wind speed from CCMP NRT that is used in surface roughness correction | m/s |
| windir | ancillary wind direction relative to North (from NCEP, meteorological convention) that is used in surface roughness correction | degrees |
| absp_oxy | atmospheric oxygen absorption coefficient | N/A |
| absp_vap | atmospheric water vapor absorption coefficient | N/A |
| absp_liq | atmospheric liquid water absorption coefficient | N/A |
| tran | total atmospheric transmittance | value between 0 and 1 |
| tbup | atmospheric upwelling brightness temperature | K |
| tbdw | atmopsheric downwelling brightness temperature | K |
| lst | land surface temperature | K |
| sm | soil moisture | m3 m-3 |
| solar_flux | ancillary solar flux (from NOAA SWPC) | 1e-22 W m-2 Hz-1 |
| rain | IMERG rain rate. Resampled to SMAP resolution | mm hr-1 |
| sst_err | ancillary sea-surface temperature error (from CMC) | K |
| xobs_sat_ccmp | number of observations in ancillary CCMP pixel | N/A |
| gice_est | estimated sea ice fraction weighted by antenna gain | value bewteen 0 and 1 |
| dtb_sea_ice_corr | TB sea-ice correction applied at TB SUR 0 | K |
| dtb_sea_ice_error_est | error associated with dtb_sea_ice_corr | K |
| icezone | sea-ice contamination zone | value between 0 and 7 |
| iceflag_amsr2 | AMSR2 iceflag associated with measurement.  Used in sea-ice correction. | value between 1 and 3 |
