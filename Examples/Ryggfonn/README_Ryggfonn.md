## Simulation of an Avalanche at Ryggfonn

A medium-size mixed-snow avalanche was measured at the test site Ryggfonn in Grasdalen, western Norway on 2021-04-11. The release area could be determined with some precision and the release depth estimated from laser scan data before and after the release. Similarly, the deposit depth was inferred from a drone survey of the entire avalanche area. The velocity of the avalanche was measured with Doppler radar and also from drone imagery. Impact pressure could, however, not be collected for this event. A description of the collected data with numerous images can be found in [(Gauer et al., 2022)](https://ngi.brage.unit.no/ngi-xmlui/bitstream/handle/11250/3093249/20200017-06-TN.pdf?sequence=1&isAllowed=y).

The main directory contains the following files:

* The file `Ryggfonn_2021-04-11.qgz` can be used to initialize QGIS
  for this site.
* Three run configuration files for MoT-Voellmy describing
    * no entrainment with constant friction parameters,
    * no entrainment with spatially variable friction parameters as
      proposed by RAMMS::AVALANCHE,
    * a run with the RAMMS friction parameters and the TJEM entrainment
      module enabled.

The subfolder `Input` contains raster and vector files:

* a DTM with resolution 5×5 m²,
* two shapefiles describing the main release area and the main
  release area combined with the secondary one,
* a raster file with the main release area and a release depth of 1.0 m,
* two GeoPackage files with contour lines at 100 and 10 m intervals,
* two raster files `rgf_MM_S300...asc` with the spatially distributed _μ_
  and _k_ values at 5×5 m² resolution, as proposed by RAMMS::AVALANCHE
  v.1.8.0, which can be used in MoT-Voellmy without changes, and
* `b0.asc` and `tauc.asc` containing fictitious data for the spatial
  distribution of initial snow depth (_b_<sub>0</sub>) and shear strength
  of the snow cover (_τ_<sub>c</sub>).

`b0.asc` and `tauc.asc` were created from the DTM with the simple `gawk` commands:

    gawk ' BEGIN { b00 = 0.5 ; db = 0.0003 }
           { if (NR < 7) print
             if (NR > 6) {
                 for (i=1; i<NF; i++) printf("%4.2f  ", b00 + db*$i) 
                 printf("%4.2f\n", b00 + db*$NF)
             }
           } ' dhm.asc > b0.asc

    gawk ' BEGIN { tau0 = 600.0 ; dtau = -0.2 }
           { if (NR < 7) print ;
             if (NR > 6) {
                 for (i=1; i<NF; i++) printf("%5.1f ", tau0 + dtau*$i)
                 printf("%5.1f\n", tau0 + dtau*$NF)
             }
           } ' dhm.asc > tauc.asc

`gawk` is usually installed in Linux distributions and macOS, and it can be downloaded for MS Windows from [Sourceforge](https://gnuwin32.sourceforge.net/downlinks/gawk-bin-zip.php). The commands above are to be run from the directory `Input`. The interpolation functions are __not__ based on actual measurements; the functions or their parameters can readily be modified.

Note that the parameters were __not__ set to reproduce the event of 2021-04-11 in the best possible way—the purpose here is to show how MoT-Voellmy can be run. The users are invited to play with the parameters to obtain a feeling for how the model reacts to variations. For comparison, the output files from a simulation with RAMMS:AVALANCHE v.1.8.0 (contributed by Peter Gauer, NGI) for the initial conditions and parameters corresponding to Run_02 are collected in the subdirectory `RAMMS_1.8.0/`.