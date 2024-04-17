files contained in this repo (or available as shared file on Google drive because they exceed the 100Mb Github limit)

conus/conusjunctions.zip (NOT IN REPO - available at https://drive.google.com/file/d/1ZrIUyqb0x6_-4_SRHzLTMVY7WcPqLKrs/view?usp=sharing)
 zipped conusjunctions.txt
 ascii file containing the output of junctionextraction.c
 7 columns (x pixel location, y pixel location, contributing area, AVB slope ratio, AVB angle in degrees, BA slope ratio, BA angle in degrees) and 20,896,511 rows
conus/conusmeanjunctions.txt
 ascii file output by junctionextraction.c quantifying mean junction angle in 30 bins of slope ratio for AVB and BA properties 
conus/conusmeanjunctions2p5km.txt
 ascii file of geometric mean of junction angles at 2.5 km/pixel (1848 columns x 1060 rows) for CONUS
conus/conusnetwork.tif (NOT IN REPO - available at https://drive.google.com/file/d/1oQGhC3CW3JlttKxthR611AG-ud-a-21p/view?usp=sharing)
 tiff image of valley network extracted by junctionextraction.c
conus/conussurfgeo.zip
 zipped conussurfgeo.txt 
 raster of size 92401x58059 representing the surficial geology of CONUS in same lcc projection as DEM
 gdal projection definition: "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs"
 xmin, ymin, xmax, ymax: -2361180 -1495630 2258880 1407330
 0 - no data;
 1 - biogenic (e.g., calcareous sedimentary rock in S FL);
 2 - aeolian;
 3 - glacial;
 4 - lacustrine;
 5 - Ogallala;
 6 - Plio-Q alluvium;
 7 - Tertiary;
 8 - volcanic rock;
 9 - water
conus/conussurfgeo2p5km.txt
 ascii file of the presence/absence of late Cenozoic alluvial deposits at 2.5 km/pixel (1848 columns x 1060 rows) for CONUS
conus/junctionextraction.c
 C code that extracts junction angles for CONUS at 50 m/pixel resultion. The input file conus.txt is 43Gb (11Gb compressed) and available upon request from jdpellet@arizona.edu. 
 author compiled and tested on a 64-bit Windows machine in cygwin using "gcc -o juctionextraction.exe junctionextraction.c -lm" 
 author executed in cygwin using "./junctionextraction.exe"

myfractaltree/junctionextraction.c
 C code that extracts junction angles for a test landscape using the 1201x800 DEM myfractaltreetopo.txt. 
 author compiled and tested on a 64-bit Windows machine in cygwin using "gcc -o juctionextraction.exe junctionextraction.c -lm" 
 author executed in cygwin using "./junctionextraction.exe"
myfractaltree/myfractaltreejunctions.txt
 ascii file containing data for every junction (output of junctionextraction.c) for the test landscape
 this file has 7 columns (x pixel location, y pixel location, orientation of junction direction 1, orientation of junction direction 2, orientation of junction direction 3, angle between 1 and 3, angle between 2 and 3) and 13 rows
myfractaltree/myfractaltreearea.txt
 ascii file of contributing area raster (1201x800) output by junctionextraction.c
myfractaltree/myfractaltreemask.txt
 ascii file of mask raster (1201x800) output by junctionextraction.c (encodes tributary directions for each junction)
myfractaltree/myfractaltreetopo.txt
 ascii DEM (1201x800) input to junctionextraction.c

soaz/junctionextraction.c
 C code that extracts junction angles for a portion of southern Arizona at 30 m/pixel resultion using the 6833x5333 DEM soaz.txt. 
 author compiled and tested on a 64-bit Windows machine in cygwin using "gcc -o juctionextraction.exe junctionextraction.c -lm" 
 author executed in cygwin using "./junctionextraction.exe"
soaz/soaz.zip (NOT IN REPO - available at https://drive.google.com/file/d/1vxM_QshV6K_LHLEQDBMtKIFWZmM89k6d/view?usp=sharing)
 zipped version of soaz.txt - 30 m/pixel resolution 6833x5333 DEM for a portion of southern Arizona
soaz/soazmeanjunctions.txt
 ascii file output by junctionextraction.c quantifying mean junction angle in 30 bins of slope ratio for AVB and BA properties (plotted in Fig. 8c of paper)
soaz/soazmeanjunctions2p5km.txt
 ascii file of mean junction angles at 2.5 km/pixel (136 columns x 106 rows) for a portion of southern Arizona
soaz/soazsurfgeo2p5km.txt
 ascii file of presence/absence of late Cenozoic alluvial deposits at 2.5 km/pixel (136 columns x 106 rows) for a portion of southern Arizona
soaz/soazjunctions.txt
 ascii file containing data for every junction (output of junctionextraction.c) for a portion of southern Arizona
 this file has 7 columns (x pixel location, y pixel location, contributing area, AVB slope ratio, AVB angle in degrees, BA slope ratio, BA angle in degrees) and 106,506 rows