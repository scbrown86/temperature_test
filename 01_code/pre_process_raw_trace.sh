#!/bin/bash

conda activate nco_stable

cd /home/dafcluster4/Desktop/TraCE_Data/

ncremap -g /home/dafcluster4/Documents/GitHub/temperature_test/02_data/00_maps/sahul_0.5r.nc -G latlon=114,110#snwe=-50.0,7.0,105.0,160.0 # 0.5x0.5 Sahul grid

# generate a map file to downscale the non temp and precip variables
ncremap -a neareststod -V Z3 --preserve=mean -R '--rgn_dst --rnr_thr=0.0' -g /home/dafcluster4/Documents/GitHub/temperature_test/02_data/00_maps/sahul_0.5r.nc -s ./raw/monthly/others/trace.36.400BP-1990CE.cam2.h0.Z3.2160101-2204012.nc -m /home/dafcluster4/Documents/GitHub/temperature_test/02_data/00_maps/trace_to_sahul_0.5r_neareststod.nc -o /home/dafcluster4/Documents/GitHub/temperature_test/02_data/00_maps/temp_output.nc

# remap TraCE21 data to 0.5 degree grid
output_dir="/home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/"
map_location="/home/dafcluster4/Documents/GitHub/temperature_test/02_data/00_maps/trace_to_sahul_0.5r_neareststod.nc"

for file in ./raw/monthly/others/*; do
    echo "$file"
    var=$(echo $file | cut -f 7 -d ".") # get variable from filename
    echo -e $var
    oname="$(basename $file .nc)"
    outname="${output_dir}${oname}.Sahul.nc"
    ncremap -v $var -m $map_location -i $file -o $outname
    infile=$outname
    unset oname outname
    oname="$(basename $infile .nc)"
    outname="${output_dir}${oname}.1900_1989CE.nc"
    # First step is 1900-01-16, last step is 1989-12-16
    if [[ "$var" == "T" || "$var" == "U" || "$var" == "V" || "$var" == "Z3" ]]; then
        cdo -setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -seltimestep,4201/5280 -sellevidx,20,26 "$infile" "$outname"
    elif [[ "$var" == "RELHUM" ]]; then
        cdo --reduce_dim -setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -seltimestep,4201/5280 -sellevidx,26 "$infile" "$outname"
    else
        cdo -setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -seltimestep,4201/5280 "$infile" "$outname"
    fi
    unset mapfile oname outname
done

cd "/home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/"

rm -rf *.Sahul.nc

cd /home/dafcluster4/Desktop/TraCE_Data/

# convert PRECC and PRECL (m/s) to volume (kg/m2/s), assuming 1000 kg/m
# Precipitation rate (kg/m²/s) = Precipitation rate (m/s) × density

cdo chunit,'m/s','kg/m2/s' -mulc,1000 -setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -seltimestep,4201/5280 -sellonlatbox,105.0,161.25,-50.0,10.0 -selname,PRECC raw/monthly/pr/trace.36.400BP-1990CE.cam2.h0.PRECC.2160101-2204012.nc /home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/trace.36.400BP-1990CE.cam2.h0.PRECC.2160101-2204012.Sahul.1900_1989CE.nc

cdo chunit,'m/s','kg/m2/s' -mulc,1000 -setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -seltimestep,4201/5280 -sellonlatbox,105.0,161.25,-50.0,10.0 -selname,PRECL raw/monthly/pr/trace.36.400BP-1990CE.cam2.h0.PRECL.2160101-2204012.nc /home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/trace.36.400BP-1990CE.cam2.h0.PRECL.2160101-2204012.Sahul.1900_1989CE.nc

cdo setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -sellonlatbox,105.0,161.25,-50.0,10.0 -selname,TSMN -seltimestep,4201/5280 ./raw/monthly/ts/trace.36.400BP-1990CE.cam2.h0.TS.2160101-2204012.nc /home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/trace.36.400BP-1990CE.cam2.h0.TS.2160101-2204012.Sahul.1900_1989CE.nc

cdo setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -sellonlatbox,105.0,161.25,-50.0,10.0 -selname,TSMN -seltimestep,4201/5280 ./raw/monthly/ts/trace.36.400BP-1990CE.cam2.h0.TSMN.2160101-2204012.nc /home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/trace.36.400BP-1990CE.cam2.h0.TSMN.2160101-2204012.Sahul.1900_1989CE.nc

cdo setreftime,1900-01-16,,1month -settaxis,1900-01-16,,1month -setcalendar,365_day -sellonlatbox,105.0,161.25,-50.0,10.0 -selname,TSMX -seltimestep,4201/5280 ./raw/monthly/ts/trace.36.400BP-1990CE.cam2.h0.TSMX.2160101-2204012.nc /home/dafcluster4/Documents/GitHub/temperature_test/02_data/01_inputs/trace.36.400BP-1990CE.cam2.h0.TSMX.2160101-2204012.Sahul.1900_1989CE.nc
