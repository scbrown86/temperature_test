## quick test of CHELSA_paleo 

Test on a small subset of TraCE-21ka data to try and recreate some previous analysis. 

There are a few input files that need to be generated between 1900 and 1989 for the test to run

#### INPUT DATA - CLIMATE DATA
files to be stored in a subdirectory /clim
- huss.nc : a netCDF file containing relative humidity at the surface of n timesteps
- pr.nc : a netCDF file containing precipitation rate at the surface of n timesteps
- ta_high.nc : a netCDF file containing air temperatures at the higher (elevation) pressure level used for the lapse rate calculation (e.g. 600.5 hPa [z=20]) of n timesteps
- ta_low.nc : a netCDF file containing air temperatures at the lower (elevation) pressure level used for the lapse rate calculation (e.g. 992.5 hPa [z=26]) of n timesteps
- tasmax.nc : a netCDF file containing daily maximum near-surface air temperature of n timesteps
- tasmin.nc : a netCDF file containing daily minimum near-surface air temperature of n timesteps
- tas.nc : a netCDF file containing daily mean near-surface air temperature of n timesteps
- uwind.nc : a netCDF file containing the zonal wind component (u) of n timesteps
- vwind.nc : a netCDF file containing the meridional wind component (v) of n timesteps
- zg_high.nc : a netCDF file containing geopotential height (in meters) at the higher (elevation) pressure level used for the lapse rate calculation (e.g. 600.5 hPa [z=20]) of n timesteps
- zg_low.nc : a netCDF file containing geopotential height (in meters) at the lower (elevation) pressure level used for the lapse rate calculation (e.g. 992.5 hPa [z=26]) of n timesteps

#### INPUT DATA - OROGRAPHIC DATA
files to be stored in a subdirectory /orog
- oro.nc : a netCDF file containing the orography at the coarse (GCM) resolution of n timesteps (modified to work with a single timestep) (i.e. TraCE-21 input orography)
- oro_high.nc : a netCDF file containing the orography at the high (target) resolution of n timesteps (modified to work with a single timestep)

#### INPUT DATA - STATIC DATA
files to be stored in a subdirectory /static
-merc_template.nc : a netCDF file containing the orography at high (target) resolution in World Mercator projection (here we use a 5-km grid)

`EPSG:3395`

`Proj4 string = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'`

#### Notes:
CHELSA_paleo has been sourced from https://gitlabext.wsl.ch/karger/chelsa_paleo and modified to only require a single (i.e. static) mask regardless of the number of timesteps being simulated. 

```python
# ds1 = import_ncdf(self.INPUT + 'orog/oro_high.nc').Get_Grid(self.time)
ds1 = import_ncdf(self.INPUT + 'orog/oro_high.nc').Get_Grid(0) #<<< use first time step only
```

It has also been modified to run using a maximum of 2 threads per step which permits us to run 12 steps (i.e. 1 year) in parallel when processing.

```python
def Load_Tool_Libraries(Verbose):
    saga_api.SG_UI_Msg_Lock(True)
    if os.name == 'nt':  # Windows
        os.environ['PATH'] = os.environ['PATH'] + ';' + os.environ['SAGA_32'] + '/dll'
        saga_api.SG_Get_Tool_Library_Manager().Add_Directory(os.environ['SAGA_32'] + '/tools', False)
    else:  # Linux
        saga_api.SG_Get_Tool_Library_Manager().Add_Directory('/usr/local/lib/saga/',
                                                             False)  # Or set the Tool directory like this!
    saga_api.SG_UI_Msg_Lock(False)
   
    saga_api.SG_OMP_Set_Max_Num_Threads(2)  #<<< cores hardcoded to 2

    if Verbose == True:
        print 'Python - Version ' + sys.version
        print saga_api.SAGA_API_Get_Version()        
        print 'number of maximum threads used: ' + str(saga_api.SG_OMP_Get_Max_Num_Threads())
        print 'number of loaded libraries: ' + str(saga_api.SG_Get_Tool_Library_Manager().Get_Count())
        print

    return saga_api.SG_Get_Tool_Library_Manager().Get_Count()
```

CHELSA_paleo has been installed in a conda environment named `CHELSA_paleo` which will need to be activated before running the test.

_most_ of the pre-processing occurs in `R` but there are some steps performed in bash scripts which require the `nco_stable` environment to be activated.
