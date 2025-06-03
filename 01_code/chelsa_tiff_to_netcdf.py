conda activate CHELSA_paleo
export SINGULARITY_IMG="/home/dafcluster4/chelsa_paleo/singularity/chelsa_paleo.sif"
singularity exec "$SINGULARITY_IMG" python

import os
import re
from osgeo import gdal, gdalconst
import numpy as np
from tqdm import tqdm
import sys
sys.path.append('01_code/00_function')
from chelsa_to_netcdf.py import *


