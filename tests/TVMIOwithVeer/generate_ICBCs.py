#!/usr/bin/env python
import os, sys
import numpy as np
import pandas as pd
import xarray as xr

from windtools.openfoam import InputFile # for reading in setUp
from mmctools.coupling.sowfa import BoundaryCoupling
from mmctools.coupling.domain import Domain

# SETUP

casedir = '.'
output_fields = {
    'U': ('u','v','w'),
}
output_format = {
    'binary': True, # changing this to ascii or ascii+gzip will require changing
    'gzip': False,  # the format parameter in 0.original/U to "ascii"
}

wspd = 10.
angmin, angmax = -30., 30.

#==============================================================================
#
# Execution starts here
#

outdir = os.path.join(casedir,'constant','boundaryData')

# Note: The start time here is arbitrary since the simulation time doesn't
#       actually correspond to a real datetime.
sowfastart = '2000-01-01 00:00'

sowfa = InputFile(os.path.join(casedir,'setUp'))
sowfadom = Domain(
    xmin=sowfa['xMin'], xmax=sowfa['xMax'],
    ymin=sowfa['yMin'], ymax=sowfa['yMax'],
    zmin=sowfa['zMin'], zmax=sowfa['zMax'],
    nx=sowfa['nx'], ny=sowfa['ny'], nz=sowfa['nz'],
)
print(sowfadom)

# Create dataset

ds = xr.Dataset(
    coords={
        'x':[sowfadom.xmin, sowfadom.xmax],
        'y':[sowfadom.ymin, sowfadom.ymax],
        'height': sowfadom.z,
    }
)

## - Setup wind vector
ang = angmin + ds.coords['height']/sowfadom.zmax * (angmax-angmin)
ds['u'] = wspd * np.cos(np.pi/180*ang)
ds['v'] = wspd * np.sin(np.pi/180*ang)
ds['w'] = 0 * ds['u']
print('Generated boundary dataset:')
print(ds)

## - Add extra time dimension at the end (there is no time dependence)
#    Since boundary fields are constant, only need to write out one time
ds = ds.expand_dims(datetime=pd.to_datetime([sowfastart]))

# Write out boundaryData

BoundaryCoupling(outdir,
    ds.sel(x=sowfadom.xmin),
    name='west',
    dateref=sowfastart
).write(output_fields, **output_format)

BoundaryCoupling(outdir,
    ds.sel(x=sowfadom.xmax),
    name='east',
    dateref=sowfastart
).write(output_fields, **output_format)

BoundaryCoupling(outdir,
    ds.sel(y=sowfadom.ymin),
    name='south',
    dateref=sowfastart
).write(output_fields, **output_format)

BoundaryCoupling(outdir,
    ds.sel(y=sowfadom.ymax),
    name='north',
    dateref=sowfastart
).write(output_fields, **output_format)

# Write out cell-centered velocity field to include in 0.original/U
#
# Note: We assume here that the cells are predictably ordered after running
#       blockMesh. The user should verify the equivalence of the cell-center
#       locations expected by this script (output as 'assumed_blockmesh_cc.csv')
#       with output from the writeCellCenters utility, called after running
#       blockMesh. 

## - Convert domain points to cell centers
x1 = (sowfadom.x[1:] + sowfadom.x[:-1]) / 2
y1 = (sowfadom.y[1:] + sowfadom.y[:-1]) / 2
z1 = (sowfadom.z[1:] + sowfadom.z[:-1]) / 2
x,y,z = np.meshgrid(x1,y1,z1,indexing='ij')
x = x.ravel(order='F')
y = y.ravel(order='F')
z = z.ravel(order='F')

csvfile = os.path.join(casedir,'constant','assumed_blockmesh_cc.csv')
pd.DataFrame({'x':x,'y':y,'z':z}).to_csv(csvfile,index=False)
print('Wrote',csvfile)

## - Generate list of cell-centered velocities for IC
ang = angmin + z/sowfadom.zmax * (angmax-angmin)
u = wspd * np.cos(np.pi/180*ang)
v = wspd * np.sin(np.pi/180*ang)
w = np.zeros(x.shape)

U0file = os.path.join(casedir,'U0')
with open(U0file, 'w') as f:
    f.write('internalField nonuniform List<vector>\n')
    f.write('{:d}\n'.format(len(x)))
    f.write('(\n')
    for ui,vi,wi in zip(u,v,w):
        f.write('({:g} {:g} {:g})\n'.format(ui,vi,wi))
    f.write(')\n;\n')
print('Wrote',U0file)

