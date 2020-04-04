# TVMIO with Veer

This small test case tests the timeVaryingMappedInletOutlet (TVMIO) boundary
condition, simulating constant wind speed with linearly varying veer. The west
and east boundaries are pure inflow and outflow, respectively, whereas the 
south and north boundaries are both 50/50% inflow/outflow; upper and lower
boundaries are slip walls. For this case, a general setup will be used with
TVMIO specified on all horizontal boundaries. The initial velocity field
corresponds to the boundary solution that is both stationary and homogeneous.
There is neutral stratification with negligible Coriolis force.

## Dependencies
These instructions assume starting from scratch a clean environment. If you
have an existing conda environment that you're happy with, feel free to skip
step 5.

1. `cd /path/to/local/python/libraries` : Where the mmctools library will live
2. ``export PYTHONPATH="`pwd`:$PYTHONPATH"`` : Add location to Python path
3. `git clone https://github.com/a2e-mmc/mmctools.git` : Get mmctools library
4. `cd mmctools && checkout dev` : Switch to the devlopment branch
5. `conda create -n mmc python=3.7` : Create a new environment (optional)
6. `conda activate mmc` : Switch to the new (or an existing) environment
7. `conda update --file environment.yaml` : Install dependencies

## Running the case
This is intended to be a very small test case, to be run quickly on 4 cores.

1. `./generate_ICBCs.py`:

    - Create reference inflow data in `<casedir>/constant/boundaryData`
    - Create initial velocity field in `<casedir>/U0`
    - Output expected cell-center locations for verification in
      `<casedir>/constant/assumed_blockmesh_cc.csv`

2. `./runscript.preprocess.py`:

    - Generates `blockMesh`
    - Runs `writeCellCenters` and verifies the output against the output from
      `generate_ICBCs.py` (this is an extra step in this case)
    - Renumbers the mesh for computational efficiency--note that after this
      the cell-center verification will fail--and then decomposes/checks the
      resulting mesh

3. `sbatch runscript.solve.1`

    - Volume output will be saved 
