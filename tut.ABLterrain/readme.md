# Instructions and general case notes

- First script requests more cores than `$nCores` for RAM workaround for `moveDynamicMesh` call
- For restarted runs: copy script `2_solve_startAt0` into `2_solve_startAt<time> and modify `startTime`; similarly for `controlDict`
- Uniform-inflow cyclic and inflow-outflow boundary conditions are not realistic
- `0.original/T` has custom strong inversion boundary conditions
- Terrain `stl` file not included
- Final mesh is assumed to be in the last saved time step from the `moveDynamicMesh` execution.


# To-do

- Allow different wind directions
- Allow different stability states
- Add WRF boundary condition
- Finish making the `bc` variable general
- Incorporate Eliot's check whether or not run was successful

