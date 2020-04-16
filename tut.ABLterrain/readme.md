# Instructions and general case notes

- First script requests more cores than `$nCores` for RAM workaround for `moveDynamicMesh` call. This is necessary if the `stl` file is large.
- For restarted runs: copy script `2_solve_startAt0` into `2_solve_startAt<time> and modify `startTime`; similarly for `controlDict`
- Uniform-inflow cyclic and inflow-outflow boundary conditions are not realistic
- `0.original/T` has custom strong inversion boundary conditions
- Terrain `stl` file not included
- The terrain sure vertical bounding box should be below `zMin`. Use `surfaceTransformPoints` to translate stl if needed
- Final mesh is assumed to be in the last saved time step from the `moveDynamicMesh` execution.


# To-do

- Put the default BC as inflow/outflow and update second script (set fields ABL)
- Allow different wind directions
- Allow different stability states
- Add WRF boundary condition
- Finish making the `bc` variable general
- Incorporate Eliot's check whether or not run was successful

