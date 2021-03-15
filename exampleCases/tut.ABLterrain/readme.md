# Instructions and general case notes

- Complex-terrain wind plant run that uses the precursor data
- A manual setup of `changeDictionaryDict` is needed depending on flow direction
- Case relies on an appropriate precursor
- No turbine files included
- Adjust the refinement regions as needed. The idea of the `ground` one is that the the refinement is near the `lower` patch and consequently the `wallModelAverageType` should be modified
- First script requests more cores than `$nCores` for RAM workaround for `moveDynamicMesh` call. This is necessary if the `stl` file is large
- `0.original/T.<type>` has custom boundary conditions. Copy that as appropriate
- Terrain `stl` file not included
- Final mesh is assumed to be in the last saved time step from the `moveDynamicMesh` execution
- If mapping from precursor, the valocity needs to be scaled by calling `setFieldsABL` with the appropriate `scaleGiven` flag
- `setFieldsABLDict` contains `useWallDist` as `true` by default. Adjust the `temperatureInitType` as needed
- For restarted runs, simply re-submit the second script, without changing

