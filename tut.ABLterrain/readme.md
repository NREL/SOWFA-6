# Instructions and general case notes

- First script requests more cores than `$nCores` for RAM workaround for `moveDynamicMesh` call. This is necessary if the `stl` file is large.
- For precursor-mapped runs, any `coded` boundary condition gets overwritten by `changeDictionary` call to set  timeVarying-type of BC
- `0.original/T.<type>` has custom boundary conditions. Copy that as appropriate
- Terrain `stl` file not included
- Old usage: The terrain vertical bounding box should be below `zMin`. Use `surfaceTransformPoints` to translate stl if needed. See `readme` file inside `constant/triSurface`
- New terrain handling: use `~/utilities/shiftFlatBlendToz0.sh` to shift the terrain. The script works by calling `shiftFlatBlendToz0 terrainFile.stl -1000 0`, where the numbers are the $x$ and $y$ of the inlet plane, which we want to be flush at $z=0$
- Final mesh is assumed to be in the last saved time step from the `moveDynamicMesh` execution
- `setFieldsABLDict` contains `useWallDist` as `true` by default. Adjust the `temperatureInitType` as needed
- For the mesh to conform to the terrain, the boundaries cannot be cyclic

# To-do

- if stable/unstable, be mindful of the qwall in the change dictionaries. 
- remove the cyclic patches on decomposepardict if running non-cyclic
- Make sure the restarted script will go through fine because of the time averaging `startTimeAvg`. Fix it.
- Put the default BC as inflow/outflow and update second script (set fields ABL)
- See if `$Rwall` and `$qwall` need to be included in the `changeDictionaryDict`s
- Allow different wind directions
- Allow different stability states
- Incorporate Eliot's check whether or not run was successful

