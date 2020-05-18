# Instructions and general case notes

- Wind plant run that uses the precursor data
- Case relies on an appropriate precursor
- No turbine files included
- Adjust the refinement regions as needed. The idea of the `ground` one is that the the refinement is near the `lower` patch and consequently the `wallModelAverageType` changes


# To-do

- Add the split run to save averages (could start averages from the beginning which is fine for now, but with turbines this is not ideal)
- Change `changeDictionrary` calls to `foamDictionary`. `changeDictionary` is deprecated. See use in post #18 of https://www.cfd-online.com/Forums/openfoam-solving/114435-restarting-simulations-openfoam-updated-boundary-conditions.html
- Incorporate Eliot's check whether or not run was successful

# bug

- changedictionary doesn't seem to be working. The blockmesh dict needs to have patch as opposed to cyclic?
