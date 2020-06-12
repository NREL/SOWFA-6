# Instructions and general case notes

- Wind plant run that uses the precursor data
- Case relies on an appropriate precursor
- No turbine files included
- Adjust the refinement regions as needed. The idea of the `ground` one is that the the refinement is near the `lower` patch and consequently the `wallModelAverageType` changes
- For a crashed/restarted run, just re-submit 2_solve. No changing of the script is necessary whatsoever


# To-do

- Add change of bottom wall model depending on the refinement region placement (4th item above)
- If the flow is aligned (N, S, E, W), then cyclic can still be used on the aligned boundaries. Right now everything is being changed to TVMIO.
- Add the statistics and statistics frequency back to ABLProperties
- Include source term handling from the precursor
