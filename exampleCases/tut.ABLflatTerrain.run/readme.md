# Instructions and general case notes

- Wind plant run that uses the precursor data
- For WRF internal coupling, this is the "precursor" case to use
- For WRF cases, specify the correct `qwall` BC. The syntaxes are commented out on `0.forWRFcoupled/qwall.wrf`
- A manual setup of `changeDictionaryDict` is needed depending on flow direction
- Unless WRF internal coupled, this case relies on an appropriate precursor
- No turbine files included
- Adjust the refinement regions as needed. The idea of the `ground` one is that the the refinement is near the `lower` patch and consequently the `wallModelAverageType` should be modified
- For a crashed/restarted run, just re-submit 2_solve. No changing of the script is necessary whatsoever
- Reference pressure handling has changed (Feb/2021)
