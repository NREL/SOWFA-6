# Instructions and general case notes

- Precursor script to generate initial flow fields to be used in more complex cases. The precursor is executed on a flat domain and the domain size and typical times are dependent on the stability state 
- For precursor runs, the BC are cyclic and the flowfield is initialized through `setFieldsABL`.
- For a crashed/restarted run, just re-submit 2_solve. No changing of the script is necessary whatsoever
- At the end, run `foamOK` to quickly see if the runs were ok


# To-do

- Put (0,0) in the middle. Adjust terrainSlices accordingly
- Add the statistics and statistics frequency back to ABLProperties

