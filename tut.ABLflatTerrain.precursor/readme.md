# Instructions and general case notes

- For precursor runs, the BC are cyclic and the flowfield is initialized through `setFieldsABL`.
- For restarted runs: copy script `2_solve_startAt0` into `2_solve_startAt<time>` and modify `startTime`; similarly for `controlDict`
- Separate run for saving boundary data incorporated into the same script to avoid having another job waiting in the queue


# To-do

- `boundaryData` file that needs to used is manual as of now
- Put (0,0) in the middle. Adjust terrainSlices accordingly
- Restarted runs may not work as currently is
- Incorporate Eliot's check whether or not run was successful

