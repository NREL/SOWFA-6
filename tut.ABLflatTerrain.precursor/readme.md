# Instructions and general case notes

- For precursor runs, the BC are cyclic and the flowfield is initialized through `setFieldsABL`.
- For restarted runs: copy script `2_solve_startAt0` into `2_solve_startAt<time>` and modify `startTime`; similarly for `controlDict`
- Separate run for saving boundary data incorporated into the same script to avoid having another job waiting in the queue


# To-do

- add a preprocessing step of copying boundaryData from, say, 20000.15432 to 20000
- Make the setup mesh general (using the four blocks). It makes it easier for mapping flat-terrain into complex-terrain cases. Update the blockmeshdict with the appropriate patches, and not cyclic. Test it.
- Put (0,0) in the middle. Adjust terrainSlices accordingly
- Restarted runs may not work as currently is
- Incorporate Eliot's check whether or not run was successful

