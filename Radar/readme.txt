This README is a brief description of the 2D CFAR algorithm implemented in radar_target_generation_and_detection.m

Implementation Steps

1. Implementation steps for the 2D CFAR process
1) Select training and guard cell size for both range and doppler.
2) Loop over all cells in the range doppler map, First calculate the noise level among all the training cells. Then set the threshold as noise level + offset. It the CUT level is greater than the threshold, mark it as a hit, otherwise, mark the CUT as miss.

2. Selection of Training, Guard cells and offset
1) Range: 5 Training Cells, 2 Guard Cells
2) Doppler: 4 Training Cells, 2 Guard Cells
3) Offset: 10 db.

3. Steps taken to suppress the non-thresholded cells at the edges
Mark all non-thresholded cells at the edges as 0.