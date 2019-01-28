# Install notes for Needles

## Download the 3 following files in a directory of your choice:

-	http://repo.mouseimaging.ca/repo/DSURQE_40micron_nifti/DSURQE_40micron_average.nii
-	http://repo.mouseimaging.ca/repo/DSURQE_40micron_labels.nii
-	http://repo.mouseimaging.ca/repo/DSURQE_40micron/DSURQE_40micron_R_mapping.csv

## Clone or download the WGs repository

Within the repository, edit the `./WGs/Physiology/needles_param.json` file and set the path_atlas to the folder in which you did download the 3 files mentionned above. For example: 
	` "path_atlas": "/datadisk/BrainAtlas/ATLASES/DSURQE_40micron",`

## Setup the Matlab Path and run

For general use without worrying about the path, run `RunNeedles.m`  file in Matlab

For advanced use with a non-frozen ibllib library, set your paths to:
-	`./WGs/Physiology/`
-	the ibllib Matlab library
