LOCATE - LOCally Adaptive Threshold Estimation 

Estimating the threshold adaptively for the white matter lesion probability maps obtained from BIANCA (refer to https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BIANCA/Userguide for more details)

Example training call:
LOCATE_training(train_image_directory_name);

Example testing call:
LOCATE_testing(test_image_directory_name, train_image_directory_name);

Example Leave-one-subject-out testing call:
LOCATE_LOSO_testing(train_image_directory_name);

LOCATE requires the training and test data to be prepared in a specific way for its execution. Kindly refer to the user manual for more information on data preperation and additional options on LOCATE training and testing.
