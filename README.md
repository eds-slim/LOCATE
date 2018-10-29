LOCATE - LOCally Adaptive Threshold Estimation (https://www.biorxiv.org/content/early/2018/10/08/437608)

Estimating the threshold adaptively for the white matter lesion probability maps obtained from BIANCA (refer to https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BIANCA/Userguide for more details)

To run LOCATE on your machine, clone the git repository in your working directory. To do this, open your terminal and type the following command in your working directory:

git clone https://git.fmrib.ox.ac.uk/vaanathi/LOCATE-BIANCA

Please ensure that you have FSL installed in your machine and the following environment variables/paths are correctly set in your system by typing the following commands in your terminal and check if they provide similar output:

Example training call:
LOCATE_training(train_image_directory_name);

Example testing call:
LOCATE_testing(test_image_directory_name, train_image_directory_name);

Example Leave-one-subject-out testing call:
LOCATE_LOSO_testing(train_image_directory_name);

LOCATE requires the training and test data to be prepared in a
specific way for its execution. Kindly refer to the user manual
(LOCATE_User_Manual_V1.1_20052018) for more information on data preperation and additional options on LOCATE training and testing.


