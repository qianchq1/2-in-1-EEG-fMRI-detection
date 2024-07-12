# 2-in-1-EEG-fMRI-detection
Reconstruct oscillation signal from 2-in-1 detector into separate EEG and fMRI signals
(1) Download the matlab scripts and put them into the same folder as the data directory. Each data directory has a numerical number.
(2) Run Sup6_Par2.m script in Matlab environment to obtain time-dependent fMRI image series and raw EEG signals. 
(3) Run PlotReconEEG.m in Matlab environment to clean up the baseline, obtaining EEG signals with refined baseline. The final result is saved as EEG.mat and uploaded to https://doi.org/10.6084/m9.figshare.26082115 
(4) Run ReadOsc_EPIs_WriteNifti.m in Matlab environment to obtain an image that includes the whole field of view. The image is saved as EPI_WholeImg.nii and uploaded to https://doi.org/10.6084/m9.figshare.26082115
(5) Run NiftiConvert2Cut.m in Matlab environment to cut a square the FOV and create an image (EPI_OneThird.nii) that can be analyzed by AFNI. The final image is saved as EPI_OneThird.nii and uploaded to https://doi.org/10.6084/m9.figshare.26082115
