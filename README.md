# Click-ExM data process and example
Matlab code and example files for manuscript *Click-ExM enables expansion microscopy for all biomolecules, by De-en Sun, et al*.

Last edited by Yujie Shi (yujieshi@scripps.edu) Oct-06-2020.

## Introduction
The script `Process_click_EXM.m` in this repo is created for expansion microscopy data processing. This includes rigid registration, calculation of expansion factor, 
B-spline non-rigid registration, calculation of distortion field and RMS error. For detailed explanation of the code, see the annotations inside the script. 
All dependencies (`registerImages.m`, `nr_registerImages.m`, `RMSE_cal.m`, and `RMSE_MS.m`) are included in this repo. The code was written and used under *Matlab R2018a*.

## Description of the method
#### Rigid registration and expansion factor calculation
The pre-expansion images were first stretched by a scaling factor *F<sub>esti</sub>* of 3.8~4.5 with bicubic interpolation to approximately 
match the scale of the post-expansion images. 
This estimated expansion factor was decided based on the macroscopic measurement of gel size. 
The stretched pre-expansion image was used as the ‘fixed image’ for rigid registration, and the post-expansion images as the ‘moving image’. 
Rigid registration was performed using the monomodal intensity-based registration with `imregtform` function. 
Affine transformation was used in the function to correct the shear deformation of the gel. 
The *X-* and *Y-* scaling elements (*x<sub>corr</sub>* and *y<sub>corr</sub>*) in the transformation matrix calculated by the imregtform function were used to correct the rough 
scaling factor (*F<sub>esti</sub>*). 
The final corrected expansion factor (*F<sub>corr</sub>*) was therefore calculated as 
*F<sub>corr</sub>* = *F<sub>esti</sub>*&times;&radic;(*x<sub>corr</sub>*&times;*y<sub>corr</sub>*)
#### Non-rigid B-spline registration and distortion analysis
The registered pre- and post-expansion images were then subjected to B-spline non-rigid registration. 
To exclude regions with no features, masks of the rigid registered pre- and post-expansion images were generated by Gaussian blur to suppress the background. 
The displacement field and B-spline registered images were then acquired with the `imregdemons` function, 
and the distortion vector field was visualized with the `quiver` function. 
RMSE was quantified by calculating the difference of distance between each pair of matching features before and after B-spline registration, 
and plotted as a function of the distance between the matching features.

## Instruction for use
1. Download all files into a same folder.
2. Run scirpt `Process_click_EXM.m` in Matlab. This will process the included example files `pre_Click_ExM_tubulin.tif` and `post_Click_ExM_tubulin.tif` included in this folder.
3. The rigid registration and non-rigid registration will finish automatically.
4. Select a small region to calculate RMS error when prompted.
5. Wait for RMSE calculation (may take a long time if a large region is selected). A "Process finished" message box will show up when all calculations are finished.
6. All results will be written into a newly created folder called `processed_data`

## Explanation of expected output
##### There will be five `.svg` files and four `.xlsx` files in the newly created `precessed_data` folder:
##### The `.svg` files shows registration results:
* `(i) compare pre and post (registered).svg` shows the rigid-registered post-ExM image and the original pre-ExM image;
* `(ii) rigid registration and exp fac.svg` merges the rigid-registered post-ExM image with the pre-ExM image; the estimated and corrected expansion factors are shown in the title;
* `(iii) non-rigid registration and distortion field.svg` shows the post-ExM image before and after non-rigid registration and the distortion field;
* `(iv) Cropped area for RMSE analysis (Distortion field).svg` shows the cropped area for RMSE analysis selected manually by the user;
* `(v) RMSE plot.svg` shows the RMSE plot, with mean value &plusmn; standard deviation.
##### The `.xlsx` files saves RMSE related values, which can be used to make RMSE plot using other softwares
- `RMSE mean.xlsx` is the calculated mean RMSE value
- `RMSE mean+std.xlsx` is the calculated mean RMSE value + standard deviation
- `RMSE mean-std.xlsx` is the calculated mean RMSE value - standard deviation
- `RMSE measurement_length.xlsx` is the measurement lengths correlated to the mean RMSE values


