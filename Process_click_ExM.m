%% Data analysis of Click_ExM pre and post expansion images
% 
% This code is to process data for Click_ExM. Contains the following five parts:
% (A) Pre-adjustments of pre and post ExM images:
%        1) read 8-bit grayscale pre and post images, transform to double
%        2) use the estiamted expansion factor (macroscopic measurement of gel dimensions, 4.5)
%             to stretch pre-ExM image by bicubic interpolation to match post_ExM size
% (B) Rigid registration and correct expansion factor
%        1) perform rigid registration (affine)
%        2) use the stretch parameters to correct the estimated expansion factor
% (C) Non-rigid registration
%        1) Gaussian blur mask to exclude regions with no features
%        2) Nonrigid registration
% (D) Distortion analysis
%        1) Distortion field
%        2) RMS error calculation
% (E) Save processed data
%        (i)   gray, rigid registered pre and post ['(i) compare pre and post (registered).svg']
%        (ii)  green-magneta, before and after rigid reg, show expansion factor ['(ii) rigid registration and exp fac.svg']
%        (iii) non-rigid registration and distortion field ['(iii) non-rigid registration and distortion field.svg']
%        (iv) RMSE cropped area
%        (v)  smoothed RMSE plot
%        (vi) RMSE calculated values as xlsx
%
% Attached dependencies: 
% registerImages.m [rigid registration function];
% nr_registerImages.m [non-rigid registration function];
% RMSE_cal.m [calculate RMSE];
% RMSE_MS.m [calculate mean and std of RMSE];
%
% Software requirements: MATLAB R2018a or later
%
% 05-Jun-2020
% Yujie Shi
% yujieshi@scripps.edu

%% (A) Pre-adjustments of pre and post ExM images
clc,clear
pre = imread('pre_Click_ExM_tubulin.tif');     % read pre-expansion image (8-bit TIFF required) 
post = imread('post_Click_ExM_tubulin.tif');     % read post-expansion image (8-bit TIFF required)
pre_um_per_pix = 18.52/77;     % transform pixel to um according to the scale of pre and post images
post_um_per_pix = 87.30/363;
total_corr = pre_um_per_pix/post_um_per_pix;     % correct post and pre scale
exp_fac_esti = 4.3;     % estimated stretch factor by macroscopic measurement


moving = im2double(post);     % use post-expansion image as moving
fixed_pre = im2double(pre);     % stretch pre-expansion image to prepare for registration as fixed
[pre_x, pre_y] = size(fixed_pre);
[post_x, post_y] = size(moving);
exp_fac_esti_corr = exp_fac_esti*total_corr;
[w0,h0] = meshgrid(1:1:pre_x, 1:1:pre_x);
[w1,h1] = meshgrid(1:1/exp_fac_esti_corr:pre_x, 1:1/exp_fac_esti_corr:pre_x);
fixed = interp2(w0,h0,fixed_pre,w1,h1,'cubic');     % bicubic interpolation for stretching
[fixed_x, fixed_y] = size(fixed);
figure
imshowpair(fixed,moving, 'montage');     % show pre and post images before registration


%% (B) Rigid registration and correct expansion factor
reg_final = registerImages(moving,fixed);     % registration process
reg_transformation = reg_final.Transformation.T;     % transformation matrix of rigid registration
x_stretch_corr = reg_transformation(1,1);
y_stretch_corr = reg_transformation(2,2);
stretch_corr = sqrt(x_stretch_corr * y_stretch_corr);
exp_factor = exp_fac_esti_corr/stretch_corr/total_corr;     % corrected expansion factor
AA = figure;
set(gcf,'Position', [0 0 1500 750]);
subplot(1,2,1)
imshowpair(fixed,moving);     % original alignment of fixed and moving
title({['Before registration']; ['(Estimated Expansion Factor = ', char(string(exp_fac_esti)), ')'];
    ['green: pre_E_x_M (fixed)']; ['magneta: post_E_x_M (moving)']})
subplot(1,2,2) % plot after registration
imshowpair(fixed,reg_final.RegisteredImage);     % registration result
title({['After registration'];['(Corrected Expansion Factor = ', char(string(exp_factor)), ')'];
    ['green: pre_E_x_M (fixed)'];['magneta: post_E_x_M (moving)']})
imwrite(reg_final.RegisteredImage, 'registered_post_ExM.tif')

%% (C) Non-rigid registration of registered post_ExM image
nr_moving_ori = im2double(imread('registered_post_ExM.tif'));     % rigid-reg post ExM img as moving
nr_moving_gaussian = imgaussfilt(nr_moving_ori, 2.5);     % gaussian blur of moving, sigma = 2.5
nr_moving_mask = imbinarize(nr_moving_gaussian, 'global');     % generate bg mask for moving
nr_moving_subbg = nr_moving_ori .* nr_moving_mask;     % moving sub bg
nr_fixed_gaussian = imgaussfilt(fixed, 3);     % pre Exm as fixed, gaussian blur, sigma = 3
nr_fixed_mask = imbinarize(fixed, 'global');     % generate bg mask for fixed
se = strel('disk', 4);     % dilate to connect fragmented mask with a disk structural element
nr_fixed_mask_dil = imdilate(nr_fixed_mask, se);    
nr_fixed_subbg = fixed .* nr_fixed_mask_dil;     % fixed sub bg
reg_nr_final = nr_registerImages(nr_moving_subbg, nr_fixed_subbg);     % non-rigid registration
NR_dis = figure;
set(gcf,'Position', [0 0 1500 750]);
subplot(1,3,1)
imshowpair(fixed,reg_final.RegisteredImage);     % rigid registration result
title({['Rigid registration']; ['']; ['green: pre_E_x_M (fixed)'];
    ['magneta: post_E_x_M (moving)']})
subplot(1,3,2)
imshowpair(fixed,reg_nr_final.RegisteredImage)     % non-rigid registration result
title({['Non-rigid registration']; ['']; ['green: pre_E_x_M (fixed)'];
    ['magneta: post_E_x_M (nr_ moving)']})
subplot(1,3,3)
imshowpair(nr_moving_ori,reg_nr_final.RegisteredImage);     % for distortion field plot
title({['distortion field']; ['']; ['green: post_E_x_M (before nr register)'];
    ['magneta: post_E_x_M (after nr register)']})

%% (D) Distortion analysis I
% distortion field
x = 1:20:fixed_x;     % generate mesh for distortion field
y = 1:20:fixed_y;
[im_x,im_y] = meshgrid(x,y);
im_u = reg_nr_final.DisplacementField(x,y,1);
im_v = reg_nr_final.DisplacementField(x,y,2);
hold on
quiver(im_x,im_y,-im_u, -im_v, 1, 'Color', [1 1 0]);     % plot distortion field

% generate align img with distortion field for generating area for RMSE calculation
figure
imshowpair(nr_moving_ori,reg_nr_final.RegisteredImage);     % for distortion field plot
h=getimage(gcf);
Dis_1=figure;
imshow(h,'border','tight')
hold on
quiver(im_x,im_y,-im_u, -im_v, 1, 'Color', [1 1 0]);     % plot distortion field
saveas(Dis_1,'Dis_field_wrongsize.tif','tif')
Dis_2=imread('Dis_field_wrongsize.tif');
[wrong_x, wrong_y, ~] = size(Dis_2);
Dis_field=imresize(Dis_2,  fixed_x/wrong_x, 'bicubic');     %Dis_field is an image with distortion field plotted

% RMSE
Dis_x = reg_nr_final.DisplacementField(:,:,1);
Dis_y = reg_nr_final.DisplacementField(:,:,2);

figure
I=imshow(nr_moving_subbg);
waitfor(msgbox('Please select the region desired for RMSE calculation with your mouse in the next window'));
[I2, cropping]=imcrop(I);
figure
imshow(I2);

Dis(:,:,1) = imcrop(Dis_x,cropping);
Dis(:,:,2) = imcrop(Dis_x,cropping);
mask = imcrop(nr_fixed_mask_dil, cropping);

%% (E) Save processed data I [Everything before RMSE]
mkdir processed_data

% (i) gray, rigid registered pre and post
I_rr_align = figure;
imshowpair(fixed,reg_final.RegisteredImage, 'montage');     % (line 59) registration result
title ('registered pre and post compare')
saveas(I_rr_align, './processed_data/(i) compare pre and post (registered).svg', 'svg');
close


% (ii) green-magneta, before and after rigid reg, show expansion factor
saveas(AA, './processed_data/(ii) rigid registration and exp fac.svg', 'svg');     % (line 53) compares before and after rigid registration, green-mag, exp_factor

% (iii) Non-rigid registration and distortion field
saveas(NR_dis, './processed_data/(iii) non-rigid registration and distortion field.svg', 'svg');     % (line 75) non-rigid reg, distortion field

% (iv) RMSE
crop_Dis_field = figure;
D_crop = imcrop(Dis_field, cropping);
imshow (D_crop);
title ('selected region for RMSE calculation')
saveas(crop_Dis_field , './processed_data/(iv) Cropped area for RMSE analysis (Distortion field).svg', 'svg');
% saveas(crop_Dis_field , './processed_data/(iv) Cropped area for RMSE analysis (Distortion field).tif', 'tif');
% cropped area for RMSE calculation (shown as part in distortion field (line 97) )

close all;     % generate more memory for RMSE_cal


%% (D) Distortion analysis II (RMSE Calculation)
[d,RMSE] = RMSE_cal(Dis, mask);

% RMSE calculate mean and std, transform to um and pre-expansion scale
[uni_d, RMSE_M, RMSE_S]=RMSE_MS(d,RMSE);
uni_d_u = uni_d .* post_um_per_pix ./ exp_factor;     % change pixel to um, scale to pre-exp using scale information obtained from ImageJ
RMSE_M_u = RMSE_M .* post_um_per_pix ./ exp_factor;
RMSE_S_u = RMSE_S .* post_um_per_pix ./ exp_factor;

% RMSE plot by medfilter to smooth the plot
lo = RMSE_M_u - RMSE_S_u;     % mean-std
hi = RMSE_M_u + RMSE_S_u;     % mean+std
me = RMSE_M_u;     % mean
len = length(RMSE_M_u);     % length of RMSE
lo_smooth = medfilt1(lo, 100);
hi_smooth = medfilt1(hi,100);
me_smooth = medfilt1(me,100);

% final calculated RMSE related values for plot
d_RMSE = uni_d_u(1: len-50);
lo_RMSE = lo_smooth(1:len-50);
m_RMSE = me_smooth(1:len-50);
h_RMSE = hi_smooth(1:len-50);

RMSE_plot = figure;
hold on
plot(d_RMSE, lo_RMSE, 'b--');
plot(d_RMSE, h_RMSE, 'b--');
plot(d_RMSE, m_RMSE, 'k-');
title('RMSE v.s. measurment length')
xlabel('measurement length (um), preexpansion scale');
ylabel('RMSE (um), preexpansion scale');


%% (E) Save processed data II [RMSE]

% RMSE plot
saveas(RMSE_plot, './processed_data/(v) RMSE plot.svg', 'svg');     % smoothed RMSE_plot (line 132)

% RMSE calculated values
xlswrite('./processed_data/RMSE measurement_length.xlsx', d_RMSE)
xlswrite('./processed_data/RMSE mean-std.xlsx', lo_RMSE)
xlswrite('./processed_data/RMSE mean+std.xlsx', h_RMSE)
xlswrite('./processed_data/RMSE mean.xlsx', m_RMSE)

close all
waitfor(msgbox('Processing finished. Check results in the processed_data folder.'));





