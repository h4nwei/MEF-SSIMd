clc,clear,close all;

%% change direction
prev_dir = pwd; file_dir = fileparts(mfilename('fullpath')); cd(file_dir);
addpath(genpath(pwd));

%% model calculation
Q = zeros(2,1); % test two fused images
imgSeqColor = uint8(load_images('./images/horse',1));
imgSeqColor = uint8(reorderByLum(imgSeqColor));
[s1, s2, s3, s4] = size(imgSeqColor);
imgSeq = zeros(s1, s2, s4);
for i = 1:s4
    imgSeq(:, :, i) =  rgb2gray( squeeze( imgSeqColor(:,:,:,i) ) ); % color to gray conversion
end

fI1 = imread('./images/Horse_Shutao.jpg');
fI1 = double(rgb2gray(fI1));
[Q(1), QMap1] = mef_ms_ssim_d(imgSeq, fI1);

fI2= imread('./images/Horse_SPDMEF.jpg');
fI2 = double(rgb2gray(fI2));
[Q(2),  QMap2] = mef_ms_ssim_d(imgSeq, fI2);

figure;
subplot(2,2,1), imshow(fI1/255), title('Li12');
subplot(2,2,2), imshow(QMap1), title('quality map');
subplot(2,2,3), imshow(fI2/255), title('SPD-MEF');
subplot(2,2,4), imshow(QMap2), title('quality map');
