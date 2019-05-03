clc,clear,close all;

%% change direction
prev_dir = pwd; file_dir = fileparts(mfilename('fullpath')); cd(file_dir);
addpath(genpath(pwd));

imgSeqFiles = {'0'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; '13'; '14'; '15'; '16';'17';'18';'19'};
load dImage;
imgPath = '/home/leshier/hanwei/TIP_code_clean/data/DeghostingDB/';

Q = zeros(20,9,2);
result=zeros(20,9);

q_map=cell(20,9);


for i = 1:size(imgSeqFiles,1)
     display(i);
    filePath = sprintf('/home/leshier/hanwei/TIP_code_clean/data/hdrSequence_png/%s', imgSeqFiles{i});
    imgSeqColorNor = uint8(load_images(filePath,1));
    imgSeqColorNor = uint8(reorderByLum(imgSeqColorNor));
    [~, ~,s3, s4] = size(imgSeqColorNor);
    for j = 1:9
        display(j);
        fI = imread(sprintf('%s/%s',imgPath, dImage{i,j}));
        f = double(rgb2gray(fI));
        [s1,s2]=size(f);
        imgSeqColor=zeros(s1,s2,s3,s4);
        imgSeq = zeros(s1, s2, s4);
        
        for l=1:s4
            temp=imresize(imgSeqColorNor(:,:,:,l),[s1,s2], 'bicubic');
            imgSeqColor(:,:,:,l)=temp;
            imgSeq(:,:,l)=rgb2gray( squeeze(temp ));
        end
        tic
        [oQ, qMap] = mef_ms_ssim_d(imgSeq, f);
        toc
        Q(i, j, 1) = oQ;
         
    end
end