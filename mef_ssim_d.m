function [oQ, qMap] = mef_ssim_d(imgSeq, fI, C, p, window, structureThres)


imgSeq = double(imgSeq);
fI = double(fI);

[s1, s2, s3] = size(imgSeq);
wSize = size(window,1);
sWindow = ones(wSize) / wSize^2;
xIdxMax = s1-wSize+1;
yIdxMax = s2-wSize+1;
stMap = zeros(xIdxMax, yIdxMax); %static region quality map
dyMap = zeros(xIdxMax, yIdxMax, s3);%dynamic region quality map
for refIdx=1:s3
    %iterative generate pseudo sequence
    numExd =  2*s3-1;
    imgSeqExd = uint8(zeros(s1, s2, numExd));
    imgSeqExd(:,:,1:s3) = uint8(imgSeq);
    count = 0;
    for i = 1 : s3
        if i ~= refIdx
            count = count + 1;
            temp = imhistmatch(uint8(imgSeqExd(:,:,refIdx)), uint8(imgSeqExd(:,:,i)), 256);
            temp( temp<0 ) = 0;
            temp( temp>255 ) = 255;
            imgSeqExd(:,:,count+s3) = temp;
        end
    end
    
    %compute binary and ideal reference patch
    [out_params] = generate_intermediate(imgSeqExd, C, p, sWindow, structureThres, refIdx);
    lMu = out_params.lMu;
    ed = out_params.ed;
    %     patchIndex = out_params.patchIndex;
    sMap = out_params.sMap;
    maxEd = out_params.maxEd;
    indexMap = out_params.indexMap;
    indM = out_params.indM;
    
    % main loop 
    stepSize = 1;
    xIdx = 1 : stepSize : xIdxMax;
    xIdx = [xIdx xIdx(end)+1 : xIdxMax];
    yIdx = 1 : stepSize : yIdxMax;
    yIdx = [yIdx yIdx(end)+1 : yIdxMax];
    
    offset = wSize-1;
    for row = 1 : length(xIdx)
        for col = 1 : length(yIdx)
            i = xIdx(row);
            j = yIdx(col);
            if indexMap(i, j)==1 && refIdx==1
                blocks = imgSeqExd(i:i+offset, j:j+offset, indM(i,j,:));
                blocks = double(blocks);
                rBlock = zeros(wSize, wSize);
                for k = 1 : s3
                    rBlock = rBlock  + sMap(i, j, k) * (blocks(:,:,k) - lMu(i, j, k) ) / ed(i, j, k);
                end
                if norm(rBlock(:)) > 0
                    rBlock = rBlock / norm(rBlock(:)) * maxEd(i, j);
                end
                fBlock = fI(i:i+offset, j:j+offset);
                rVec = rBlock(:);
                fVec = fBlock(:);
                mu1 = sum( window(:) .* rVec );
                mu2 = sum( window(:) .* fVec );
                sigma1Sq = sum( window(:) .* (rVec - mu1).^2 );
                sigma2Sq = sum( window(:) .* (fVec - mu2).^2 );
                sigma12 = sum(  window(:) .* (rVec - mu1) .* (fVec - mu2)  );
                stMap(i,j) = ( 2 * sigma12 + C ) ./ ( sigma1Sq + sigma2Sq + C );
            elseif indexMap(i, j)==0
                blocks = imgSeqExd(i:i+offset, j:j+offset, indM(i,j,:));
                blocks = double(blocks);
                rBlock = zeros(wSize, wSize);
                for k = 1 : s3
                    rBlock = rBlock  + sMap(i, j, k) * (blocks(:,:,k) - lMu(i, j, k) ) / ed(i, j, k);
                end
                if norm(rBlock(:)) > 0
                    rBlock = rBlock / norm(rBlock(:)) * maxEd(i, j);
                end
                fBlock = fI(i:i+offset, j:j+offset);
                rVec = rBlock(:);
                fVec = fBlock(:);
                mu1 = sum( window(:) .* rVec );
                mu2 = sum( window(:) .* fVec );
                sigma1Sq = sum( window(:) .* (rVec - mu1).^2 );
                sigma2Sq = sum( window(:) .* (fVec - mu2).^2 );
                sigma12 = sum(  window(:) .* (rVec - mu1) .* (fVec - mu2)  );
                dyMap(i,j,refIdx) = ( 2 * sigma12 + C ) ./ ( sigma1Sq + sigma2Sq + C );
            end
        end
    end
end
staticNum = sum(indexMap(:)==1);
dynamicNum = size(indexMap, 1) * size(indexMap, 2) - staticNum;
Qs = sum(sum(stMap/staticNum));
Qs
if dynamicNum~=0
    QdExd = zeros(s3, 1);
    for i = 1 :s3
        QdExd(i) = sum(sum(dyMap(:,:,i)/dynamicNum));
    end    
    [Qd,index]=max(QdExd);
    Qd
    qMap=dyMap(:,:,index) + stMap;    
    oQ = (Qs + Qd) / 2;
else
    oQ = Qs;
end


