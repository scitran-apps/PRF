% Read the mask images (the apertures, they saw words inside the apertures)
% This file is the same for everybody
BAR = load('~/soft/PRF/local/sub-23MAGNO6134T2/stimuli/maskimages.mat');
images = BAR.maskimages;

% Load bcbl stimuli
stimFileBcbl{1} = load('~/soft/PRF/local/sub-23MAGNO6134T2/stimuli/subj6134_T2_run1_exp103.mat');
stimFileBcbl{2} = load('~/soft/PRF/local/sub-23MAGNO6134T2/stimuli/subj6134_T2_run2_exp103.mat');
stimFileBcbl{3} = load('~/soft/PRF/local/sub-23MAGNO6134T2/stimuli/subj6134_T2_run3_exp103.mat');

for ns=1:length(stimFileBcbl)
    % And how many images from the beginning will be discarded
    params{ns}.removeImagesN   = 4;
    % Maddi se aseguro de que termine bien, por lo tanto, esto dura: 161 * 1.82 segundos
    params{ns}.scanDuration = 300;
    params{ns}.totalScans   = 165;
    params{ns}.tr           = 1.82;
    config{ns}.config.rescaleStimuliImages = true;
    config{ns}.config.imageSideSize = 100;  
    
    
    % Every 1.82 an image was acquired, and every 1.82 seconds a different image
    % was shown. Select the image that was shown every 1.82. We know that the
    % response is later, though
    
    % Create the 2D+time dataset of stimulus images over TRs.
    % We need to sub-sample to have the same number of images as the fMRI timepoints
    % Create the image index that will select the images shown when the acquisition
    % took place. Do step by step for clarity:
%     imIndx{ns} = stimFileBcbl{ns}.frameorder(2,:);
%     imIndx{ns} = imIndx{ns}(imagesPerSecond{ns}:imagesPerSecond{ns}:4500);
%     imIndx{ns} = imIndx{ns}(params{ns}.tr:params{ns}.tr:params{ns}.scanDuration+1);
%     imIndx{ns} = imIndx{ns}((params{ns}.removeImagesN + 1):params{ns}.totalScans);
%     imIndx{ns}(imIndx{ns}==0)=1;

    % This is what has been shown
    allImgRefs = stimFileBcbl{ns}.frameorder(2,:);
    % Every tr=1.82 an image has been shown, what is the index of this image?
    % We need to obtain 161 indices, and selected the corresponding number from
    % the existing 4500 images. 
    tmpInd = round(1:length(allImgRefs)/params{ns}.totalScans:length(allImgRefs));
    % Remove the images that correspond to the acquisitions we removed
    ImgRefs = allImgRefs(tmpInd((params{ns}.removeImagesN + 1):end));
    % ImgRefs(ImgRefs==0)=1;
   
    % Create the new trAdjImages stimuli set
    trAdjImages{ns} = zeros([size(images,1),size(images,2),size(ImgRefs,2)]);
    trAdjImages{ns}(:,:,ImgRefs~=0) = images(:,:,ImgRefs(ImgRefs~=0));
    
    temp = zeros(config{ns}.config.imageSideSize, config{ns}.config.imageSideSize, size(stimFileBcbl,3));
    for p=1:size(trAdjImages{ns}, 3)
        temp(:,:,p) = imresize(trAdjImages{ns}(:,:,p),[config{ns}.config.imageSideSize config{ns}.config.imageSideSize],'cubic');
    end
    trAdjImages{ns} = temp;

    % ensure that all values are between 0 and 1
    trAdjImages{ns}(trAdjImages{ns} < 0) = 0;
    trAdjImages{ns}(trAdjImages{ns} > 1) = 1;
 
end
stimulusbcbl = trAdjImages;

stim = stimulusbcbl{1};



    % Prepare the new file.
    vidObj = VideoWriter('local/stimuli.avi');
    open(vidObj);
    % Create an animation.
    axis equal tight off
    set(gca,'nextplot','replacechildren');
    for k = 1:size(stim,3)
        imagesc(stim(:,:,k)); colormap gray;
       % Write each frame to the file.
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
    end
  
    % Close the file.
    close(vidObj);