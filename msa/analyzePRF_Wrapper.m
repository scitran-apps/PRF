function [] = analyzePRF_Wrapper(configFile)
% This function will take the files from FW and prepare them to be input to the
%       analyzePRF code.
%       Example call:
%           analyzePRF_Wrapper('~/soft/PRF/config.json')
%
% This function is based on the  top-level script sent by John Winawer,
% used for the HCP retinotopy dataset, which Noah Benson ran on the NYU
% cluster: https://osf.io/ps532/
% The main command inside this script is the following line, which takes
% time series and stimuli as input (and TR and a code for the kind of seed)
% and produces an array of pRF parameters as output.
%       a1 = analyzePRF(stimulus,data,tr,struct('seedmode',2));
% A modest improvement might be to allow the input data to be a path to a
% NIFTI file or files rather than a data array, and the output to be NIFTIs
% rather than a Matlab struct. That could also be handled by a simple wrapper.
% This script shows how analyzePRF_HCP7TRET was used to analyze the
% HCP 7T Retinotopy time-series data.
% This script has been modified so that it can read:
%   Stimuli: the output logs from vistadisp in .mat format.
%   Data   : the output files from the fmriPrep pipeline gear in FW
%
% TESTING
%{
% configFile = '~/soft/PRF/config.json';
subBase = '~/soft/PRF/local/sub-14magno7806';
configFile = fullfile(subBase, 'config.json');
stimDir    = fullfile(subBase, 'config.json');
analyzePRF_Wrapper(configFile)








%}
%
% Code dependencies:
% - analyzePRF
% - knkutils (http://github.com/kendrickkay/knkutils/)
% - MatlabCIFTI
% - workbench
%
% 2018: GLU, Vistalab garikoitz@gmail.com


%% Check if there is a configFile, if there is read it
if ~exist('configFile', 'var') || ~exist(configFile, 'file')
    configFile = '/flywheel/v0/config.json';
    if ~exist(configFile, 'file')
        error('There is no config file.')
    end
end
% Read the config file
config = jsonread(configFile);


%% Deal with the stimuli part
% Load "images_******.mat", this loads a variable called "images".
% "images" is a 3D matrix of size 768x768x2560.
% It is a matrix of images: 2560 images that are 768x768
load(config.inputs.Stimuli.images.location.path)

% Load "params_********.mat"
% This loads a bunch of variables, one of which is called "params"
% params.seq is a vector of length 4500 (15 frames per second for 300 seconds).
% The values of params.seq range from 1 to 2560. It tells you what image was
% shown at what time point.
% **Note** that the time series that you are working with has 144 time points.
% Each TR is 2 seconds, so this amounts to 288 seconds. This is less than 300
% seconds because the first 6 frames were clipped from the data (scanner data is
% noisy at the start of each run), so 6*2secs = 12secs

load(config.inputs.Stimuli.params.location.path)

% TODO: there is no seq or no structure which is 4500 long. It is called
% stimulus.seq

% Try to get the parameters from the file.
% images       : all possible images (size(images, 3) = 2580 in this example)
% stimulus.seq : says what image was shown
% params.tr    : time between whole brain image acquisitions

% We can calculate how many image images per second where displayed
imagesPerSecond = size(stimulus.seq, 2) / params.scanDuration;
% And how many images from the beginning will be discarded
removeImagesN   = params.prescanDuration / params.tr;

% Create the 2D+time dataset of stimulus images over TRs.
% We need to sub-sample to have the same number of images as the fMRI timepoints
% Create the image index that will select the images shown when the acquisition
% took place. Do step by step for clarity:
imIndx = stimulus.seq;
imIndx = imIndx(imagesPerSecond:imagesPerSecond:4500);
imIndx = imIndx(params.tr:params.tr:300);
imIndx = imIndx((removeImagesN + 1):150);
% Create the new trAdjImages stimuli set
trAdjImages = images(:,:,imIndx);

% See if we do need to make the stimuli file smaller (for copmutational purposes)
if config.config.rescaleStimuliImages
    %{
    % We do not want to use imresize because it will interpolate the values
    % Instead we will find the closest possible value and reduce the filesize manually
    factorX = floor(size(trAdjImages, 1) / config.config.imageSizeX);
    factorY = floor(size(trAdjImages, 2) / config.config.imageSizeY);
    if (factorX < 1); factorX=1; end; % Just in case the desired is bigger, do nothing
    if (factorY < 1); factorY=1; end;
    trAdjImages = trAdjImages(1:factorX:end, 1:factorY:end, :);
    %}
    % To the exact size
    temp = zeros(config.config.imageSideSize, config.config.imageSideSize, size(stimulus,3));
    for p=1:size(trAdjImages, 3)
        temp(:,:,p) = imresize(trAdjImages(:,:,p),[config.config.imageSideSize config.config.imageSideSize],'cubic');
    end
    trAdjImages = temp;

    % ensure that all values are between 0 and 1
    trAdjImages(trAdjImages < 0) = 0;
    trAdjImages(trAdjImages > 1) = 1;
end

% For testing, check the images are ok
%{
figure; set(gcf,'Position', [100 100 1000 700])
for ii = 1:69
    subplot(7,10,ii);
    imagesc(trAdjImages(:,:,ii));
    title(ii);
    axis image tight off;
    colormap(gray);
end
%}


% reshape stimuli into a "flattened" format: 69 stimuli x 100*100 positions
% (note: I do not get why is this, I understand is a analyzePRF requirement)
stimulus = reshape(trAdjImages, ...
                   config.config.imageSideSize*config.config.imageSideSize, ...
                   size(trAdjImages, 3))';


%% Deal with the data part







%% Launch the analysis

































% It is assumed that the time-series data have been organized and placed
% into /path/to/HCP7TRET/. Also, it is assumed that group-average subjects
% 999997, 999998, and 999999 have been computed and saved.
%
% We use a modified version of analyzePRF. This is provided as a
% frozen static version called "analyzePRF_HCP7TRET", distinguishing it
% from the version on github (http://github.com/kendrickkay/analyzePRF/).
% Relative to the github version, modifications have been made in the files:
%   analyzePRF.m



% define
path2func = '/Users/glerma/Downloads/5b6cc2ce7daf4d0010cbc2b1/fmriprep/sub-ex11352/ses-201511281621/func/';
% runs = {'tfMRI_RETCCW_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries.nii' ...
%         'tfMRI_RETCW_7T_PA_Atlas_MSMAll_hp2000_clean.dtseries.nii' ...
%         'tfMRI_RETEXP_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries.nii' ...
%         'tfMRI_RETCON_7T_PA_Atlas_MSMAll_hp2000_clean.dtseries.nii' ...
%         'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries.nii' ...
%         'tfMRI_RETBAR2_7T_PA_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
runs = {'sub-ex11352_ses-201511281621_task-91fMRIRetknk_bold_space-T1w_preproc.nii.gz'};
    
% subjs = matchfiles('/path/to/HCP7TRET/??????');  % match 6-digit subject IDs
subjs = matchfiles('/Users/glerma/Downloads/5b6cc2ce7daf4d0010cbc2b1/fmriprep/sub-ex?????');  % match 6-digit subject IDs
tr = 1;                % temporal sampling rate in seconds
pxtodeg = 16.0/200;    % conversion from pixels to degrees
wbcmd = 'wb_command';  % path to workbench command

% define which subject to analyze (1 through 184)
wh = 1;

% define which model fit to perform (1 through 3)
typ = 1;  % 1 is all runs, 2 is first half of each run, 3 is second half of each run

% load stimulus apertures
    aperturefiles = {'RETCCWsmall.mat' ...
                     'RETCWsmall.mat' ...
                     'RETEXPsmall.mat' ...
                     'RETCONsmall.mat' ...
                     'RETBARsmall.mat' ...
                     'RETBARsmall.mat'};
a1 = loadmulti(aperturefiles,'stim',4);
stimulus = splitmatrix(a1,4);
clear a1;

% load data
data = {};
for p=1:6
  data{p} = double(getfield(ciftiopen([subjs{wh} '/' runs{p}],wbcmd),'cdata'));
end

% deal with subsetting
switch typ
case 1
case 2
  stimulus = cellfun(@(x) x(:,:,1:150),stimulus,'UniformOutput',0);
  data =     cellfun(@(x) x(:,1:150),  data,    'UniformOutput',0);
case 3
  stimulus = cellfun(@(x) x(:,:,151:300),stimulus,'UniformOutput',0);
  data =     cellfun(@(x) x(:,151:300),  data,    'UniformOutput',0);
end

% fit the models
% a1 = analyzePRF(stimulus,data,tr,struct('seedmode',2));



% This is the stim file that outputs the knk tool
stimulus  = load('~/soft/PRF/local/sub-14magno7806/stimuli/20190303192310_subj7806_run1_exp103.mat');
% imgBOLD   = niftiRead('~/Downloads/5b6cc2ce7daf4d0010cbc2b1/fmriprep/sub-ex11352/ses-201511281621/func/sub-ex11352_ses-201511281621_task-91fMRIRetknk_bold_space-T1w_preproc.nii.gz')
% data{1}   = imgBOLD.data(:,:,:,7:end);
% BAR       = load('~/Downloads/apertures/RETBARsmall.mat')

% Testing spanish data, automate this
% allData = gifti('~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_space-fsnative_hemi-L.func.gii');
% !mri_convert ~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_space-fsnative_hemi-L.func.gii ~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_space-fsnative_hemi-L.func.mgh
% !mri_convert ~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-02_space-fsnative_hemi-L.func.gii ~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-02_space-fsnative_hemi-L.func.mgh
% !mri_convert ~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-03_space-fsnative_hemi-L.func.gii ~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-03_space-fsnative_hemi-L.func.mgh
allData = MRIread('~/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_space-fsnative_hemi-L.func.mgh');




% Remove the initial images that do not apply
data    = {};
data{1} = squeeze(allData.vol(:,:,:,5:end));
% Do we need to reshape to make the code below work? Do it here or delete this line 
% data{1} = data{1};

% Read the stimuli
BAR = load('~/soft/PRF/local/sub-14magno7806/stimuli/maskimages.mat');
images = BAR.maskimages;

% And how many images from the beginning will be discarded
removeImagesN   = 4;



% MAddi se aseguro de que termine bien, por lo tanto, esto dura: 161 * 1.82 segundos
params.scanDuration = 300;
params.totalScans   = 165;
params.tr           = 1.82;
% And how many images from the beginning will be discarded
params.removeImagesN   = 4;
config.config.rescaleStimuliImages = true;
config.config.imageSideSize = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% THIS IS DONE ABOVE, DELETE AFTER TESTS %%%%%%%%
%%%%%%%%%%%% CAREFUL, ABOVE IS USING THE OLD STIM FILES %%%%
%%%%%%%%%%%% BELOW ADAPTED TO KNK %%%%%%%%%%%%%%%%%%%%%%%%%%


% TODO: there is no seq or no structure which is 4500 long. It is called
% stimulus.seq

% Try to get the parameters from the file.
% images       : all possible images (size(images, 3) = 2580 in this example)
% stimulus.seq : says what image was shown
% params.tr    : time between whole brain image acquisitions

% We can calculate how many image images per second where displayed
imagesPerSecond = size(stimulus.frameorder, 2) / params.scanDuration;



% Create the 2D+time dataset of stimulus images over TRs.
% We need to sub-sample to have the same number of images as the fMRI timepoints
% Create the image index that will select the images shown when the acquisition
% took place. Do step by step for clarity:
imIndx = stimulus.frameorder(2,:);
imIndx = imIndx(imagesPerSecond:imagesPerSecond:4500);
imIndx = imIndx(params.tr:params.tr:params.scanDuration+1);
imIndx = imIndx((params.removeImagesN + 1):params.totalScans);
imIndx(imIndx==0)=1;
% Create the new trAdjImages stimuli set
trAdjImages = images(:,:,int16(abs(imIndx)));


% See if we do need to make the stimuli file smaller (for copmutational purposes)
if config.config.rescaleStimuliImages
    %{
    % We do not want to use imresize because it will interpolate the values
    % Instead we will find the closest possible value and reduce the filesize manually
    factorX = floor(size(trAdjImages, 1) / config.config.imageSizeX);
    factorY = floor(size(trAdjImages, 2) / config.config.imageSizeY);
    if (factorX < 1); factorX=1; end; % Just in case the desired is bigger, do nothing
    if (factorY < 1); factorY=1; end;
    trAdjImages = trAdjImages(1:factorX:end, 1:factorY:end, :);
    %}
    % To the exact size
    temp = zeros(config.config.imageSideSize, config.config.imageSideSize, size(stimulus,3));
    for p=1:size(trAdjImages, 3)
        temp(:,:,p) = imresize(trAdjImages(:,:,p),[config.config.imageSideSize config.config.imageSideSize],'cubic');
    end
    trAdjImages = temp;

    % ensure that all values are between 0 and 1
    trAdjImages(trAdjImages < 0) = 0;
    trAdjImages(trAdjImages > 1) = 1;
end

% For testing, check the images are ok
%{
figure; set(gcf,'Position', [100 100 1000 700])
for ii = 1:77
    subplot(7,11,ii);
    imagesc(trAdjImages(:,:,ii));
    title(ii);
    axis image tight off;
    colormap(gray);
end
%}


% reshape stimuli into a "flattened" format: numscans stimuli x 100*100 positions
% (note: I do not get why is this, I understand is a analyzePRF requirement)
% stimulus = reshape(trAdjImages, ...
%                    config.config.imageSideSize*config.config.imageSideSize, ...
%                    size(trAdjImages, 3))';


% Read labels, and only to things within the label
lh_V1_label_fname = '/Users/glerma/soft/PRF/local/sub-14magno7806/freesurfer/sub-14MAGNO7806/label/lh.V1_exvivo.thresh.label';
setenv('SUBJECTS_DIR','/Users/glerma/soft/PRF/local/sub-14magno7806/freesurfer')
lh_V1_label = myFSread_label('sub-14MAGNO7806', lh_V1_label_fname, true);
lh_V1_label_ind = lh_V1_label(:,1);

lhwhite_fname = '/Users/glerma/soft/PRF/local/sub-14magno7806/freesurfer/sub-14MAGNO7806/surf/lh.white';


               
               
% myStim = double(BAR.stim(:,:,1:144));
myData = double(data{1});

myDataV1 = myData(lh_V1_label_ind+1, :);
myStim = trAdjImages;
tr     = 1.82;
a1     = analyzePRF(myStim,myDataV1,tr,struct('seedmode',0));

% Just in case
% save('/Users/glerma/soft/PRF/local/sub-14magno7806/results/a1.mat','a1')
% load('/Users/glerma/soft/PRF/local/sub-14magno7806/results/a1.mat')
% Save it as a FS file

a1toSurf = allData;
a1toSurf.vol     = zeros(size(a1toSurf.vol));
a1toSurf.vol     = a1toSurf.vol(:,:,:,1:6);
a1toSurf.nframes = 6;
% assign the values of interest
pxtodeg = 16.0/200;    % conversion from pixels to degrees
a1toSurf.vol(:,lh_V1_label_ind + 1,1,1) = a1.ang;
tmpecc = a1.ecc*pxtodeg;     % convert to degrees;
% whenever eccentricity is exactly 0, we set angle to NaN since it is ill-defined.
tmpecc(tmpecc==0) = NaN;
a1toSurf.vol(:,lh_V1_label_ind + 1,1,2) = tmpecc;
a1toSurf.vol(:,lh_V1_label_ind + 1,1,3) = a1.expt;
a1toSurf.vol(:,lh_V1_label_ind + 1,1,4) = a1.rfsize*pxtodeg;  % convert to degrees
a1toSurf.vol(:,lh_V1_label_ind + 1,1,5) = a1.R2;
a1toSurf.vol(:,lh_V1_label_ind + 1,1,6) = a1.gain;


MRIwrite(a1toSurf, '/Users/glerma/soft/PRF/local/sub-14magno7806/results/a1.mgh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







































% prepare outputs
quants = {'ang' 'ecc' 'gain' 'meanvol' 'R2' 'rfsize'};
allresults = zeros(91282,length(quants),length(subjs),3,'single');  % 91282 x 6 x 184 x 3
allresults(:,1,wh,typ) = a1.ang;
allresults(:,2,wh,typ) = a1.ecc*pxtodeg;     % convert to degrees
allresults(:,3,wh,typ) = a1.gain;
allresults(:,4,wh,typ) = a1.meanvol;
allresults(:,5,wh,typ) = a1.R2;
allresults(:,6,wh,typ) = a1.rfsize*pxtodeg;  % convert to degrees

% one final modification to the outputs:
% whenever eccentricity is exactly 0, we set angle to NaN since it is ill-defined.
allresults = squish(permute(allresults,[1 3 4 2]),3);  % 91282*184*3 x 6
allresults(allresults(:,2)==0,1) = NaN;
allresults = permute(reshape(allresults,[91282 184 3 6]),[1 4 2 3]);






end
