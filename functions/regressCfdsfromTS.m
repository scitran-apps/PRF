function newts=regressCfdsfromTS(datamat, confoundsmat, varargin)
% Regresses Friston 24 confounds out from fmriprep output. 
% Based on: https://www.mail-archive.com/hcp-users@humanconnectome.org/msg07291.html
% It accepts matrices directly or nifti files, and it can write back the
% corrected file if writeNifti is set to 1
%{
datamat: can be a matrix or a path to a nifti/mgh file, if it is nifti file it will
rewrite with same name _REGRESSED
confoundsmat: can be a matrix or a path name to the .tsv output file
Optional parameters:

datamat      = '/Users/glerma/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_space-fsnative_hemi-L.func.mgh';
datamat      = '/Users/glerma/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_space-T1w_desc-preproc_bold.nii.gz';
confoundsmat = '/Users/glerma/soft/PRF/local/sub-14magno7806/fmriprep/sub-14MAGNO7806/ses-20190303/func/sub-14MAGNO7806_ses-20190303_task-ret_run-01_desc-confounds_regressors.tsv';
writeNifti   = true;

    newts = regressCfdsfromTS(datamat, confoundsmat, 'writeNifti', true);


%}
    %% Parse inputs
    p = inputParser;

    addRequired(p, 'datamat');
    addRequired(p, 'confoundsmat');
    addOptional(p, 'writeNifti',false, @islogical);

    parse(p,datamat,confoundsmat,varargin{:});

    writeNifti = p.Results.writeNifti;
    
    %% Do the thing

    % Obtain the matrices
    if ischar(datamat)
        mriFile = MRIread(datamat);
        volFile = mriFile.vol;
        mriSize = size(volFile);
        if mriSize(1) == 1
            ts = squeeze(volFile);
        else      
            ts = reshape(volFile, [mriSize(1)*mriSize(2)*mriSize(3),mriSize(4)]);
        end
    end
    
    % If we pass the .tsv from fmriprep, obtain the friston 24, and write them
    % just in case
    % Confoundsmat will be size(ts,2) in size, and with 24 regressor columns
    if ischar(confoundsmat)
        [FILEPATH,NAME] = fileparts(confoundsmat);
        outputPathName = [FILEPATH filesep NAME '_friston24.txt'];
        confoundsmat = createNewRegressors(confoundsmat,outputPathName);
    end
    
    
    % Check it just in case, if it doesn't work throw an error and ask to check
    if ~isequal(size(ts,2), size(confoundsmat,1))
        error('Data and regressors do not have same size, revise')
    end
    
    
    % Do the regression
    demeanmat = confoundsmat - mean(confoundsmat);
    newts = ts - (demeanmat * (pinv(demeanmat) * ts'))';
    
    % Write the file back if asked to
    if writeNifti && ischar(datamat)
        mriFile.vol = reshape(newts, mriSize);
        [FILEPATH,NAME,EXT] = fileparts(datamat);
        if strcmp(EXT,'.gz')
            NAME = NAME(1:end-4);
            EXT  = '.nii.gz';
        end
        writeNiftiName = [FILEPATH filesep NAME '_REGRESSED' EXT];
        MRIwrite(mriFile, writeNiftiName);
    end
end
