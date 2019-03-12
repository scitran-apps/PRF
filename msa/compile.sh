#!/bin/bash
# module load matlab/2017a

cat > build.m <<END

addpath(genpath('/data/localhome/glerma/soft/AFQ'));
addpath(genpath('/data/localhome/glerma/soft/afq-pipeline'));
addpath(genpath('/data/localhome/glerma/soft/vistasoft'));
addpath(genpath('/black/localhome/glerma/soft/spm8'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/toolbox/Beamforming'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/toolbox/DARTEL'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/toolbox/MEEGtools'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/toolbox/FieldMap'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/toolbox/dcm_meeg'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/external/bemcp'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/external/ctf'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/external/eeprobe'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/external/fieldtrip'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/external/mne'));
rmpath(genpath('/black/localhome/glerma/soft/spm8/external/yokogawa'));
addpath(genpath('/data/localhome/glerma/soft/jsonlab'));
addpath(genpath('/data/localhome/glerma/soft/encode'));
addpath(genpath('/data/localhome/glerma/soft/JSONio'));
addpath(genpath('/data/localhome/glerma/soft/app-life'));

mcc -m -R -nodisplay -a /data/localhome/glerma/soft/encode/mexfiles -a /data/localhome/glerma/soft/vistasoft/mrDiffusion -d compiled AFQ_StandAlone_QMR.m
exit
END
Matlabr2017a -nodisplay -nosplash -r build && rm build.m