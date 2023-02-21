clear all; close all; clc

% Folders (& spm path)
w.datadir        = '/neurospin/unicog/protocols/IRMf/ABseq16fMRI_PlantonAlRoumi_2020/Data/'; 
w.scriptsdir     = '/neurospin/unicog/protocols/IRMf/ABseq16fMRI_PlantonAlRoumi_2020/Scripts/DataAnalysis/';
spm_path         = '/home/sp253886/MatlabTools/spm12'; addpath(spm_path)
addpath(w.scriptsdir);
addpath(fullfile(w.scriptsdir, 'ProcessingFunctions'));

w.niftidir       = fullfile(w.datadir, '2_PREPROC', filesep);

% Subject
w.subjects  = {                                  'sub_04'  'sub_05'   'sub_06'   'sub_07'   'sub_08'   'sub_09'   'sub_10'  ...
                'sub_11'   'sub_12'   'sub_13'   'sub_14'   'sub_15'   'sub_16'   'sub_17'   'sub_18'   'sub_19'   'sub_20'  ...
                'sub_21'   'sub_22'   'sub_23'   'sub_24'   'sub_25'   'sub_26'};

% NOTE: The required MRI sequences numbers are defined in the subjects loop below

spm('defaults','fmri');  
spm_jobman('initcfg');

%% =================================== %%
%%% LOOP OVER SUBJECTS: RUN FUNCTIONS %%%
%%% ================================= %%%

% Subjects to process:
subs_to_do = 1:numel(w.subjects);

% Subjects loop
for iS=subs_to_do
    
    %=======================================%
    %+++++++++++++++++++++++++++++++++++++++%
    Do_ABseq16fMRI_1_ProcessTopUp_Fieldmap(w,iS)
    Do_ABseq16fMRI_2_SPM_Preprocessing(w,iS)
    %+++++++++++++++++++++++++++++++++++++++%
    %=======================================%
end

load chirp; sound(y,Fs*2/3)
msgbox('%%% ================= Done !!! ================= %%%')
%%% ================================= %%%