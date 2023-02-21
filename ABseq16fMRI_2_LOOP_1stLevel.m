clear all; close all; clc
addpath(pwd);

% warning off
% Folders (& spm path)
w.datadir        = '/neurospin/unicog/protocols/IRMf/ABseq16fMRI_PlantonAlRoumi_2020/Data/'; 
w.scriptsdir     = '/neurospin/unicog/protocols/IRMf/ABseq16fMRI_PlantonAlRoumi_2020/Scripts/DataAnalysis/';
spm_path         = '/home/sp253886/MatlabTools/spm12'; addpath(spm_path)
spm_get_defaults('stats.maxmem', 2^34) % allow to use 16GB RAM
spm_get_defaults('stats.resmem', true); % keep temporary files in memory instead of writing on disk
addpath(w.scriptsdir);
addpath(fullfile(w.scriptsdir, 'ProcessingFunctions'));
addpath(fullfile(w.scriptsdir, 'Utils'));

w.niftidir           = fullfile(w.datadir, '2_PREPROC', filesep);  

w.subjects  = {                           'sub_04' 'sub_05' 'sub_06' 'sub_07' 'sub_08' 'sub_09' 'sub_10'...
               'sub_11' 'sub_12' 'sub_13' 'sub_14' 'sub_15' 'sub_16' 'sub_17' 'sub_18' 'sub_19' 'sub_20'...
               'sub_21' 'sub_22' 'sub_23' 'sub_24' 'sub_25' 'sub_26'};        
spm('defaults','fmri');  
spm_jobman('initcfg');

%% =================================== %%
%%% LOOP OVER SUBJECTS: RUN FUNCTIONS %%%
%%% ================================= %%%

% Subjects to process:
subs_to_do = 4:numel(w.subjects);

%--------------------% Initialize WaitBar %--------------------%
 tic; s = clock; h = waitbar(0,{strrep(mfilename, '_', '-')...
     ' ' ['Loop on subjects - Subject ' strrep(w.subjects{subs_to_do(1)},'_','-') ' in progress']...
     ' ' 'will be updated after each subject...'});   
  movegui(h,'northwest');
%---------------------------------------------------------------%


% Subjects loop
for iS=subs_to_do
    
    switch w.subjects{iS}
        case{'sub_01'}
            syntaxloca = false;  
        otherwise
            syntaxloca = true;  
    end
    
    %=======================================%
    %+++++++++++++++++++++++++++++++++++++++%
    if syntaxloca; Do_ABseq16fMRI_3_SPM_SynLoca_1stLevel(w,iS); end
    Do_ABseq16fMRI_4_SPM_ABseq_1stLevel(w,iS);
    %+++++++++++++++++++++++++++++++++++++++%
    %=======================================%
    
    
    %--------------% Update WaitBar & remaining time %--------------%
    if iS<numel(w.subjects)
        remainTime = (etime(clock,s)/iS) * (numel(w.subjects)-iS);
        endTime = datetime('now','Format','HH:mm')+seconds(remainTime);               
        waitbar(iS/numel(w.subjects),h,{strrep(mfilename, '_', '-') ' '...
            ['Loop on subjects - Subject ' strrep(w.subjects{iS+1},'_','-') ' in progress']... 
            ['Remaining time < ',num2str( remainTime /60,'%.1f' ),' minutes' ]...
            ['End time ~ ' char(endTime)] });
    end
    %---------------------------------------------------------------%

end
delete(h)


% =================================== %%
ABseq16fMRI_3_2ndLevel
% =================================== %%


load chirp; sound(y,Fs*2/3)
msgbox('%%% ================= Done !!! ================= %%%')
%%% ================================= %%%