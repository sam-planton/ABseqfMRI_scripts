function Do_ABseq16fMRI_5_SPM_SynLoca_1stLevel(w,iS)
  
cdd = pwd;

%% Parameters 
w.TR              = 1.81;  % Repetition time (s)
w.hpf             = 128;   % High-pass filter (s)
w.mthresh         = 0.3;   % implicit mask threshold
w.contrast_only   = false;
w.sessions        = {'SynLoca'};   % session directories (parent=functional)

w.with_physIO     = true;
w.with_deriv      = false;
w.with_deriv2     = false;

w.with_art        = false;
w.with_hvc        = false; % include "high variance confounds" regressors
w.new_confounds   = 'high_variance_confounds_1.0.csv';

%%
fprintf([' \n \n']);       
fprintf('========================================================================\n');
fprintf(['              ' w.subjects{iS} ' SPM SynLoca 1st level...\n']);       
fprintf('========================================================================\n');

%%  Subject Directories
w.subPath          = fullfile (w.niftidir,  w.subjects{iS}); 
w.T1Path           = fullfile (w.subPath, 'anat');
w.funcPath         = fullfile (w.subPath, 'func');
w.stimPath         = fullfile (w.datadir, 'SPM_onsets', w.subjects{iS}, filesep); 
if ~exist(w.stimPath, 'dir'); mkdir(w.stimPath);end
    
%% ==== SubFunctions to run ==== %  
SynLoca_Generate_SPM_onsets(w, iS) % copy onsets from "localizer.mat" to the subject's directory and prepare contrasts (because they are supposed to be always the same!)
DoPhysIO(w,iS)
DoFirstLevel(w, iS) % one subject full analysis
%%===============================%  
DoReCreateContrastLabelsFile(w, iS)

   
%% Create "confounds" regressors with PhysIO toolbox (including PCA of white matter and csf, multiple rp, and marking volumes with high mvt...)
function DoPhysIO(w, iS)
    
    %==============================================================%
    %  PhysIO
    %==============================================================% 

    % Parameters of EPI files	  
    tmp = load(fullfile(w.niftidir, w.subjects{iS}, 'func', 'SliceTimingInfo.mat'));
    w.nSlices         = tmp.NumberOfSlices;     %
    w.TR              = tmp.TR;  % Repetition time (s)
    w.thickness       = tmp.SliceThickness;    % Slice thickness (mm)
    w.READOUT_TIME    = tmp.total_readout_time_spm*1000 ;
    w.slice_times     = tmp.SliceTiming;
    w.refslice        = w.slice_times(floor(size(w.slice_times,2)/2));  % refSlice is the 34th on 69... /// this is not used anymore

    % ROIs for noise PCA
    wc2 = spm_select('FPList', w.T1Path, ['^wc2'  '.*\.nii$']); 
    wc3 = spm_select('FPList', w.T1Path, ['^wc3'  '.*\.nii$']); 

    % The toolbox will force coregistration, but images become misaligned
    % after coregistration !!! and even the original wc images get modified
    % (maybe just some orientation info) ! WHY ?? Same thing if
    % coregistration is done before

    % Possible solution: recreate the normalized tissue images (wc2, wc3) using
    % the voxel size of the EPI
    c2 = spm_select('FPList', w.T1Path, ['^c2'  '.*\.nii$']); 
    c3 = spm_select('FPList', w.T1Path, ['^c3'  '.*\.nii$']); 
    forwardDeformation = spm_select('FPList', w.T1Path, ['^y_.*\.nii$']); 
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {c2;c3};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1.75 1.75 1.75];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch); 


    clear matlabbatch;
    stage = 1;
    for j=1:numel(w.sessions) % Sessions loop

        % Get func images
        w.funcPathSession = fullfile (w.funcPath, w.sessions{j});
        unsmoothedEPI = spm_select('ExtFPList',  fullfile(w.funcPathSession), ['^wua' '.*\.nii$'], Inf); 

        % Get the head movements file
        rpF = spm_select('FPList',  fullfile(w.funcPathSession), ['^rp_' '.*\.txt$']);  

        % TAPAS PhysIO SPM batch
        matlabbatch{stage}.spm.tools.physio.save_dir = {w.funcPathSession};
        matlabbatch{stage}.spm.tools.physio.log_files.vendor = 'Siemens';
        matlabbatch{stage}.spm.tools.physio.log_files.cardiac = {''};
        matlabbatch{stage}.spm.tools.physio.log_files.respiration = {''};
        matlabbatch{stage}.spm.tools.physio.log_files.scan_timing = {''};
        matlabbatch{stage}.spm.tools.physio.log_files.sampling_interval = [];
        matlabbatch{stage}.spm.tools.physio.log_files.relative_start_acquisition = 0;
        matlabbatch{stage}.spm.tools.physio.log_files.align_scan = 'last';

        matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nslices = numel(w.slice_times);
        matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
        matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.TR = w.TR;
        matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
        matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nscans = size(unsmoothedEPI,1);
        matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.onset_slice = 1;

        matlabbatch{stage}.spm.tools.physio.model.output_multiple_regressors = [w.subjects{iS} '_run_' num2str(j, '%02.f') '_multiple_regressors.txt'];
        matlabbatch{stage}.spm.tools.physio.model.output_physio = [w.subjects{iS} '_run_' num2str(j, '%02.f') '_physio.mat'];
        matlabbatch{stage}.spm.tools.physio.model.orthogonalise = 'none';
        matlabbatch{stage}.spm.tools.physio.model.censor_unreliable_recording_intervals = false; %new
        matlabbatch{stage}.spm.tools.physio.model.retroicor.no = struct([]); %new
        matlabbatch{stage}.spm.tools.physio.model.rvt.no = struct([]); %new
        matlabbatch{stage}.spm.tools.physio.model.hrv.no = struct([]); %new
        
        matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.fmri_files = cellstr(unsmoothedEPI);
        matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.roi_files = {wc2; wc3}; % white matter / CSF
        matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.force_coregister = 'No'; %'Yes';
        matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.thresholds = [0.99 0.99]; % Threshold for the ROI (was using 0.90 0.90 before)
        matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_voxel_crop = [1 0]; % N voxels to crop
        matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_components = 5; %

        matlabbatch{stage}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {rpF};
        matlabbatch{stage}.spm.tools.physio.model.movement.yes.order = 6;
        matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_method = 'FD';
        matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_threshold = 0.5;  %

        matlabbatch{stage}.spm.tools.physio.model.other.no = struct([]);
        matlabbatch{stage}.spm.tools.physio.verbose.level = 3; % 3 to get all figures, but Matlab may crash if too many..
        matlabbatch{stage}.spm.tools.physio.verbose.fig_output_file = [w.subjects{iS} '_run_' num2str(j, '%02.f') '_PhysIO_figure.fig'];
        matlabbatch{stage}.spm.tools.physio.verbose.use_tabs = true;

        stage = stage +1;
    end

    w.batchdir = fullfile (w.datadir, '2_PREPROC', w.subjects{iS}); 
    save(fullfile(w.batchdir, 'SPM12_matlabbatch_PhysIO.mat'),'matlabbatch');  
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  

end

%% First Level
function DoFirstLevel(w, iS)
    
    %% Output folder
    if     w.with_hvc && ~w.with_art
        model_folder = ['With_' w.new_confounds(1:end-4)];
    elseif w.with_art && ~w.with_hvc
        model_folder = 'With_Art';
    elseif w.with_art && w.with_hvc
        model_folder = ['With_Art&hvc_' w.new_confounds(1:end-4)];
    elseif w.with_physIO
        model_folder = 'PhysIO';
    else
        model_folder = 'Standard';
    end
    w.firstDir = fullfile (w.datadir, '3_FIRST_LEVEL',  'SynLoca', w.subjects{iS}, model_folder);
    if w.with_deriv
        w.firstDir = fullfile (w.firstDir, 'with_deriv');
    elseif w.with_deriv2
        w.firstDir = fullfile (w.firstDir, 'with_deriv2');
    end

    if isfolder (w.firstDir) && w.contrast_only == false; delete ([w.firstDir '/*']); else mkdir(w.firstDir); end  
    disp(w.firstDir)   

    %% 1st level processing
    stage = 1;
    cd(w.firstDir); 
    
    %% Get explicit mask
    explicitMask = spm_select('FPList', w.T1Path, '^brexplicitMask_wc1wc2wc3');
    if isempty(explicitMask); error('Can''t load  explicit mask'); end

    %%  SPM batch
    if w.contrast_only == false
        
    %==============================================================%
    %  fMRI model specification
    %==============================================================% 

    clear matlabbatch;
    
    matlabbatch{stage}.spm.stats.fmri_spec.dir =  {w.firstDir};
    matlabbatch{stage}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{stage}.spm.stats.fmri_spec.timing.RT = w.TR;
    matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t0 = 1;   % First time bin because reference for slice timing was T=0ms    
    
    w.cond_list = []; cond_list = [];
    for j=1:numel(w.sessions) % Sessions loop // not really needed here..., one session
 
        w.funcPathSession         = fullfile (w.funcPath, w.sessions{j});
       
        if w.with_physIO == false

            % Get the head movements file (store in a matrix)
            rpF = spm_select('FPList',  fullfile(w.funcPathSession), ['^rp_' '.*\.txt$']);
            R = load(deblank(rpF));

            % Get the high_variance_confounds generated with Nilearn in Python (and add them to the confound struct)
            % Or the ones with Art
            if     w.with_hvc && ~w.with_art
                R = [R load(fullfile(w.funcPathSession, w.new_confounds))];
            elseif w.with_art && ~w.with_hvc
                f = dir(fullfile(w.funcPathSession, 'art_regression_outliers_and_movement_*'));
                load(fullfile(f.folder,f.name));
            elseif w.with_art && w.with_hvc
                f = dir(fullfile(w.funcPathSession, 'art_regression_outliers_and_movement_*'));
                load(fullfile(f.folder,f.name));
                R = [R load(fullfile(w.funcPathSession, w.new_confounds))];
            end

            %%Save confounds in a mat file
            save(fullfile(w.funcPathSession,'confounds.mat'), 'R');
            confounds_file = cellstr(spm_select('FPList', w.funcPathSession, ['confounds.*\.mat$']));

        elseif w.with_physIO == true

            % Get the movement/noise regressors file (store in a matrix)
            rpF = fullfile(w.funcPathSession, [w.subjects{iS} '_run_' num2str(j, '%02.f') '_multiple_regressors.txt']);
            R = load(deblank(rpF));

            %%Save confounds in a mat file
            save(fullfile(w.funcPathSession,'confounds.mat'), 'R');
            confounds_file = cellstr(spm_select('FPList', w.funcPathSession, ['confounds.*\.mat$']));

        end

        % Get stimulation onsets (from session number) and store conditions list
        onset_file = {};
        onset_file = cellstr(spm_select('FPList', w.stimPath, ['syntaxlocalizer_multicond.*' '.*\.mat$']));
        if strcmp(onset_file, ''); error('cannot find onset_file'); end
        w.onset_file = load(onset_file{1}); 
        if w.with_deriv
            A = w.onset_file.names';
            B = repmat({'rp'},size(w.onset_file.names',1),1);
            w.cond_list = [w.cond_list; reshape([A';B'],1,[])'; repmat({'rp'},size(R,2),1)];
        elseif w.with_deriv2
            A = w.onset_file.names';
            B = [];
            for kk=1:numel(A)
                B = [B; A(kk); {'rp'}; {'rp'}];
            end
            w.cond_list = [w.cond_list; B; repmat({'rp'},size(R,2),1)];
        else               
            w.cond_list = [w.cond_list; w.onset_file.names'; repmat({'rp'},size(R,2),1)];
        end

        %% Get EPI smoothed images         
        EPI = spm_select('ExtFPList',  fullfile(w.funcPathSession), ['^swua' '.*\.nii$'], Inf); 

        %% Fill SPM batch (sessions)
        matlabbatch{stage}.spm.stats.fmri_spec.sess(j).scans = cellstr(EPI);  
        matlabbatch{stage}.spm.stats.fmri_spec.sess(j).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        matlabbatch{stage}.spm.stats.fmri_spec.sess(j).multi = onset_file;
        matlabbatch{stage}.spm.stats.fmri_spec.sess(j).regress = struct('name', {}, 'val', {});        
        matlabbatch{stage}.spm.stats.fmri_spec.sess(j).multi_reg = confounds_file; 
        matlabbatch{stage}.spm.stats.fmri_spec.sess(j).hpf = w.hpf;      %%%%%%%%%%%%%%%%%%%

    end
    
    matlabbatch{stage}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{stage}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; %%%%%%%%%%%%%%%%%%%%% 
    if w.with_deriv
        matlabbatch{stage}.spm.stats.fmri_spec.bases.hrf.derivs  = [1 0];
    end
    if w.with_deriv2
        matlabbatch{stage}.spm.stats.fmri_spec.bases.hrf.derivs  = [1 1];
    end
    matlabbatch{stage}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{stage}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{stage}.spm.stats.fmri_spec.mthresh = w.mthresh;
    matlabbatch{stage}.spm.stats.fmri_spec.mask = {explicitMask};
    matlabbatch{stage}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
    cond_list = w.cond_list;
    save('cond_list.mat', 'cond_list')
        
    stage = stage +1;
    
    %==============================================================%
    %  Model Estimation
    %==============================================================%   
    
    matlabbatch{stage}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{stage}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{stage}.spm.stats.fmri_est.method.Classical = 1;
    
    stage = stage +1;
    
    end
    
    %==============================================================%
    %  Contrast manager
    %==============================================================%   

    %========== function to 'normalize' contrast
%     normcon = @(x) transpose(x/numel(find(x>0)));
    normcon = @(x) (x/numel(find(x>0)));

    %========== Retrieve conditions & prepare contrast weights
    load('cond_list.mat')
    nconds = numel(cond_list);
    conds = eye(nconds);
    rp_conds = double(strcmp(cond_list,'rp')');   

    %====== Create 0/1 vectors for all possible conditions
    %-- Sequence
    Audio_betas       = contains(cond_list,'Audio');
    Video_betas       = contains(cond_list,'Video');
    Motor_betas       = contains(cond_list,'Click');
    Computation_betas = contains(cond_list,'Computation');
    Sentences_betas   = contains(cond_list,'Sentences');
    Motor_Left        = contains(cond_list,'Left');
    Motor_Right       = contains(cond_list,'Right');
    SentencesControl  = contains(cond_list,'Rotated') | contains(cond_list,'False_Font');

    %==========  Prepare contrasts
    labels = {'Motor',...
              'Computation', ...
              'Sentences',...
              'Audio Sentences - Rotated Audio',...
              'Video Sentences - False Font Video',...        
              'Left Click - Right Click',...
              'Right Click - Left Click',...
              'Video - Audio',...
              'Audio - Video',...
              'Motor - Cognitive',...
              'Cognitive - Motor',...
              'Audio Computation - Audio Sentences',...
              'Audio Sentences - Audio Computation',...
              'Video Computation - Video Sentences',...
              'Video Sentences - Video Computation',...
              'Computation - Sentences',...
              'Sentences - Computation',...
              'Sentences - controls',...
              'controls - Sentences';
              };
    convalues = [normcon(Motor_betas)';
                 normcon(Computation_betas)';
                 normcon(Sentences_betas)';
                 normcon((Audio_betas & Sentences_betas) - (Audio_betas & SentencesControl))';
                 normcon((Video_betas & Sentences_betas) - (Video_betas & SentencesControl))';
                 normcon(Motor_Left - Motor_Right)';
                 normcon(Motor_Right - Motor_Left)';
                 normcon(Video_betas - Audio_betas)';
                 normcon(Audio_betas - Video_betas)';
                 normcon(Motor_betas - (Sentences_betas | Computation_betas))';
                 normcon((Sentences_betas | Computation_betas) - Motor_betas)';
                 normcon((Audio_betas & Computation_betas) - (Audio_betas & Sentences_betas))';
                 normcon((Audio_betas & Sentences_betas) - (Audio_betas & Computation_betas))';
                 normcon((Video_betas & Computation_betas) - (Video_betas & Sentences_betas))';
                 normcon((Video_betas & Sentences_betas) - (Video_betas & Computation_betas))';
                 normcon(Computation_betas - Sentences_betas)';
                 normcon(Sentences_betas - Computation_betas)';
                 normcon(Sentences_betas - SentencesControl)';
                 normcon(SentencesControl - Sentences_betas)';
                 ];

    % add beta contrasts for each unique condition
    labels_add = {}; convalues_add = [];
    kk = 1;
    for jj=1:numel(cond_list)
        if ~strcmp(cond_list{jj}, 'rp')
            labels_add{kk} = cond_list{jj};
            convalues_add(kk, :) = normcon(contains(cond_list,cond_list{jj}))';
            kk = kk + 1;
        end
    end
    labels = [labels labels_add];
    convalues = [convalues; convalues_add];
    
    %==========  SPM batch
    matlabbatch{stage}.spm.stats.con.spmmat = cellstr(fullfile(w.firstDir,'SPM.mat')); 
    
    %==========  Fill contrasts
    n_cont = 1;        
    %-------% Effects of interests F contrast
    task_conds   = eye(nconds); % All conditions
    task_conds(find(rp_conds==1),:)=[]; % Removing realignment parameters columns
    matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_CONDITIONS';
    matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = task_conds;
    matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;
    %-------% All T contrasts
    for ii = 1:numel(labels)
        %-------%
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = labels{ii};
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = convalues(ii,:);
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none'; 
        n_cont = n_cont+1;
        %-------%
    end
    disp(['Number of contrasts = ' num2str(n_cont-1)])
    

    %-------%
    matlabbatch{stage}.spm.stats.con.delete = 1;
    if w.contrast_only == false
        save(fullfile(w.firstDir, 'SPM12_1stLevel_matlabbatch.mat'),'matlabbatch');  
    end
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  

end

%% Subfunction to copy onsets from "localizer.mat" to the subject's directory
function SynLoca_Generate_SPM_onsets(w, iS)

    % Original onset file
    if ispc
        % Original onset file
        matfile = 'X:\Scripts\Experiment\Localizer\localizer.mat';
    elseif isunix
        % Original onset file
        matfile = '/neurospin/unicog/protocols/IRMf/ABseq16fMRI_PlantonAlRoumi_2020/Scripts/Experiment/Localizer/localizer.mat';
    end
    onset_info = load(matfile);

    newnames = {'Rotated_Audio',... 
                'False_Font_Video',... 
                'Right_Audio_Click',...
                'Left_Audio_Click',...
                'Right_Video_Click',...
                'Left_Video_Click',...
                'Audio_Computation',...
                'Video_Computation',...
                'Video_Sentences',...
                'Audio Sentences'};

    % Save new file
    clear names onsets durations
    names            = newnames; %onset_info.names;
    onsets           = onset_info.onsets;
    durations        = onset_info.durations;
    fname= fullfile(w.stimPath, [w.subjects{iS} '_syntaxlocalizer_multicond.mat']);
    save(fname, 'names', 'onsets', 'durations')

end

function DoReCreateContrastLabelsFile(w, iS)
    %% Output folder
    if w.with_physIO
        model_folder = 'PhysIO';
    else
        model_folder = 'Standard';
    end
    w.firstDir = fullfile (w.datadir, '3_FIRST_LEVEL',  'SynLoca', w.subjects{iS}, model_folder);
    if w.with_deriv
        w.firstDir = fullfile (w.firstDir, 'with_deriv');
    elseif w.with_deriv2
        w.firstDir = fullfile (w.firstDir, 'with_deriv2');
    end  
    tmp = load(fullfile(w.firstDir, 'SPM.mat'));
    labels = {tmp.SPM.xCon.name}.';
    save(fullfile(w.firstDir, 'Contrast_labels.mat'),'labels');
    disp([fullfile(w.firstDir, 'Contrast_labels.mat'), ' saved.'])
end

cd(cdd)

end