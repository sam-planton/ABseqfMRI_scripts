function Do_ABseq16fMRI_6_SPM_ABseq_1stLevel(w,iS)

cdd = pwd;

%% Parameters
w.TR              = 1.81;  % Repetition time (s)
w.hpf             = 100;  % High-pass filter (s)
w.mthresh         = 0.3;   % implicit mask threshold
w.sessions        = {'RUN1' 'RUN2' 'RUN3' 'RUN4' 'RUN5'};   % session directories (parent=functional)

w.contrast_only   = false;
w.with_physIO     = true;
w.with_deriv      = false;
w.with_deriv2     = false;
w.concatenate_runs= false;

%%
fprintf(' \n \n');
fprintf('========================================================================\n');
fprintf(['              ' w.subjects{iS} ' SPM ABseq 1st level...\n']);
fprintf('========================================================================\n');


%%  Subject Directories
w.subPath          = fullfile (w.niftidir,  w.subjects{iS});
w.T1Path           = fullfile (w.subPath, 'anat');
w.funcPath         = fullfile (w.subPath, 'func');
w.logPath          = fullfile (w.datadir, 'ExperimentLogs', w.subjects{iS}, filesep);
w.stimPath         = fullfile (w.datadir, 'SPM_onsets', w.subjects{iS}, filesep);
if ~exist(w.logPath , 'dir'); mkdir(w.logPath ); end
if ~exist(w.stimPath, 'dir'); mkdir(w.stimPath); end

%% ==== SubFunctions to run ==== %
ABseq16fMRI_Generate_SPM_onsets(w, iS) % Generate the onsets files using log (.mat) files from the experiment
DoPhysIO(w,iS)
DoFirstLevel_v1(w, iS) % with onsets VERSION 1: sequence*(hab+stand+1dev) - with "overlapping" regressors for dev  - NO MANUAL RESPONSE REG
DoFirstLevel_v2(w, iS) % with onsets VERSION 10: sequence*(hab+stand+1dev) - with "overlapping" regressors for DetectedDev only
%%===============================%

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

        %         % First, coregister the ROI files with the first volume
        %         firstvol = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{1}), ['^wua' '.*\.nii$'], 1);
        %         clear matlabbatch;
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {firstvol};
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.source = {wc2};
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.other = {wc3};
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        %         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        %         spm_jobman('initcfg');
        %         spm_jobman('run',matlabbatch);
        %         % ROIs for noise PCA
        %         wc2 = spm_select('FPList', w.T1Path, ['^rwc2'  '.*\.nii$']);
        %         wc3 = spm_select('FPList', w.T1Path, ['^rwc3'  '.*\.nii$']);

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
            matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_components = 5; % Number of principal components (!!!!! was using 3 in V1)

            matlabbatch{stage}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {rpF};
            matlabbatch{stage}.spm.tools.physio.model.movement.yes.order = 6; % !!!!! was using 12 in v1
            matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_method = 'FD';
            matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_threshold = 0.5;  % Censoring outlier threshold in mm (advice is 0.2mm) !!!!! was using 0.5 in v1 / 0.3 in v2

            matlabbatch{stage}.spm.tools.physio.model.other.no = struct([]);
            matlabbatch{stage}.spm.tools.physio.verbose.level = 3; % 3 to get all figures, but Matlab may crash if too many..
            matlabbatch{stage}.spm.tools.physio.verbose.fig_output_file = [w.subjects{iS} '_run_' num2str(j, '%02.f') '_PhysIO_figure.fig'];
            matlabbatch{stage}.spm.tools.physio.verbose.use_tabs = false; %% false

            stage = stage +1;
        end

        w.batchdir = fullfile (w.datadir, '2_PREPROC', w.subjects{iS});
        save(fullfile(w.batchdir, 'SPM12_matlabbatch_PhysIO.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

%% Model 1: one subject full analysis, with onsets VERSION 1: sequence*(hab+stand+1dev) - NO MANUAL RESPONSE REG
    function DoFirstLevel_v1(w, iS)

        onsets_type = 'onsets_v1';
        which_subfolder = 'Model_1';

        fprintf([' \n \n']);
        fprintf('========================================================================\n');
        fprintf(['  ' w.subjects{iS} ': running Model 1: one subject full analysis using ' onsets_type '...\n']);
        fprintf('========================================================================\n');


        %% Output folder
        if w.with_physIO
            model_folder = fullfile(which_subfolder, 'PhysIO');
        else
            model_folder = fullfile(which_subfolder, 'Standard');
        end
        w.firstDir = fullfile (w.datadir, '3_FIRST_LEVEL',  'ABseq', w.subjects{iS}, model_folder);
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

        % Get explicit mask
        explicitMask = spm_select('FPList', w.T1Path, '^brexplicitMask_wc1wc2wc3');
        if isempty(explicitMask); error('Can''t load  explicit mask'); end

        if w.contrast_only == false
            %==============================================================%
            %  fMRI model specification
            %==============================================================%

            clear matlabbatch;

            matlabbatch{stage}.spm.stats.fmri_spec.dir =  {w.firstDir};
            matlabbatch{stage}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{stage}.spm.stats.fmri_spec.timing.RT = w.TR;
            matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t0 = 1;    % First time bin because reference for slice timing was T=0ms

            w.cond_list = [];
            for j=1:numel(w.sessions) % Sessions loop

                w.funcPathSession = fullfile (w.funcPath, w.sessions{j});

                if w.with_physIO == false

                    % Get the head movements file (store in a matrix)
                    rpF = spm_select('FPList',  fullfile(w.funcPathSession), ['^rp_' '.*\.txt$']);
                    R = load(deblank(rpF));

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
                onset_file = cellstr(spm_select('FPList', w.stimPath, ['.*run_' num2str(j, '%02.f') '_' onsets_type '.mat$']));
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

                % Get EPI smoothed images
                EPI = spm_select('ExtFPList',  fullfile(w.funcPathSession), ['^swua' '.*\.nii$'], Inf);

                % Fill SPM batch (sessions)
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).scans = cellstr(EPI);
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).multi = onset_file;
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).regress = struct('name', {}, 'val', {});
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).multi_reg = confounds_file;
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).hpf = w.hpf;

            end

            matlabbatch{stage}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{stage}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
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
            2
            %==============================================================%
            %  Model Estimation
            %==============================================================%

            matlabbatch{stage}.spm.stats.fmri_est.spmmat(1) = cellstr(fullfile(w.firstDir,'SPM.mat'));
            matlabbatch{stage}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{stage}.spm.stats.fmri_est.method.Classical = 1;

            stage = stage +1;

        end

        %==============================================================%
        %  Contrast manager
        %==============================================================%

        %========== function to 'normalize' contrast
        normcon = @(x) transpose(x/numel(find(x>0)));
        normcon_sumtoone = @(x) x./sum(x);

        %========== Retrieve conditions & prepare contrast weights
        load('cond_list.mat')
        nconds = numel(cond_list);
        cond_list = strrep(cond_list, 'Pairs+Alt_bis', 'Pairs+Altbis');  % change name to avoid issues...
        rp_conds = double(strcmp(cond_list,'rp')');

        %     % SHOW DESIGN...
        %     load('SPM.mat')
        %     tmp=find(~rp_conds);
        %     toshow = tmp(1:(3*4+2));
        %     figure;plot(SPM.xX.X(:,toshow),'Linewidth',2)
        %     legend(strrep(cond_list(toshow),'_','-'))
        %     xlim([1 360])

        %====== Create 0/1 vectors for all possible conditions
        v = create_conditions_vectors(cond_list);

        %==========  SPM batch
        matlabbatch{stage}.spm.stats.con.spmmat = cellstr(fullfile(w.firstDir,'SPM.mat'));

        %==========  Contrasts
        n_cont = 1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- Effects of interest F-contrasts (mainly for plots)
        task_conds   = eye(nconds); % All conditions
        task_conds(find(rp_conds==1),:)=[]; % Removing realignment parameters columns
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_CONDITIONS';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = task_conds;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- T-contrasts
        labels = {
            % Habituation (betas)
            'Hab_ALL'
            'Hab_Repeat'
            'Hab_Alter'
            'Hab_Pairs'
            'Hab_Quad'
            'Hab_PairsAlt'
            'Hab_Shrink'
            'Hab_PairsAltBis'
            'Hab_ThreeTwo'
            'Hab_CenterMir'
            'Hab_Complex'
            % Standard (betas)
            'Stand_ALL'
            'Stand_Repeat'
            'Stand_Alter'
            'Stand_Pairs'
            'Stand_Quad'
            'Stand_PairsAlt'
            'Stand_Shrink'
            'Stand_PairsAltBis'
            'Stand_ThreeTwo'
            'Stand_CenterMir'
            'Stand_Complex'
            % Deviants (betas)
            'Dev_ALL'
            'Dev_Repeat'
            'Dev_Alter'
            'Dev_Pairs'
            'Dev_Quad'
            'Dev_PairsAlt'
            'Dev_Shrink'
            'Dev_PairsAltBis'
            'Dev_ThreeTwo'
            'Dev_CenterMir'
            'Dev_Complex'
            % AB_Habituation (betas)
            'AB_Hab_ALL'
            'AB_Hab_Repeat'
            'AB_Hab_Alter'
            'AB_Hab_Pairs'
            'AB_Hab_Quad'
            'AB_Hab_PairsAlt'
            'AB_Hab_Shrink'
            'AB_Hab_PairsAltBis'
            'AB_Hab_ThreeTwo'
            'AB_Hab_CenterMir'
            'AB_Hab_Complex'
            % AB_Standard (betas)
            'AB_Stand_ALL'
            'AB_Stand_Repeat'
            'AB_Stand_Alter'
            'AB_Stand_Pairs'
            'AB_Stand_Quad'
            'AB_Stand_PairsAlt'
            'AB_Stand_Shrink'
            'AB_Stand_PairsAltBis'
            'AB_Stand_ThreeTwo'
            'AB_Stand_CenterMir'
            'AB_Stand_Complex'
            % AB_Deviants (betas)
            'AB_Dev_ALL'
            'AB_Dev_Repeat'
            'AB_Dev_Alter'
            'AB_Dev_Pairs'
            'AB_Dev_Quad'
            'AB_Dev_PairsAlt'
            'AB_Dev_Shrink'
            'AB_Dev_PairsAltBis'
            'AB_Dev_ThreeTwo'
            'AB_Dev_CenterMir'
            'AB_Dev_Complex'
            % BA_Habituation (betas)
            'BA_Hab_ALL'
            'BA_Hab_Repeat'
            'BA_Hab_Alter'
            'BA_Hab_Pairs'
            'BA_Hab_Quad'
            'BA_Hab_PairsAlt'
            'BA_Hab_Shrink'
            'BA_Hab_PairsAltBis'
            'BA_Hab_ThreeTwo'
            'BA_Hab_CenterMir'
            'BA_Hab_Complex'
            % BA_Standard (betas)
            'BA_Stand_ALL'
            'BA_Stand_Repeat'
            'BA_Stand_Alter'
            'BA_Stand_Pairs'
            'BA_Stand_Quad'
            'BA_Stand_PairsAlt'
            'BA_Stand_Shrink'
            'BA_Stand_PairsAltBis'
            'BA_Stand_ThreeTwo'
            'BA_Stand_CenterMir'
            'BA_Stand_Complex'
            % BA_Deviants (betas)
            'BA_Dev_ALL'
            'BA_Dev_Repeat'
            'BA_Dev_Alter'
            'BA_Dev_Pairs'
            'BA_Dev_Quad'
            'BA_Dev_PairsAlt'
            'BA_Dev_Shrink'
            'BA_Dev_PairsAltBis'
            'BA_Dev_ThreeTwo'
            'BA_Dev_CenterMir'
            'BA_Dev_Complex'
            % Other
            'Message' % Instructions
            };

        convalues =  [
            % Habituation (betas)
            normcon_sumtoone(v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_Shrink+v.Hab_PairsAltBis+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex);
            v.Hab_Repeat;
            v.Hab_Alter;
            v.Hab_Pairs;
            v.Hab_Quad;
            v.Hab_PairsAlt;
            v.Hab_Shrink;
            v.Hab_PairsAltBis;
            v.Hab_ThreeTwo;
            v.Hab_CenterMir;
            v.Hab_Complex;
            % Standard (betas)
            normcon_sumtoone(v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex);
            v.Stand_Repeat;
            v.Stand_Alter;
            v.Stand_Pairs;
            v.Stand_Quad;
            v.Stand_PairsAlt;
            v.Stand_Shrink;
            v.Stand_PairsAltBis;
            v.Stand_ThreeTwo;
            v.Stand_CenterMir;
            v.Stand_Complex;
            % Deviants (betas)
            normcon_sumtoone(v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex);
            v.Dev_Repeat;
            v.Dev_Alter;
            v.Dev_Pairs;
            v.Dev_Quad;
            v.Dev_PairsAlt;
            v.Dev_Shrink;
            v.Dev_PairsAltBis;
            v.Dev_ThreeTwo;
            v.Dev_CenterMir;
            v.Dev_Complex;
            % AB_Habituation (betas)
            normcon_sumtoone(v.AB_Hab_Repeat+v.AB_Hab_Alter+v.AB_Hab_Pairs+v.AB_Hab_Quad+v.AB_Hab_PairsAlt+v.AB_Hab_Shrink+v.AB_Hab_PairsAltBis+v.AB_Hab_ThreeTwo+v.AB_Hab_CenterMir+v.AB_Hab_Complex);
            v.Hab_Repeat;
            v.Hab_Alter;
            v.Hab_Pairs;
            v.Hab_Quad;
            v.Hab_PairsAlt;
            v.Hab_Shrink;
            v.Hab_PairsAltBis;
            v.Hab_ThreeTwo;
            v.Hab_CenterMir;
            v.Hab_Complex;
            % AB_Standard (betas)
            normcon_sumtoone(v.AB_Stand_Repeat+v.AB_Stand_Alter+v.AB_Stand_Pairs+v.AB_Stand_Quad+v.AB_Stand_PairsAlt+v.AB_Stand_PairsAltBis+v.AB_Stand_Shrink+v.AB_Stand_ThreeTwo+v.AB_Stand_CenterMir+v.AB_Stand_Complex);
            v.AB_Stand_Repeat;
            v.AB_Stand_Alter;
            v.AB_Stand_Pairs;
            v.AB_Stand_Quad;
            v.AB_Stand_PairsAlt;
            v.AB_Stand_Shrink;
            v.AB_Stand_PairsAltBis;
            v.AB_Stand_ThreeTwo;
            v.AB_Stand_CenterMir;
            v.AB_Stand_Complex;
            % AB_Deviants (betas)
            normcon_sumtoone(v.AB_Dev_Repeat+v.AB_Dev_Alter+v.AB_Dev_Pairs+v.AB_Dev_Quad+v.AB_Dev_PairsAlt+v.AB_Dev_PairsAltBis+v.AB_Dev_Shrink+v.AB_Dev_ThreeTwo+v.AB_Dev_CenterMir+v.AB_Dev_Complex);
            v.AB_Dev_Repeat;
            v.AB_Dev_Alter;
            v.AB_Dev_Pairs;
            v.AB_Dev_Quad;
            v.AB_Dev_PairsAlt;
            v.AB_Dev_Shrink;
            v.AB_Dev_PairsAltBis;
            v.AB_Dev_ThreeTwo;
            v.AB_Dev_CenterMir;
            v.AB_Dev_Complex;
            % BA_Habituation (betas)
            normcon_sumtoone(v.BA_Hab_Repeat+v.BA_Hab_Alter+v.BA_Hab_Pairs+v.BA_Hab_Quad+v.BA_Hab_PairsAlt+v.BA_Hab_Shrink+v.BA_Hab_PairsAltBis+v.BA_Hab_ThreeTwo+v.BA_Hab_CenterMir+v.BA_Hab_Complex);
            v.BA_Hab_Repeat;
            v.BA_Hab_Alter;
            v.BA_Hab_Pairs;
            v.BA_Hab_Quad;
            v.BA_Hab_PairsAlt;
            v.BA_Hab_Shrink;
            v.BA_Hab_PairsAltBis;
            v.BA_Hab_ThreeTwo;
            v.BA_Hab_CenterMir;
            v.BA_Hab_Complex;
            % BA_Standard (betas)
            normcon_sumtoone(v.BA_Stand_Repeat+v.BA_Stand_Alter+v.BA_Stand_Pairs+v.BA_Stand_Quad+v.BA_Stand_PairsAlt+v.BA_Stand_PairsAltBis+v.BA_Stand_Shrink+v.BA_Stand_ThreeTwo+v.BA_Stand_CenterMir+v.BA_Stand_Complex);
            v.BA_Stand_Repeat;
            v.BA_Stand_Alter;
            v.BA_Stand_Pairs;
            v.BA_Stand_Quad;
            v.BA_Stand_PairsAlt;
            v.BA_Stand_Shrink;
            v.BA_Stand_PairsAltBis;
            v.BA_Stand_ThreeTwo;
            v.BA_Stand_CenterMir;
            v.BA_Stand_Complex;
            % BA_Deviants (betas)
            normcon_sumtoone(v.BA_Dev_Repeat+v.BA_Dev_Alter+v.BA_Dev_Pairs+v.BA_Dev_Quad+v.BA_Dev_PairsAlt+v.BA_Dev_PairsAltBis+v.BA_Dev_Shrink+v.BA_Dev_ThreeTwo+v.BA_Dev_CenterMir+v.BA_Dev_Complex);
            v.BA_Dev_Repeat;
            v.BA_Dev_Alter;
            v.BA_Dev_Pairs;
            v.BA_Dev_Quad;
            v.BA_Dev_PairsAlt;
            v.BA_Dev_Shrink;
            v.BA_Dev_PairsAltBis;
            v.BA_Dev_ThreeTwo;
            v.BA_Dev_CenterMir;
            v.BA_Dev_Complex;
            % Other
            v.Message;
            ];

        %         if sum(sum(convalues,2)) ~=0; error('Issue with contrast weights...'); end
        %         figure;imagesc(convalues)

        % Create spm contrasts
        for ii = 1:numel(labels)
            %-------%
            matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = labels{ii};
            matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = convalues(ii,:);
            matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none';
            n_cont = n_cont+1;
            %-------%
        end

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- ADDITIONAL F-contrasts (mainly for plots)
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_MAIN_CONDITIONS';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = convalues([2:11;13:22;24:33],:);
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- Effects of interest F-contrasts v2,
        %WITHOUT DEV  %%
        task_conds   = eye(nconds); % All conditions
        task_conds([find(rp_conds==1) find(v.dev_betas)' find(v.Message)],:)=[]; % Removing realignment parameters columns & deviant & message columns
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_CONDITIONS_v2';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = task_conds;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- ADDITIONAL T-contrasts (complexity)
        vectors = [ 1 2 3 4 5 6 7 8 9 10;
            4 6 6 6 12 14 13 15 17 23;
            4 6 6 6 12 15 16 18 21 28;
            4 6 6 6 12 15 16 18 21 4;
            0 1 .47 .20 .73 .47 .73 .47 .47 .47;
            16 8 4 2 2 1 2 2 1 1;
            1 2 4 8 8 16 8 8 16 16];
        % add quadratic of geochunk
        y = vectors(3,:);
        A = y;
        B = y.^2;
        B = B - mean(B);
        C = A - (sum(A.*B)./sum(B.^2).*B);
        vectors = [vectors; C];
        vnames  = { 'BasicComplexity';
            'GeoComplexity'
            'GeoChunkComplexity'
            'GeoChunkCollapse'
            'pAlt'
            'Periodicity'
            'Period'
            'GeoChunkQuadra'};
        for nn = 1 :size(vectors,1)
            for ttype = {'Hab', 'Stand', 'Dev', 'AB_Hab', 'AB_Stand', 'AB_Dev', 'BA_Hab', 'BA_Stand', 'BA_Dev'}
                vector = []; vector = vectors(nn,:) - mean(vectors(nn,:) );
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = ['T-' char(ttype) '-' vnames{nn}];
                weights =  [(v.([char(ttype) '_Repeat']))     *vector(1)+...
                    (v.([char(ttype) '_Alter']))      *vector(2)+...
                    (v.([char(ttype) '_Pairs']))      *vector(3)+...
                    (v.([char(ttype) '_Quad']))       *vector(4)+...
                    (v.([char(ttype) '_PairsAlt']))   *vector(5)+...
                    (v.([char(ttype) '_Shrink']))     *vector(6)+...
                    (v.([char(ttype) '_PairsAltBis']))*vector(7)+...
                    (v.([char(ttype) '_ThreeTwo']))   *vector(8)+...
                    (v.([char(ttype) '_CenterMir']))  *vector(9)+...
                    (v.([char(ttype) '_Complex']))    *vector(10)];
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = normcon2(weights, true)';
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none'; n_cont = n_cont+1;
            end
        end

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- ADDITIONAL T-contrasts (negative complexityvalues)
        vectors = [ 1 2 3 4 5 6 7 8 9 10;
            4 6 6 6 12 14 13 15 17 23;
            4 6 6 6 12 15 16 18 21 28;
            4 6 6 6 12 15 16 18 21 4;
            0 1 .47 .20 .73 .47 .73 .47 .47 .47;
            16 8 4 2 2 1 2 2 1 1;
            1 2 4 8 8 16 8 8 16 16];
        % add quadratic of geochunk
        y = vectors(3,:);
        A = y;
        B = y.^2;
        B = B - mean(B);
        C = A - (sum(A.*B)./sum(B.^2).*B);
        vectors = [vectors; C];
        vnames  = { 'BasicComplexity';
            'GeoComplexity'
            'GeoChunkComplexity'
            'GeoChunkCollapse'
            'pAlt'
            'Periodicity'
            'Period'
            'GeoChunkQuadra'};
        vectors = -vectors;
        for nn = 1 :size(vectors,1)
            for ttype = {'Hab', 'Stand', 'Dev', 'AB_Hab', 'AB_Stand', 'AB_Dev', 'BA_Hab', 'BA_Stand', 'BA_Dev'}
                vector = []; vector = vectors(nn,:) - mean(vectors(nn,:) );
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = ['T-' char(ttype) '-' vnames{nn} '-negat'];
                weights =  [(v.([char(ttype) '_Repeat']))     *vector(1)+...
                    (v.([char(ttype) '_Alter']))      *vector(2)+...
                    (v.([char(ttype) '_Pairs']))      *vector(3)+...
                    (v.([char(ttype) '_Quad']))       *vector(4)+...
                    (v.([char(ttype) '_PairsAlt']))   *vector(5)+...
                    (v.([char(ttype) '_Shrink']))     *vector(6)+...
                    (v.([char(ttype) '_PairsAltBis']))*vector(7)+...
                    (v.([char(ttype) '_ThreeTwo']))   *vector(8)+...
                    (v.([char(ttype) '_CenterMir']))  *vector(9)+...
                    (v.([char(ttype) '_Complex']))    *vector(10)];
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = normcon2(weights, true)';
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none'; n_cont = n_cont+1;
            end
        end

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- ADDITIONAL T-contrasts (partial complexity)
        vector = [4 6 6 6 12 15 16 18 21 28];
        vname  = 'GeoChunkComplexity';
        for nn = 1:6
            for ttype = {'Hab', 'Stand', 'Dev', 'AB_Hab', 'AB_Stand', 'AB_Dev', 'BA_Hab', 'BA_Stand', 'BA_Dev'}
                vectorpart = vector(1:4+nn);
                vectorpart = vectorpart- mean(vectorpart);
                vectorpart = [vectorpart zeros(1, 6-nn)];
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = ['T-' char(ttype) '-' vname '_part' num2str(nn)];
                weights =  [(v.([char(ttype) '_Repeat']))     *vectorpart(1)+...
                    (v.([char(ttype) '_Alter']))      *vectorpart(2)+...
                    (v.([char(ttype) '_Pairs']))      *vectorpart(3)+...
                    (v.([char(ttype) '_Quad']))       *vectorpart(4)+...
                    (v.([char(ttype) '_PairsAlt']))   *vectorpart(5)+...
                    (v.([char(ttype) '_Shrink']))     *vectorpart(7)+...
                    (v.([char(ttype) '_PairsAltBis']))*vectorpart(6)+...
                    (v.([char(ttype) '_ThreeTwo']))   *vectorpart(8)+...
                    (v.([char(ttype) '_CenterMir']))  *vectorpart(9)+...
                    (v.([char(ttype) '_Complex']))    *vectorpart(10)];
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = normcon2(weights, true)';
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none'; n_cont = n_cont+1;
            end
        end

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- CORRELATION WITH PERFORMANCE
        %         d = load(fullfile (w.datadir, 'BehavioralData_Results', 'AllSubj_BehavioralData.mat'));
        d = load(fullfile (w.datadir, 'BehavioralData_Results', 'AllSubj_BehavioralData_ABBA.mat'));
        subj_data = d.subj_data_all(strcmp(string(d.subj_data_all.Subject), w.subjects{iS}),:);
        subj_data = sortrows(subj_data,'seqID','ascend');
        subj_data = sortrows(subj_data,'version','ascend');
        if size(subj_data,1) ~= 20; error('Missing some behavioral data!'); end
        vectors = [subj_data.MissRate'; subj_data.mean_RT'; subj_data.LISAS'; subj_data.nFA'];
        vnames  = { 'MissRate';
            'mean_RT'
            'LISAS'
            'nFA'};
        for nn = 1 :size(vectors,1)
            for ttype = {'Hab', 'Stand', 'Dev'}
                vector = []; vector = vectors(nn,:) - nanmean(vectors(nn,:) );
                vector(isnan(vector)) = 0;
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = ['T-' char(ttype) '-' vnames{nn}];
                weights =  [(v.(['AB_' char(ttype) '_Repeat']))     *vector(1)+...
                    (v.(['AB_' char(ttype) '_Alter']))      *vector(2)+...
                    (v.(['AB_' char(ttype) '_Pairs']))      *vector(3)+...
                    (v.(['AB_' char(ttype) '_Quad']))       *vector(4)+...
                    (v.(['AB_' char(ttype) '_PairsAlt']))   *vector(5)+...
                    (v.(['AB_' char(ttype) '_Shrink']))     *vector(6)+...
                    (v.(['AB_' char(ttype) '_PairsAltBis']))*vector(7)+...
                    (v.(['AB_' char(ttype) '_ThreeTwo']))   *vector(8)+...
                    (v.(['AB_' char(ttype) '_CenterMir']))  *vector(9)+...
                    (v.(['AB_' char(ttype) '_Complex']))    *vector(10)+...
                    (v.(['BA_' char(ttype) '_Repeat']))     *vector(11)+...
                    (v.(['BA_' char(ttype) '_Alter']))      *vector(12)+...
                    (v.(['BA_' char(ttype) '_Pairs']))      *vector(13)+...
                    (v.(['BA_' char(ttype) '_Quad']))       *vector(14)+...
                    (v.(['BA_' char(ttype) '_PairsAlt']))   *vector(15)+...
                    (v.(['BA_' char(ttype) '_Shrink']))     *vector(16)+...
                    (v.(['BA_' char(ttype) '_PairsAltBis']))*vector(17)+...
                    (v.(['BA_' char(ttype) '_ThreeTwo']))   *vector(18)+...
                    (v.(['BA_' char(ttype) '_CenterMir']))  *vector(19)+...
                    (v.(['BA_' char(ttype) '_Complex']))    *vector(20)
                    ];
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = normcon2(weights, true)';
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none'; n_cont = n_cont+1;
            end
        end

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- F tests between sequences
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-Test_allseq_HAB';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = [v.Hab_Repeat;v.Hab_Alter;v.Hab_Pairs;v.Hab_Quad;v.Hab_PairsAlt;v.Hab_PairsAltBis;v.Hab_Shrink;v.Hab_ThreeTwo;v.Hab_CenterMir;v.Hab_Complex];
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-Test_allseq_STAND';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = [v.Stand_Repeat;v.Stand_Alter;v.Stand_Pairs;v.Stand_Quad;v.Stand_PairsAlt;v.Stand_PairsAltBis;v.Stand_Shrink;v.Stand_ThreeTwo;v.Stand_CenterMir;v.Stand_Complex];
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-Test_allseq_DEV';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = [v.Dev_Repeat;v.Dev_Alter;v.Dev_Pairs;v.Dev_Quad;v.Dev_PairsAlt;v.Dev_PairsAltBis;v.Dev_Shrink;v.Dev_ThreeTwo;v.Dev_CenterMir;v.Dev_Complex];
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-Test_allseqVS_HAB';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = [   v.Hab_Repeat - (v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_Alter - (v.Hab_Repeat+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_Pairs - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_Quad - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_PairsAlt - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_Shrink - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_PairsAltBis - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_ThreeTwo - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_CenterMir+v.Hab_Complex)/9;...
            v.Hab_CenterMir - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_Complex)/9;...
            v.Hab_Complex - (v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_PairsAltBis+v.Hab_Shrink+v.Hab_ThreeTwo+v.Hab_CenterMir)/9];
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-Test_allseqVS_STAND';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = [   v.Stand_Repeat - (v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_Alter - (v.Stand_Repeat+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_Pairs - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_Quad - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_PairsAlt - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_Shrink - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_PairsAltBis - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_ThreeTwo - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_CenterMir+v.Stand_Complex)/9;...
            v.Stand_CenterMir - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_Complex)/9;...
            v.Stand_Complex - (v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir)/9];
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-Test_allseqVS_DEV';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = [   v.Dev_Repeat - (v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_Alter - (v.Dev_Repeat+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_Pairs - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_Quad - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_PairsAlt - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_Shrink - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_PairsAltBis - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_ThreeTwo - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_CenterMir+v.Dev_Complex)/9;...
            v.Dev_CenterMir - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_Complex)/9;...
            v.Dev_Complex - (v.Dev_Repeat+v.Dev_Alter+v.Dev_Pairs+v.Dev_Quad+v.Dev_PairsAlt+v.Dev_PairsAltBis+v.Dev_Shrink+v.Dev_ThreeTwo+v.Dev_CenterMir)/9];
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;


        disp(['Number of contrasts = ' num2str(n_cont-1)])
        %-------%
        matlabbatch{stage}.spm.stats.con.delete = 1;

        %%
        if w.contrast_only == false
            save(fullfile(w.firstDir, 'SPM12_1stLevel_matlabbatch.mat'),'matlabbatch');
        end
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

        tmp = load(fullfile(w.firstDir, 'SPM.mat'));
        labels = {tmp.SPM.xCon.name}.';
        save(fullfile(w.firstDir, 'Contrast_labels.mat'),'labels');

    end

%% Model 2: one subject full analysis, with onsets VERSION 1: sequence*(hab+stand+1dev) - with DetectedDev only
    function DoFirstLevel_v2(w, iS)

        onsets_type = 'onsets_v10';
        which_subfolder = 'Model_10';

        fprintf([' \n \n']);
        fprintf('========================================================================\n');
        fprintf(['  ' w.subjects{iS} ': running Model 1: one subject full analysis using ' onsets_type '...\n']);
        fprintf('========================================================================\n');


        %% Output folder
        if w.with_physIO
            model_folder = fullfile(which_subfolder, 'PhysIO');
        else
            model_folder = fullfile(which_subfolder, 'Standard');
        end
        w.firstDir = fullfile (w.datadir, '3_FIRST_LEVEL',  'ABseq', w.subjects{iS}, model_folder);
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

        % Get explicit mask
        explicitMask = spm_select('FPList', w.T1Path, '^brexplicitMask_wc1wc2wc3');
        if isempty(explicitMask); error('Can''t load  explicit mask'); end

        if w.contrast_only == false
            %==============================================================%
            %  fMRI model specification
            %==============================================================%

            clear matlabbatch;

            matlabbatch{stage}.spm.stats.fmri_spec.dir =  {w.firstDir};
            matlabbatch{stage}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{stage}.spm.stats.fmri_spec.timing.RT = w.TR;
            matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t0 = 1;    % First time bin because reference for slice timing was T=0ms

            w.cond_list = [];
            for j=1:numel(w.sessions) % Sessions loop

                w.funcPathSession = fullfile (w.funcPath, w.sessions{j});

                if w.with_physIO == false

                    % Get the head movements file (store in a matrix)
                    rpF = spm_select('FPList',  fullfile(w.funcPathSession), ['^rp_' '.*\.txt$']);
                    R = load(deblank(rpF));

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
                onset_file = cellstr(spm_select('FPList', w.stimPath, ['.*run_' num2str(j, '%02.f') '_' onsets_type '.mat$']));
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

                % Get EPI smoothed images
                EPI = spm_select('ExtFPList',  fullfile(w.funcPathSession), ['^swua' '.*\.nii$'], Inf);

                % Fill SPM batch (sessions)
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).scans = cellstr(EPI);
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).multi = onset_file;
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).regress = struct('name', {}, 'val', {});
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).multi_reg = confounds_file;
                matlabbatch{stage}.spm.stats.fmri_spec.sess(j).hpf = w.hpf;

            end


            matlabbatch{stage}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{stage}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
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

            matlabbatch{stage}.spm.stats.fmri_est.spmmat(1) = cellstr(fullfile(w.firstDir,'SPM.mat'));
            matlabbatch{stage}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{stage}.spm.stats.fmri_est.method.Classical = 1;

            stage = stage +1;

        end

        %==============================================================%
        %  Contrast manager
        %==============================================================%

        %========== function to 'normalize' contrast
        normcon = @(x) transpose(x/numel(find(x>0)));
        normcon_sumtoone = @(x) x./sum(x,'omitnan');

        %========== Retrieve conditions & prepare contrast weights
        load('cond_list.mat')
        nconds = numel(cond_list);
        cond_list = strrep(cond_list, 'Pairs+Alt_bis', 'Pairs+Altbis');  % change name to avoid issues...
        rp_conds = double(strcmp(cond_list,'rp')');

        %     % SHOW DESIGN...
        %     load('SPM.mat')
        %     tmp=find(~rp_conds);
        %     toshow = tmp(1:(3*4+2));
        %     figure;plot(SPM.xX.X(:,toshow),'Linewidth',2)
        %     legend(strrep(cond_list(toshow),'_','-'))
        %     xlim([1 360])

        %====== Create 0/1 vectors for all possible conditions
        v = create_conditions_vectors(cond_list);

        %==========  SPM batch
        matlabbatch{stage}.spm.stats.con.spmmat = cellstr(fullfile(w.firstDir,'SPM.mat'));

        %==========  Contrasts
        n_cont = 1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- Effects of interest F-contrasts (mainly for plots)
        task_conds   = eye(nconds); % All conditions
        task_conds(find(rp_conds==1),:)=[]; % Removing realignment parameters columns
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_CONDITIONS';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = task_conds;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- T-contrasts
        labels = {
            % Habituation (betas)
            'Hab_ALL'
            'Hab_Repeat'
            'Hab_Alter'
            'Hab_Pairs'
            'Hab_Quad'
            'Hab_PairsAlt'
            'Hab_Shrink'
            'Hab_PairsAltBis'
            'Hab_ThreeTwo'
            'Hab_CenterMir'
            'Hab_Complex'
            % Standard (betas)
            'Stand_ALL'
            'Stand_Repeat'
            'Stand_Alter'
            'Stand_Pairs'
            'Stand_Quad'
            'Stand_PairsAlt'
            'Stand_Shrink'
            'Stand_PairsAltBis'
            'Stand_ThreeTwo'
            'Stand_CenterMir'
            'Stand_Complex'
            % Deviants (betas)
            'DetectedDev_ALL'
            'DetectedDev_Repeat'
            'DetectedDev_Alter'
            'DetectedDev_Pairs'
            'DetectedDev_Quad'
            'DetectedDev_PairsAlt'
            'DetectedDev_Shrink'
            'DetectedDev_PairsAltBis'
            'DetectedDev_ThreeTwo'
            'DetectedDev_CenterMir'
            'DetectedDev_Complex'
            % AB_Habituation (betas)
            'AB_Hab_ALL'
            'AB_Hab_Repeat'
            'AB_Hab_Alter'
            'AB_Hab_Pairs'
            'AB_Hab_Quad'
            'AB_Hab_PairsAlt'
            'AB_Hab_Shrink'
            'AB_Hab_PairsAltBis'
            'AB_Hab_ThreeTwo'
            'AB_Hab_CenterMir'
            'AB_Hab_Complex'
            % AB_Standard (betas)
            'AB_Stand_ALL'
            'AB_Stand_Repeat'
            'AB_Stand_Alter'
            'AB_Stand_Pairs'
            'AB_Stand_Quad'
            'AB_Stand_PairsAlt'
            'AB_Stand_Shrink'
            'AB_Stand_PairsAltBis'
            'AB_Stand_ThreeTwo'
            'AB_Stand_CenterMir'
            'AB_Stand_Complex'
            % AB_Deviants (betas)
            'AB_DetectedDev_ALL'
            'AB_DetectedDev_Repeat'
            'AB_DetectedDev_Alter'
            'AB_DetectedDev_Pairs'
            'AB_DetectedDev_Quad'
            'AB_DetectedDev_PairsAlt'
            'AB_DetectedDev_Shrink'
            'AB_DetectedDev_PairsAltBis'
            'AB_DetectedDev_ThreeTwo'
            'AB_DetectedDev_CenterMir'
            'AB_DetectedDev_Complex'
            % BA_Habituation (betas)
            'BA_Hab_ALL'
            'BA_Hab_Repeat'
            'BA_Hab_Alter'
            'BA_Hab_Pairs'
            'BA_Hab_Quad'
            'BA_Hab_PairsAlt'
            'BA_Hab_Shrink'
            'BA_Hab_PairsAltBis'
            'BA_Hab_ThreeTwo'
            'BA_Hab_CenterMir'
            'BA_Hab_Complex'
            % BA_Standard (betas)
            'BA_Stand_ALL'
            'BA_Stand_Repeat'
            'BA_Stand_Alter'
            'BA_Stand_Pairs'
            'BA_Stand_Quad'
            'BA_Stand_PairsAlt'
            'BA_Stand_Shrink'
            'BA_Stand_PairsAltBis'
            'BA_Stand_ThreeTwo'
            'BA_Stand_CenterMir'
            'BA_Stand_Complex'
            % BA_Deviants (betas)
            'BA_DetectedDev_ALL'
            'BA_DetectedDev_Repeat'
            'BA_DetectedDev_Alter'
            'BA_DetectedDev_Pairs'
            'BA_DetectedDev_Quad'
            'BA_DetectedDev_PairsAlt'
            'BA_DetectedDev_Shrink'
            'BA_DetectedDev_PairsAltBis'
            'BA_DetectedDev_ThreeTwo'
            'BA_DetectedDev_CenterMir'
            'BA_DetectedDev_Complex'
            % Other
            'Message' % Instructions
            };

        convalues =  [
            % Habituation (betas)
            normcon_sumtoone(v.Hab_Repeat+v.Hab_Alter+v.Hab_Pairs+v.Hab_Quad+v.Hab_PairsAlt+v.Hab_Shrink+v.Hab_PairsAltBis+v.Hab_ThreeTwo+v.Hab_CenterMir+v.Hab_Complex);
            v.Hab_Repeat;
            v.Hab_Alter;
            v.Hab_Pairs;
            v.Hab_Quad;
            v.Hab_PairsAlt;
            v.Hab_Shrink;
            v.Hab_PairsAltBis;
            v.Hab_ThreeTwo;
            v.Hab_CenterMir;
            v.Hab_Complex;
            % Standard (betas)
            normcon_sumtoone(v.Stand_Repeat+v.Stand_Alter+v.Stand_Pairs+v.Stand_Quad+v.Stand_PairsAlt+v.Stand_PairsAltBis+v.Stand_Shrink+v.Stand_ThreeTwo+v.Stand_CenterMir+v.Stand_Complex);
            v.Stand_Repeat;
            v.Stand_Alter;
            v.Stand_Pairs;
            v.Stand_Quad;
            v.Stand_PairsAlt;
            v.Stand_Shrink;
            v.Stand_PairsAltBis;
            v.Stand_ThreeTwo;
            v.Stand_CenterMir;
            v.Stand_Complex;
            % Deviants (betas)
            %             normcon_sumtoone(sum([v.DetectedDev_Repeat;v.DetectedDev_Alter;v.DetectedDev_Pairs;v.DetectedDev_Quad;v.DetectedDev_PairsAlt;v.DetectedDev_PairsAltBis;v.DetectedDev_Shrink;v.DetectedDev_ThreeTwo;v.DetectedDev_CenterMir;v.DetectedDev_Complex], 'omitnan'));
            normcon_sumtoone(v.detecteddev_betas')
            v.DetectedDev_Repeat;
            v.DetectedDev_Alter;
            v.DetectedDev_Pairs;
            v.DetectedDev_Quad;
            v.DetectedDev_PairsAlt;
            v.DetectedDev_Shrink;
            v.DetectedDev_PairsAltBis;
            v.DetectedDev_ThreeTwo;
            v.DetectedDev_CenterMir;
            v.DetectedDev_Complex;
            % AB_Habituation (betas)
            normcon_sumtoone(v.AB_Hab_Repeat+v.AB_Hab_Alter+v.AB_Hab_Pairs+v.AB_Hab_Quad+v.AB_Hab_PairsAlt+v.AB_Hab_Shrink+v.AB_Hab_PairsAltBis+v.AB_Hab_ThreeTwo+v.AB_Hab_CenterMir+v.AB_Hab_Complex);
            v.Hab_Repeat;
            v.Hab_Alter;
            v.Hab_Pairs;
            v.Hab_Quad;
            v.Hab_PairsAlt;
            v.Hab_Shrink;
            v.Hab_PairsAltBis;
            v.Hab_ThreeTwo;
            v.Hab_CenterMir;
            v.Hab_Complex;
            % AB_Standard (betas)
            normcon_sumtoone(v.AB_Stand_Repeat+v.AB_Stand_Alter+v.AB_Stand_Pairs+v.AB_Stand_Quad+v.AB_Stand_PairsAlt+v.AB_Stand_PairsAltBis+v.AB_Stand_Shrink+v.AB_Stand_ThreeTwo+v.AB_Stand_CenterMir+v.AB_Stand_Complex);
            v.AB_Stand_Repeat;
            v.AB_Stand_Alter;
            v.AB_Stand_Pairs;
            v.AB_Stand_Quad;
            v.AB_Stand_PairsAlt;
            v.AB_Stand_Shrink;
            v.AB_Stand_PairsAltBis;
            v.AB_Stand_ThreeTwo;
            v.AB_Stand_CenterMir;
            v.AB_Stand_Complex;
            % AB_Deviants (betas)
            normcon_sumtoone(sum([v.AB_DetectedDev_Repeat;v.AB_DetectedDev_Alter;v.AB_DetectedDev_Pairs;v.AB_DetectedDev_Quad;v.AB_DetectedDev_PairsAlt;v.AB_DetectedDev_PairsAltBis;v.AB_DetectedDev_Shrink;v.AB_DetectedDev_ThreeTwo;v.AB_DetectedDev_CenterMir;v.AB_DetectedDev_Complex], 'omitnan'));
            v.AB_DetectedDev_Repeat;
            v.AB_DetectedDev_Alter;
            v.AB_DetectedDev_Pairs;
            v.AB_DetectedDev_Quad;
            v.AB_DetectedDev_PairsAlt;
            v.AB_DetectedDev_Shrink;
            v.AB_DetectedDev_PairsAltBis;
            v.AB_DetectedDev_ThreeTwo;
            v.AB_DetectedDev_CenterMir;
            v.AB_DetectedDev_Complex;
            % BA_Habituation (betas)
            normcon_sumtoone(v.BA_Hab_Repeat+v.BA_Hab_Alter+v.BA_Hab_Pairs+v.BA_Hab_Quad+v.BA_Hab_PairsAlt+v.BA_Hab_Shrink+v.BA_Hab_PairsAltBis+v.BA_Hab_ThreeTwo+v.BA_Hab_CenterMir+v.BA_Hab_Complex);
            v.BA_Hab_Repeat;
            v.BA_Hab_Alter;
            v.BA_Hab_Pairs;
            v.BA_Hab_Quad;
            v.BA_Hab_PairsAlt;
            v.BA_Hab_Shrink;
            v.BA_Hab_PairsAltBis;
            v.BA_Hab_ThreeTwo;
            v.BA_Hab_CenterMir;
            v.BA_Hab_Complex;
            % BA_Standard (betas)
            normcon_sumtoone(v.BA_Stand_Repeat+v.BA_Stand_Alter+v.BA_Stand_Pairs+v.BA_Stand_Quad+v.BA_Stand_PairsAlt+v.BA_Stand_PairsAltBis+v.BA_Stand_Shrink+v.BA_Stand_ThreeTwo+v.BA_Stand_CenterMir+v.BA_Stand_Complex);
            v.BA_Stand_Repeat;
            v.BA_Stand_Alter;
            v.BA_Stand_Pairs;
            v.BA_Stand_Quad;
            v.BA_Stand_PairsAlt;
            v.BA_Stand_Shrink;
            v.BA_Stand_PairsAltBis;
            v.BA_Stand_ThreeTwo;
            v.BA_Stand_CenterMir;
            v.BA_Stand_Complex;
            % BA_Deviants (betas)
            normcon_sumtoone(sum([v.BA_DetectedDev_Repeat;v.BA_DetectedDev_Alter;v.BA_DetectedDev_Pairs;v.BA_DetectedDev_Quad;v.BA_DetectedDev_PairsAlt;v.BA_DetectedDev_PairsAltBis;v.BA_DetectedDev_Shrink;v.BA_DetectedDev_ThreeTwo;v.BA_DetectedDev_CenterMir;v.BA_DetectedDev_Complex], 'omitnan'));
            v.BA_DetectedDev_Repeat;
            v.BA_DetectedDev_Alter;
            v.BA_DetectedDev_Pairs;
            v.BA_DetectedDev_Quad;
            v.BA_DetectedDev_PairsAlt;
            v.BA_DetectedDev_Shrink;
            v.BA_DetectedDev_PairsAltBis;
            v.BA_DetectedDev_ThreeTwo;
            v.BA_DetectedDev_CenterMir;
            v.BA_DetectedDev_Complex;
            % Other
            v.Message;
            ];

        %         figure;imagesc(convalues)

        % Remove "missing" conditions
        torem = find(isnan(sum(convalues,2)));
        convalues(torem, :) = [];
        labels(torem) = [];


        % Create spm contrasts
        for ii = 1:numel(labels)
            %-------%
            matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = labels{ii};
            matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = convalues(ii,:);
            matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none';
            n_cont = n_cont+1;
            %-------%
        end

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- ADDITIONAL F-contrasts (mainly for plots)
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_MAIN_CONDITIONS';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = convalues([2:11;13:22;24:33],:);
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- Effects of interest F-contrasts v2,
        %WITHOUT DEV  %%
        task_conds   = eye(nconds); % All conditions
        task_conds([find(rp_conds==1) find(v.dev_betas)' find(v.Message)],:)=[]; % Removing realignment parameters columns & deviant & message columns
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.name    = 'F-ALL_CONDITIONS_v2';
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.weights = task_conds;
        matlabbatch{stage}.spm.stats.con.consess{n_cont}.fcon.sessrep = 'none'; n_cont = n_cont+1;

        %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_- ADDITIONAL T-contrasts (complexity)
        vectors = [ 1 2 3 4 5 6 7 8 9 10;
            4 6 6 6 12 14 13 15 17 23;
            4 6 6 6 12 15 16 18 21 28;
            4 6 6 6 12 15 16 18 21 4;
            0 1 .47 .20 .73 .47 .73 .47 .47 .47;
            16 8 4 2 2 1 2 2 1 1;
            1 2 4 8 8 16 8 8 16 16];
        % add quadratic of geochunk
        y = vectors(3,:);
        A = y;
        B = y.^2;
        B = B - mean(B);
        C = A - (sum(A.*B)./sum(B.^2).*B);
        vectors = [vectors; C];
        vnames  = { 'BasicComplexity';
            'GeoComplexity'
            'GeoChunkComplexity'
            'GeoChunkCollapse'
            'pAlt'
            'Periodicity'
            'Period'
            'GeoChunkQuadra'};
        for nn = 1 :size(vectors,1)
            for ttype = {'Hab', 'Stand', 'DetectedDev', 'AB_Hab', 'AB_Stand', 'AB_DetectedDev', 'BA_Hab', 'BA_Stand', 'BA_DetectedDev'}
                vector = []; vector = vectors(nn,:) - mean(vectors(nn,:) );
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.name    = ['T-' char(ttype) '-' vnames{nn}];
                weights = [v.([char(ttype) '_Repeat']);
                    v.([char(ttype) '_Alter']);
                    v.([char(ttype) '_Pairs']);
                    v.([char(ttype) '_Quad']);
                    v.([char(ttype) '_PairsAlt']);
                    v.([char(ttype) '_Shrink']);
                    v.([char(ttype) '_PairsAltBis']);
                    v.([char(ttype) '_ThreeTwo']);
                    v.([char(ttype) '_CenterMir']);
                    v.([char(ttype) '_Complex']);];
                weights(weights>0) = 1;
                weights = weights.*vector';
                weights = sum(weights, 'omitnan');
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.weights = normcon2(weights)';
                matlabbatch{stage}.spm.stats.con.consess{n_cont}.tcon.sessrep = 'none'; n_cont = n_cont+1;
            end
        end


        disp(['Number of contrasts = ' num2str(n_cont-1)])
        %-------%
        matlabbatch{stage}.spm.stats.con.delete = 1;

        %%
        if w.contrast_only == false
            save(fullfile(w.firstDir, 'SPM12_1stLevel_matlabbatch.mat'),'matlabbatch');
        end
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

        tmp = load(fullfile(w.firstDir, 'SPM.mat'));
        labels = {tmp.SPM.xCon.name}.';
        save(fullfile(w.firstDir, 'Contrast_labels.mat'),'labels');

    end

end