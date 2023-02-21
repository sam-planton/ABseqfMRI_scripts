function Do_ABseq16fMRI_4_SPM_Preprocessing(w,iS)

cdd = pwd;
if isunix
    addpath('/home/sp253886/MatlabTools/Art')
end

%% Parameters of EPI files
tmp = load(fullfile(w.niftidir, w.subjects{iS}, 'func', 'SliceTimingInfo.mat'));
w.nSlices         = tmp.NumberOfSlices;     %
w.TR              = tmp.TR;  % Repetition time (s)
w.thickness       = tmp.SliceThickness;    % Slice thickness (mm)
w.READOUT_TIME    = tmp.total_readout_time_spm*1000 ;
w.slice_times     = tmp.SliceTiming;
w.refslice        = w.slice_times(floor(size(w.slice_times,2)/2));  % refSlice is the 34th on 69... /// this is not used anymore

%% Launch preprocessing functions

fprintf(' \n \n');
fprintf('=======================================================================\n');
fprintf(['                 ' w.subjects{iS} ' SPM preprocessing...\n']);
fprintf('=======================================================================\n');

w.sessions    =  {'RUN1' 'RUN2' 'RUN3' 'RUN4' 'RUN5' 'SynLoca'};   	% session directories (parent=functional)

w.subPath          = fullfile (w.niftidir,  w.subjects{iS});
w.funcPath         = fullfile (w.subPath, 'func');
w.T1Path           = fullfile (w.subPath, 'anat');
w.fieldmapPath     = fullfile (w.subPath, 'fieldmap');
cd(w.subPath);

% Prefix used to identify unprocessed images (subject code...)
% We take it (1st 8 characters) from the anat file the first time
% (and save it) or recover it from a text file other times
tmpname = fullfile(w.subPath, 'subPrefix.txt');
if exist(tmpname,'file')
    fid = fopen(tmpname,'r');
    w.imgPrefix = fscanf(fid,'%s');
else
    anat = spm_select('FPList', w.T1Path, ['^' '.*'  '.*\.nii$']); [~,name] = fileparts(anat);
    w.imgPrefix = name(1:8);
    fid = fopen(tmpname,'wt'); fprintf(fid, '%s', w.imgPrefix); fclose(fid);
end
disp(w)

%%======================================================%%
%%===== Do Preprocessing step by step  =================%%
%%======================================================%%
DoFieldMap(w,iS)
DoSliceTiming(w,iS);
DoRealignUnwarp(w,iS);
DoCoregister(w,iS);
DoSegment(w,iS);
DoNormalise(w,iS);
DoSurfaceRenderAndBrainExtract(w,iS)
DoSmooth(w,iS)
DoExplicitMask(w,iS)
%%======================================================%%
%%======================================================%%
%%======================================================%%

%% SPM batchs as subfunctions
    function DoFieldMap(w,iS)
        clear matlabbatch;

        fprintf(' \n \n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' Fieldmap...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Precalc fieldmap files (computed with FSL topup)
        precalcfield = spm_select('ExtFPList', w.fieldmapPath, ['^' w.subjects{iS} '_topup_fieldmap.nii']);
        magfield     = spm_select('ExtFPList', w.fieldmapPath, ['^' w.subjects{iS} '_topup_magnitude.nii'], 1:1);
        if isempty(precalcfield)
            error(['File not found !! ' w.fieldmapPath '/^' w.subjects{iS} '_topup_fieldmap.nii' ])
        end

        % Get the T1 template
        path_FielpMap = which('Fieldmap');
        [path , ~, ~] = fileparts(path_FielpMap);
        template      = fullfile(path, 'T1.nii');

        % Get T1 structural file
        anat = spm_select('FPList', w.T1Path, ['^' w.imgPrefix '.*'  '.*\.nii$']);

        % Get the first EPI file
        EPIfile = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{1}), ['^' w.imgPrefix  '.*\.nii$'], 1:1);

        % SPM batch
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = cellstr(precalcfield);
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap = cellstr(magfield);
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [NaN NaN]; % We do not have short/long echo times !!?? [w.SHORT_ECHO_TIME w.LONG_ECHO_TIME]
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0; % 1= Magnitude Image is choosed to generate Mask Brain
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = 1; % Check if this is correct
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = w.READOUT_TIME; % readout time
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;  % non-EPI bases fieldmap
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;    % Jacobian use do not use
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Huttonish';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 15;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = cellstr(template);
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;

        % Loop for sessions
        for j=1:numel(w.sessions)
            % Get the first EPI file
            f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.imgPrefix  '.*\.nii$'], 1:1);
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(j).epi = cellstr(f);
        end
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1; % ===== NEW ====== % This will coregister the field map data to the selected EPI for each run/session
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
        %     matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = cellstr(anat); % For visualization, seems to crash if there is more than 1 session
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

        save(fullfile(w.subPath, 'SPM12_matlabbatch_1_Fieldmap.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
    end

    function DoSliceTiming(w,iS)
        clear matlabbatch;

        fprintf(' \n \n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' SliceTiming...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Loop for sessions
        matlabbatch{1}.spm.temporal.st.scans = {};
        for j=1:numel(w.sessions)
            % Get EPI raw files (filename starts with subject name)
            f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.imgPrefix  '.*\.nii$'], Inf);
            matlabbatch{1}.spm.temporal.st.scans{j} = cellstr(f);
        end
        matlabbatch{1}.spm.temporal.st.nslices = w.nSlices;
        matlabbatch{1}.spm.temporal.st.tr = w.TR;
        matlabbatch{1}.spm.temporal.st.ta = 0;
        matlabbatch{1}.spm.temporal.st.so = w.slice_times;
        matlabbatch{1}.spm.temporal.st.refslice = 0; % Using the first instead of the middle "w.refslice" ??
        matlabbatch{1}.spm.temporal.st.prefix = 'a';

        save(fullfile(w.subPath, 'SPM12_matlabbatch_2_SliceTiming.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

    function DoRealignUnwarp(w,iS)
        clear matlabbatch

        fprintf(' \n \n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' RealignUnwarp...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        EPI = {};
        % Loop for sessions
        for j=1:numel(w.sessions)
            expReg = ['^vdm5.*' num2str(j) '\.nii$'];
            vdm = cellstr(spm_select('FPList', w.fieldmapPath, expReg));  % vdm (voxel depplacement map) file or phase map

            % Get EPI sliced-time corrected files
            f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^a' w.imgPrefix  '.*\.nii$'], Inf);
            matlabbatch{1}.spm.spatial.realignunwarp.data(j).scans     = cellstr(f);
            matlabbatch{1}.spm.spatial.realignunwarp.data(j).pmscan    = vdm;
        end

        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0; % Register to 1st
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 4; % 4th degree Bspline
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0]; % Should I use "[0 1 0]", wrapping along the Y axis (because acquisition in PA direction) ???
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 6;    %  6th degree Bspline (slower)
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0]; % Should I use "[0 1 0]", wrapping along the Y axis (because acquisition in PA direction) ???
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

        save(fullfile(w.subPath, 'SPM12_matlabbatch_3_Realign&Unwarp.mat'),'matlabbatch');

        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

        DoPlotRealignFig(w,iS)

    end

    function DoRealign(w,iS)
        clear matlabbatch

        fprintf(' \n \n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' Realign...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Loop for sessions
        for j=1:numel(w.sessions)
            % Get EPI sliced-time corrected files
            f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^a' w.imgPrefix  '.*\.nii$'], Inf);
            matlabbatch{1}.spm.spatial.realign.estwrite.data{j} = cellstr(f);
        end

        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3; % The separation (in mm) between the points sampled in the reference image.  Smaller sampling distances gives more accurate results, but will be slower.
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % Register to 1st might be fine (and shorter) than register to mean (=1)...
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0]; % Should I use "[0 1 0]", wraping along the Y axis (because acquisition in PA direction) ???
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0]; % Should I use "[0 1 0]", wraping along the Y axis (because acquisition in PA direction) ???
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

        save(fullfile(w.subPath, 'SPM12_matlabbatch_3_Realign.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

        DoPlotRealignFig(w,iS)

    end

    function DoCoregister(w,iS)
        clear matlabbatch

        fprintf(' \n \n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' Coregister...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Get T1 structural file
        anat = spm_select('FPList', w.T1Path, ['^' w.imgPrefix '.*'  '.*\.nii$']);

        % Get mean realigned EPI
        meanRealign = spm_select('FPList', fullfile(w.funcPath, w.sessions{1}), ['^meanua' w.imgPrefix  '.*\.nii$']); % Mean is in session 1 folder

        clear matlabbatch
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(meanRealign);
        matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(anat);
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

        save(fullfile(w.subPath, 'SPM12_matlabbatch_4_Coregister.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

    function DoSegment(w,iS)
        clear matlabbatch

        fprintf(' \n \n');
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' Segment...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Select the (coregistered) T1 3D image
        coregAnat = cellstr(spm_select('FPList', w.T1Path, ['^' w.imgPrefix '.*'  '.*\.nii$']));

        % Get template of each cerebral tissue
        tmpGM          = {fullfile(spm('Dir'),'tpm', 'TPM.nii,1')};
        tmpWM          = {fullfile(spm('Dir'),'tpm', 'TPM.nii,2')};
        tmpCSF         = {fullfile(spm('Dir'),'tpm', 'TPM.nii,3')};
        tmpBone        = {fullfile(spm('Dir'),'tpm', 'TPM.nii,4')};
        tmpSoftTissue  = {fullfile(spm('Dir'),'tpm', 'TPM.nii,5')};
        tmpAirBck      = {fullfile(spm('Dir'),'tpm', 'TPM.nii,6')};

        clear matlabbatch;
        matlabbatch{1}.spm.spatial.preproc.channel.vols = coregAnat;
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; % Bias regularisation light
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;   % 60 mm cutoff
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];   % Save bias corrected
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = tmpGM;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];  %native
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 0];  %create the normalized version (wc*)
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = tmpWM;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];  %native
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 0];  %create the normalized version (wc*)
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = tmpCSF;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];  %native
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 0];  %create the normalized version (wc*)
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = tmpBone;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];  %native
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = tmpSoftTissue;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];  %native
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = tmpAirBck;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];   %  Forward + inverse deformation field

        save(fullfile(w.subPath, 'SPM12_matlabbatch_5_Segment.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

    function DoNormalise(w,iS)
        clear matlabbatch

        fprintf([' \n \n']);
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' Normalise...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Get Field Deformation image
        forwardDeformation = spm_select('FPList', w.T1Path, ['^y_' w.imgPrefix '.*\.nii$']);

        % Get bias corrected structural image (? difference with original anat ?)
        coregAnat = spm_select('FPList', w.T1Path, ['^m' w.imgPrefix '.*'  '.*\.nii$']);

        % Get Sliced EPI images of all runs
        EPI = {};
        % Loop on sessions
        for j=1:numel(w.sessions)
            % Get EPI Realigned files
            f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^ua' w.imgPrefix  '.*\.nii$'], Inf);
            EPI = vertcat(EPI, cellstr(f));
        end

        % Get c1  c2  and c3
        c1 = spm_select('FPList', w.T1Path, ['^c1' '.*\.nii$']);
        c2 = spm_select('FPList', w.T1Path, ['^c2' '.*\.nii$']);
        c3 = spm_select('FPList', w.T1Path, ['^c3' '.*\.nii$']);
        c1c2c3 = vertcat(c1, c2, c3);

        clear matlabbatch;
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {coregAnat};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};
        matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cellstr(EPI);
        matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
        matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [w.thickness w.thickness w.thickness];
        matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 6; % 6th degree B-Spline (slower)
        matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

        % The following step is not necessary because wc* files were created at the segmentation step ??
        % Or we need it to get the same voxel size as in the EPI to make the explicit mask ??
        %     matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};
        %     matlabbatch{3}.spm.spatial.normalise.write.subj.resample = cellstr(c1c2c3) ;
        %     matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
        %     matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [w.thickness w.thickness w.thickness];
        %     matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
        %     matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';

        save(fullfile(w.subPath, 'SPM12_matlabbatch_6_Normalize.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

    function DoSmooth(w,iS)
        clear matlabbatch;

        fprintf([' \n \n']);
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' Smooth...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        EPI = [];
        % Get Normalized EPI files of all sessions
        for j=1:numel(w.sessions)
            % Get EPI normalized files without dummy files
            f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^wua' w.imgPrefix  '.*\.nii$'], Inf);
            EPI = vertcat(EPI, cellstr(f));
        end
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(EPI);
        matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5]; % Smothing kernel FWHM
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';

        save(fullfile(w.subPath, 'SPM12_matlabbatch_7_Smooth.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

    function DoExplicitMask(w,iS)
        clear matlabbatch;

        fprintf([' \n \n']);
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' ExplicitMask...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Get normalized tissus (grey and white matter, CSF)
        wc1 = spm_select('FPList', w.T1Path, ['^wc1'  '.*\.nii$']);
        wc2 = spm_select('FPList', w.T1Path, ['^wc2'  '.*\.nii$']);
        wc3 = spm_select('FPList', w.T1Path, ['^wc3'  '.*\.nii$']);
        P = [wc1; wc2; wc3];

        matlabbatch{1}.spm.util.imcalc.input = cellstr(P);
        matlabbatch{1}.spm.util.imcalc.output = fullfile(w.T1Path, 'explicitMask_wc1wc2wc3_0.3.nii');
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 +i3)>0.3';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

        % Get non-normalized tissus (grey and white matter, CSF)
        c1 = spm_select('FPList', w.T1Path, ['^c1'  '.*\.nii$']);
        c2 = spm_select('FPList', w.T1Path, ['^c2'  '.*\.nii$']);
        c3 = spm_select('FPList', w.T1Path, ['^c3'  '.*\.nii$']);
        P = [c1; c2; c3];

        matlabbatch{2}.spm.util.imcalc.input = cellstr(P);
        matlabbatch{2}.spm.util.imcalc.output = fullfile(w.T1Path, 'explicitMask_c1c2c3_0.3.nii');
        matlabbatch{2}.spm.util.imcalc.outdir = {''};
        matlabbatch{2}.spm.util.imcalc.expression = '(i1 + i2 +i3)>0.3';
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = 1;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 4;

        matlabbatch{3}.spm.util.imcalc.input = cellstr(P);
        matlabbatch{3}.spm.util.imcalc.output = fullfile(w.T1Path, 'explicitMask_c1_0.01.nii');
        matlabbatch{3}.spm.util.imcalc.outdir = {''};
        matlabbatch{3}.spm.util.imcalc.expression = '(i1)>0.1';
        matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{3}.spm.util.imcalc.options.mask = 0;
        matlabbatch{3}.spm.util.imcalc.options.interp = 1;
        matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
        save(fullfile(w.subPath, 'SPM12_matlabbatch_8_Mask.mat'),'matlabbatch');
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

        % Coregister (& re-binarize) MNI-space masks with func to make sure it can be used
        % with "high_variance_confounds" python function, CosmoMVPA... (reslicing only
        % to correct voxel size didn't work...)
        clear matlabbatch;
        f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{1}), ['^swua' w.imgPrefix  '.*\.nii$'], 1);
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {f};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(fullfile(w.T1Path, 'explicitMask_wc1wc2wc3_0.3.nii'));
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        matlabbatch{2}.spm.util.imcalc.input =  cellstr(fullfile(w.T1Path, 'rexplicitMask_wc1wc2wc3_0.3.nii'));
        matlabbatch{2}.spm.util.imcalc.output = fullfile(w.T1Path, 'brexplicitMask_wc1wc2wc3_0.3.nii');
        matlabbatch{2}.spm.util.imcalc.outdir = {''};
        matlabbatch{2}.spm.util.imcalc.expression = 'i1>0.1';
        matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = 1;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

        % Coregister with non normalized non smoothed (for MVPA)
        clear matlabbatch;
        f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{1}), ['^ua' w.imgPrefix  '.*\.nii$'], 1);
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {f};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(fullfile(w.T1Path, 'explicitMask_c1_0.01.nii'));
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        matlabbatch{2}.spm.util.imcalc.input =  cellstr(fullfile(w.T1Path, 'rexplicitMask_c1_0.01.nii'));
        matlabbatch{2}.spm.util.imcalc.output = fullfile(w.T1Path, 'brexplicitMask_c1_0.01.nii');
        matlabbatch{2}.spm.util.imcalc.outdir = {''};
        matlabbatch{2}.spm.util.imcalc.expression = 'i1>0.1';
        matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = 1;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
    end

%% Other functions

    function DoPlotRealignFig(w,iS)

        fprintf([' \n \n']);
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' PlotRealignFig...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Compare with SPM figure to check if it is correct !!

        % ! Here we do one figure combining all runs (after realignement) !
        % Load data from rp_*.txt file
        data = []; amplitude = []; stdev = [];
        for j=1:numel(w.sessions)
            rpF = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^rp_' '.*\.txt$']);
            data = [data;load(deblank(rpF))];
            amplitude(j,:) = max(load(deblank(rpF))) - min(load(deblank(rpF))); % per run... if needed
            stdev(j,:)     = std(load(deblank(rpF)));                           % per run... if needed
        end

        % Figure
        printfig = figure;
        set(printfig, 'Name', ['Motion parameters: subject ' w.subjects{iS} ], 'Visible', 'on');

        subplot(2,1,1);
        plot(data(:,1:3)); legend('x translation', 'y translation','z translation'); grid on;
        title({['Motion parameters: translation (mm)'],...
            ['Average amplitude = ' num2str(mean(max(data(:,1:3))-min(data(:,1:3))), '%.02f' ) ', average stdev = ' num2str(mean(std(data(:,1:3))), '%.02f' )]},...
            'interpreter', 'none');
        subplot(2,1,2);
        plot(data(:,4:6)*180/pi); legend('pitch', 'roll','yaw'); grid on;
        title({['Motions parameters: rotation (dg)'],...
            ['Average amplitude = ' num2str(mean(max(data(:,4:6)*180/pi)-min(data(:,4:6)*180/pi)), '%.02f' ) ', average stdev = ' num2str(mean(std(data(:,4:6)*180/pi)), '%.02f' )]},...
            'interpreter', 'none');

        % Save and close figure
        mydate = date;
        filename = [w.niftidir filesep 'motion_' w.subjects{iS} '_' date '.png'];
        set(gcf, 'InvertHardcopy','off', 'PaperPosition',[0, 0, 22, 16],'PaperOrientation', 'portrait')
        print(printfig, '-dpng', '-r300', filename);
        close(printfig);
    end

    function DoSurfaceRenderAndBrainExtract(w,iS)
        clear matlabbatch;

        fprintf([' \n \n']);
        fprintf('-------------------------------------------------------------------------\n');
        fprintf(['                 ' w.subjects{iS} ' SurfaceRenderAndBrainExtract...\n']);
        fprintf('-------------------------------------------------------------------------\n');

        % Get normalized (bias corrected) structural image
        normAnat = spm_select('FPList', w.T1Path, ['^wm' '.*'  '.*\.nii$']);
        if size(normAnat,1)>1; normAnat = strrep(normAnat(1,:), ' ', ''); end % if there is more than 1 image (from a previous preprocessing) keep the first (and remove space characters)
        [a1,b1] = fileparts(normAnat);

        % Get normalized tissues (grey and white matter, CSF)
        wc1 = spm_select('FPList', w.T1Path, ['^wc1'  '.*\.nii$']);
        wc2 = spm_select('FPList', w.T1Path, ['^wc2'  '.*\.nii$']);
        wc3 = spm_select('FPList', w.T1Path, ['^wc3'  '.*\.nii$']);

        % Get non-normalized structural image
        Anat = spm_select('FPList', w.T1Path, ['^' w.imgPrefix '.*'  '.*\.nii$']);
        [a2,b2] = fileparts(Anat);

        % Get non-normalized tissues (grey and white matter, CSF)
        c1 = spm_select('FPList', w.T1Path, ['^c1'  '.*\.nii$']);
        c2 = spm_select('FPList', w.T1Path, ['^c2'  '.*\.nii$']);
        c3 = spm_select('FPList', w.T1Path, ['^c3'  '.*\.nii$']);

        % Get Field Deformation image
        forwardDeformation = spm_select('FPList', w.T1Path, ['^y_' w.imgPrefix '.*\.nii$']);


        %======% Normalized brain extraction with threshold and render %======%
        thresh = 0.6; % Skull stripping threshold
        P = {normAnat; wc1; wc2};
        exp = ['i1 .*((i2+i3)>' num2str(thresh) ')'];

        matlabbatch{1}.spm.util.imcalc.input = (P);
        matlabbatch{1}.spm.util.imcalc.output = fullfile(w.T1Path, [b1 '_brain_extracted.nii']);
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = exp;
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

        matlabbatch{2}.spm.util.render.extract.data(1) =  {fullfile(w.T1Path, [b1 '_brain_extracted.nii'])};
        matlabbatch{2}.spm.util.render.extract.mode = 3;
        matlabbatch{2}.spm.util.render.extract.thresh = 0.5;

        %======%  Extraction of non-normalized brain, including CSF... %======%
        P = {c1; c2; c3; Anat};
        exp =  '(i1 + i2 + i3) .* i4';

        matlabbatch{3}.spm.util.imcalc.input = (P);
        matlabbatch{3}.spm.util.imcalc.output = fullfile(w.T1Path, ['Brain.nii']);
        matlabbatch{3}.spm.util.imcalc.outdir = {''};
        matlabbatch{3}.spm.util.imcalc.expression = exp;
        matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{3}.spm.util.imcalc.options.mask = 0;
        matlabbatch{3}.spm.util.imcalc.options.interp = 1;
        matlabbatch{3}.spm.util.imcalc.options.dtype = 4;

        %======%  Normalization of extracted brain (including CSF) %======%
        matlabbatch{4}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
        matlabbatch{4}.spm.spatial.normalise.write.subj.resample = {fullfile(w.T1Path, ['Brain.nii'])};
        matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
        matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
        matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w';

        % Run batch
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);

    end

cd(cdd)
end