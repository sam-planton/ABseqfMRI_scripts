%% Process AP & PA spin-echo images with FSL TopUp (to be used with SPM Fieldmap)
function Do_ABseq16fMRI_1_ProcessTopUp_Fieldmap(w,iS)

        
    % Script parts from: F. Meyniel "fmspm12batch_AddTopupCorrection_job1sub.m"
    % Also see https://lcni.uoregon.edu/kb-articles/kb-0003

    % unix(cmd) % DID NOT INITIALLY WORK IN MATLAB (ONLY IF COPY-PASTED IN A TERMINAL):
    % "/neurospin/local/matlab/R2017a/sys/os/glnxa64/libstdc++.so.6: version`GLIBCXX_3.4.21' not found (required by /usr/lib/fsl/5.0/libnewimage.so)"
    % Using the following command before launching Matlab, forcing Matlab to use system libstdc++, seems to solve the issue:
    % alias matlab='LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /neurospin/local/bin/matlab -desktop'

    fprintf('=======================================================================\n');
    fprintf([w.subjects{iS} ' Fieldmap fsl TopUp...  (can take ~20min) \n']);       
    fprintf('=======================================================================\n');

    load(fullfile(w.niftidir, w.subjects{iS}, 'func', 'SliceTimingInfo.mat'))

    fdir = fullfile(w.niftidir, w.subjects{iS}, 'fieldmap');
    fdir = string(fdir);

    % Identify the AP PA file for calibrating the deformation, with the nii extention ,1
    B0_AP = spm_select('ExtFPList', fdir, '.*AP.*\.nii', 1);
    B0_PA = spm_select('ExtFPList', fdir, '.*PA.*\.nii', 1);

    % check that files exist
    if strcmp(B0_AP, ''); error('cannot find B0_AP file'); end
    if strcmp(B0_PA, ''); error('cannot find B0_PA file'); end

    % ESTIMATE THE DEFORMATION MAP WITH FSL & TOPUP
    % =========================================================================
    cmd = sprintf('cd %s; fslmerge -t b0_APPA %s %s', ...
        fdir, B0_AP(1:end-2), B0_PA(1:end-2));
    system(cmd)
    % clipboard('copy',cmd)

    % Create a text file with the direction of phase encoding.
    % A>P is -1; P>A is 1
    % (could be checked with Romain's script topup_param_from_dicom).
    cmd = sprintf('cd %s; echo $''0 -1 0 %6.5f \n0 1 0 %6.5f'' > acq_param.txt', ...
        fdir, total_readout_time_fsl, total_readout_time_fsl);
    system(cmd)

    % Compute deformation with Topup (this takes ~20 minutes and only 1 CPU)
    % The result that will be used is APPA_DefMap
    % For sanity checks, sanitycheck_DefMap is the deformation field and
    % sanitycheck_unwarped_B0 are the corrected images
    % NB: in the b02b0.cnf, subsamp starts =2, which implies that the size
    % of the image cannot be an odd number. Here we force subsamplin1845rd?ceaACDg =1.
    fprintf('\n Computing the APPA deformation with Topup... (~20min)')
    cmd = sprintf(['cd %s; topup ', ...
        '--imain=b0_APPA --datain=acq_param.txt --config=b02b0.cnf ', ...
        '--out=APPA_DefMap --fout=%s_topup_fieldmap --iout=sanitycheck_unwarped_B0 ', ...
        '--subsamp=1,1,1,1,1,1,1,1,1'], ...
        fdir, w.subjects{iS});
    system(cmd)

    %============ TOPUP INFO ============%
    % --out
    % The value for the --out parameter is used as the basename for the output files 
    % produced by topup. Given a set of files, let us use the example in --imain above, 
    % topup will estimate a field, that is common to all scans, and movement parameters 
    % for scans 2-n. The movement parameters will encode the positions of scans 2-n 
    % relative to scan 1'
    % ....

    % --fout
    % Specifies the name of an image file containing the estimated field in Hz. The actual
    % information is the same as that in the *_fieldcoef.nii.gz file from the --out parameter.
    % They are different in that the --fout image contains the actual voxel-values 
    % (as opposed to spline coefficients) which makes it easier to use with applications
    % that cannot read the topup output format. It is for example useful if you want to 
    % feed it as a fieldmap into FEAT. 

    % --iout
    % Specifies the name of a 4D image file that contains unwarped and movement corrected
    % images. Each volume in the --imain will have a corresponding corrected volume in --iout.
    % This output uses traditional interpolation and Jacobian modulation (see applytopup) 
    % and is used mainly as a sanity check that things have worked and that topup has estimated
    % a reasonable field. It can also be useful for making an undistorted mask for use in 
    % eddy by running BET on the first volume (or the average of all volumes) of the output. 


    % IN https://lcni.uoregon.edu/kb-articles/kb-0003:
    %  --fout=my_fieldmap --iout=se_epi_unwarped
    % Finally, if you will be using this with Feat, you'll need a single magnitude image, 
    %     and a brain extracted version of that image. You can get that this way:
    cmd = sprintf(['cd %s; fslmaths ', ...
                    'sanitycheck_unwarped_B0 -Tmean %s_topup_magnitude'],fdir, w.subjects{iS});
    system(cmd)

    % Convert from .nii.gz to .nii
    disp('DECOMPRESS NII.GZ USING gzip')
    cmd = sprintf('gzip -d %s', fullfile(fdir, [w.subjects{iS} '_topup_magnitude.nii.gz'])); system(cmd);
    cmd = sprintf('gzip -d %s', fullfile(fdir, [w.subjects{iS} '_topup_fieldmap.nii.gz'])); system(cmd);

end