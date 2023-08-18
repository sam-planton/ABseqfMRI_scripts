%% ======= GETTING READY ======= %%
close all; clearvars; clc;
 
addpath('utils')

subjectID=inputdlg('Subject ID?','Output',1,{'sub_00'}); subjectID = subjectID{1};

% Directories & output file
rootdir = pwd;
% idcs   = strfind(rootdir,filesep);
% outdir = fullfile(rootdir(1:idcs(end)-1), 'stimuli_presentation_lists', subjectID);
outdir = fullfile(rootdir, 'stimuli_presentation_lists', subjectID);
Filename = [outdir filesep subjectID '_presentation_lists.mat'];
if ~exist(outdir, 'dir'); mkdir(outdir); end
if exist(Filename)
    replace=questdlg(['File ' Filename, ' already exists. Do you want to replace it?']);
    if strcmp( replace, 'No' )
        error('Exit')
    end
end

% A run is constituted of n_habituation repetitions of the original
% sequence without any violation.
% Then there are n_original of the original sequence and n_violation
% violated sequences presented in a shuffled way.
n_habituation = 10;
n_original    = 12;
n_violations  = 24;
% this is the maximal number of consecutive violated sequences
n_max_consecutive_violations = 3;

%% ======= DEFINE THE SEQUENCES AND THE POSITIONS OF THE VIOLATIONS  ======= %%

% /!\ Check the script: utils\sequences_and_violations_list.m /!\ 

% ======= Get the list of sequences and violations ===== %
sequences_and_violations = sequences_and_violations_list();
nsemiruns = size(sequences_and_violations,1);
% ===================================================== %

%% ======= GENERATE THE DIFFERENT BLOCKS ======= %%
% saves a 'presentation_list_runXX.mat' file for each run in the subject
% folder (and a run_contents.csv that summarizes the list of blocks)
% 8 runs, each containing 2 semiruns (2 lists)

% Block indices: Hab block (10 trials), then 3 test blocks (of 12 trials)
blockID_list = [ones(n_habituation,1); ones((n_original+n_violations)/3,1)*2; ones((n_original+n_violations)/3,1)*3; ones((n_original+n_violations)/3,1)*4];

sequence_perm = generate_balanced_set(nsemiruns);
all_semiruns = 1:nsemiruns;
table_order = table(all_semiruns',sequences_and_violations(sequence_perm,1),'Variablenames',{'run_number','Stim'});
writetable(table_order,fullfile(outdir, 'run_contents.csv'));

semirun_number = 1; run_number =1;
while semirun_number<nsemiruns
    
    % FIRST SEMIRUN
    seq_num = sequence_perm(semirun_number);
    sequence = sequences_and_violations{seq_num,1}; % the considered sequence for block i
    violation_pos = sequences_and_violations{seq_num,2}; % the positions of the violations
    % the following function generates the list of the sequences that will be presented and the position of the violations. If it is 0, then it means no violation
    % We exclude 3 standard trials that will be re-added later (1st of eachtest block)
    [list_seq, list_violation_position] = sequence_block_stimlist(sequence,n_original-3,n_violations,violation_pos,n_habituation,n_max_consecutive_violations);
    % Re-add one standard at the beginning of each block
    tmpN = (n_original+n_violations)/3-1;
    newlist_violation_position = [list_violation_position(1:n_habituation); 0; list_violation_position(n_habituation+1:n_habituation+tmpN);0; list_violation_position(n_habituation+tmpN+1:n_habituation+tmpN+tmpN);0; list_violation_position(n_habituation+tmpN+tmpN+1:n_habituation+tmpN+tmpN+tmpN)];    
    newlist_seq = [list_seq(1:n_habituation,:); list_seq(1,:); list_seq(n_habituation+1:n_habituation+tmpN,:);list_seq(1,:); list_seq(n_habituation+tmpN+1:n_habituation+tmpN+tmpN,:);list_seq(1,:); list_seq(n_habituation+tmpN+tmpN+1:n_habituation+tmpN+tmpN+tmpN,:)];    
    list_seq_for_table = arrayfun(@(x) mat2str(x),newlist_seq);
    semirun1_presentation_list = table(list_seq_for_table,newlist_violation_position,'Variablenames',{'Presented_sequence','Position_Violation'});
    condition = repmat(string(get_seq_name(sequence)),size(list_seq_for_table,1),1);
    trialnum = (1:size(list_seq_for_table,1))'; semirun = repmat(semirun_number,1,size(list_seq_for_table,1))';
    block = blockID_list; run = repmat(run_number, size(list_seq_for_table,1),1); deviant = newlist_violation_position>0;
    semirun1_presentation_list = [semirun1_presentation_list table(condition, run, semirun, block, trialnum, deviant)];
    semirun_number = semirun_number + 1;
    
    % SECOND SEMIRUN
    seq_num = sequence_perm(semirun_number);
    sequence = sequences_and_violations{seq_num,1}; % the considered sequence for block i
    violation_pos = sequences_and_violations{seq_num,2}; % the positions of the violations
    % the following function generates the list of the sequences that will be presented and the position of the violations. If it is 0, then it means no violation
    % We exclude 3 standard trials that will be re-added later (1st of eachtest block)
    [list_seq, list_violation_position] = sequence_block_stimlist(sequence,n_original-3,n_violations,violation_pos,n_habituation,n_max_consecutive_violations);
    % Re-add one standard at the beginning of each block
    tmpN = (n_original+n_violations)/3-1;
    newlist_violation_position = [list_violation_position(1:n_habituation); 0; list_violation_position(n_habituation+1:n_habituation+tmpN);0; list_violation_position(n_habituation+tmpN+1:n_habituation+tmpN+tmpN);0; list_violation_position(n_habituation+tmpN+tmpN+1:n_habituation+tmpN+tmpN+tmpN)];    
    newlist_seq = [list_seq(1:n_habituation,:); sequence; list_seq(n_habituation+1:n_habituation+tmpN,:);sequence; list_seq(n_habituation+tmpN+1:n_habituation+tmpN+tmpN,:);sequence; list_seq(n_habituation+tmpN+tmpN+1:n_habituation+tmpN+tmpN+tmpN,:)];    
    list_seq_for_table = arrayfun(@(x) mat2str(x),newlist_seq);
    semirun2_presentation_list = table(list_seq_for_table,newlist_violation_position,'Variablenames',{'Presented_sequence','Position_Violation'});
    condition = repmat(string(get_seq_name(sequence)),size(list_seq_for_table,1),1);
    trialnum = (1:size(list_seq_for_table,1))'; semirun = repmat(semirun_number,1,size(list_seq_for_table,1))';
    block = blockID_list; run = repmat(run_number, size(list_seq_for_table,1),1); deviant = newlist_violation_position>0;
    semirun2_presentation_list = [semirun2_presentation_list table(condition, run, semirun, block, trialnum, deviant)];
    semirun_number = semirun_number + 1;
    
    presentation_list = [semirun1_presentation_list; semirun2_presentation_list];
    fname = fullfile(outdir, ['run' num2str(run_number, '%02.f') '_presentation_list.mat']);
    save(fname, 'presentation_list')
    disp(fname)
    
    run_number = run_number + 1;
    
end

function sequence_name = get_seq_name(bin_sequence)
    strseq = strrep(num2str(bin_sequence),' ','');
    if     strcmp(strseq, '1111111111111111') || strcmp(strseq, '0000000000000000')
        sequence_name = 'Repeat';
    elseif strcmp(strseq, '1010101010101010') || strcmp(strseq, '0101010101010101')
        sequence_name = 'Alternate';
    elseif strcmp(strseq, '1100110011001100') || strcmp(strseq, '0011001100110011')
        sequence_name = 'Pairs';
    elseif strcmp(strseq, '1111000011110000') || strcmp(strseq, '0000111100001111')
        sequence_name = 'Quadruplets';
    elseif strcmp(strseq, '1100101011001010') || strcmp(strseq, '0011010100110101')
        sequence_name = 'Pairs+Alt';
    elseif strcmp(strseq, '1111000011001010') || strcmp(strseq, '0000111100110101')
        sequence_name = 'Shrinking';
    elseif strcmp(strseq, '1011100001001110') || strcmp(strseq, '0100011110110001')
        sequence_name = 'Complex';
    end
end