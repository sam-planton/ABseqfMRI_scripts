close all; clearvars; clc;

addpath('utils')

% =========== Get the list of sequences and violations ========= %
% /!\ Check the script: utils\sequences_and_violations_list.m /!\ 
sequences_and_violations = sequences_and_violations_list();
% ============================================================== %

% Global parameters 
sound_duration = 0.050;% in seconds
SOA            = 0.250; % in seconds
audiofreq      = 48000; % Hz

% Where you save the stimuli
rootdir = pwd;
path_to_save = fullfile(rootdir, 'stimuli');
if ~exist(path_to_save, 'dir'); mkdir(path_to_save); end

% should we send a trigger in the left ear ?
trigger = false;   
    
for nseq = 1:size(sequences_and_violations,1)
    sequence = sequences_and_violations{nseq,1}; % the considered sequence for block i
    violation_pos = sequences_and_violations{nseq,2}; % the positions of the violations

    viol_seq = [];
    viol_position = [];
    for k = violation_pos
        vio = sequence;
        vio(k) = 1 - sequence(k);
        viol_seq = [viol_seq; vio];
        viol_position = [viol_position; k];
    end

    list_stims_sequences = [sequence; viol_seq]; % original sequence + violated version

    % Generate the sound files
    InitializePsychSound
    generate_sequences_sound_files(list_stims_sequences, sound_duration, audiofreq, SOA, path_to_save, trigger)

end