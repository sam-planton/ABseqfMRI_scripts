% global parameters 
sound_duration = 0.05;% in seconds
SOA = 0.250; % in seconds
ISI = 0.500; % in seconds
% A run is constituted of n_habituation repetitions of the original
% sequence without any violation.
% Then there are n_original of the original sequence and n_violation
% violated sequences presented in a shuffled way.

n_habituation = 10;
n_original = 12;
n_violations = 24;
% this is the maximal number of consecutive violated sequences
n_max_consecutive_violations = 4;

% Where you save the stimuli
path_to_save = '/Volumes/ABseq_stim/stims_14102019/';

% should we send a trigger in the left ear ?
trigger = true;

% Initialize PsychPortAudio
InitializePsychSound;

%% DEFINE THE SEQUENCES AND THE POSITIONS OF THE VIOLATIONS
sequences_and_violations = cell(14,2);
% 2 versions for each sequence
s1 = 'AAAAAAAAAAAAAAAA'-'A';
s2 = 'ABABABABABABABAB'-'A';
s3 = 'AABBAABBAABBAABB'-'A';
s4 = 'AAAABBBBAAAABBBB'-'A';
s5 = 'AABBABABAABBABAB'-'A';
s6 = 'AAAABBBBAABBABAB'-'A';
s7 = 'ABAAABBBBABBAAAB'-'A';
s8 = abs('AAAAAAAAAAAAAAAA'-'B');
s9 = abs('ABABABABABABABAB'-'B');
s10 = abs('AABBAABBAABBAABB'-'B');
s11= abs('AAAABBBBAAAABBBB'-'B');
s12 = abs('AABBABABAABBABAB'-'B');
s13 = abs('AAAABBBBAABBABAB'-'B');
s14 = abs('ABAAABBBBABBAAAB'-'B');

violation_positions = {[9, 12, 13, 15];
                       [9, 12, 14,15];
                       [10, 11, 14, 15];
                       [9, 12, 13, 15];
                       [10, 11, 14, 15];
                       [10, 11, 14, 15];
                       [9, 12, 14, 15];
                       [9, 12, 13, 15];
                       [9, 12, 14,15];
                       [10, 11, 14, 15];
                       [9, 12, 13, 15];
                       [10, 11, 14, 15];
                       [10, 11, 14, 15];
                       [9, 12, 14, 15]};
                   
sequences_and_violations(:,1) = {s1;s2;s3;s4;s5;s6;s7;s8;s9;s10;s11;s12;s13;s14};
sequences_and_violations(:,2) = violation_positions;

%% SESSION PARAMETERS


for session_number = 13:19

    sequence_perm = generate_balanced_set();
    mkdir([path_to_save,sprintf('%i',session_number)]);
    all_runs = 1:14;
    table_order = table(all_runs',sequences_and_violations(sequence_perm,1),'Variablenames',{'run_number','Stim'});
    writetable(table_order,[path_to_save,sprintf('%i',session_number),'/run_contents.csv']);


    %% GENERATE THE DIFFERENT BLOCKS
    for i = 1:14
        seq_num = sequence_perm(i);
        sequence = sequences_and_violations{seq_num,1}; % the considered sequence for block i
        violation_pos = sequences_and_violations{seq_num,2}; % the positions of the violations
        % the following function generates the list of the sequences that will
        % be presented and the position of the violations. If it is 0, then it
        % means no violation
        [list_seq, list_violation_position] = sequence_block_stimlist(sequence,n_original,n_violations,violation_pos,n_habituation,n_max_consecutive_violations);
        % save the presented sequence, the position of the violation
        list_seq_for_table = arrayfun(@(x) mat2str(x),list_seq);
        Table_to_save = table(list_seq_for_table,list_violation_position,'Variablenames',{'Presented_sequence','Position_Violation'});
        writetable(Table_to_save,[path_to_save,sprintf('%i',session_number),'/',sprintf('info_run%i',i),'.csv']);
        % GENERATE THE SOUNDTRACK
        generate_run_sequences(list_seq,sound_duration,ISI, SOA,path_to_save,sprintf('%i',session_number),num2str(i),trigger)
    end

end



