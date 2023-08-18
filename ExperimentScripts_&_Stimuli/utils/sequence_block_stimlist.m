function [ list_seq_and_viol,nonviol_viol ] = sequence_block_stimlist( sequence,n_repetitions,n_violations,violation_po,n_habituation,n_max_consecutive_violations )
%SEQUENCE_BLOCK_STIMLIST Generates the list of sequences+their
%violations that will be presented to the participants. nonviol_viol
%contains the list of violation positions. If 0, then the sequence is not
%violated.


list_stims_original_sequences = repmat(sequence,n_repetitions,1);
n_violations_per_position = n_violations/length(violation_po);

viol_seq = [];
viol_position = [];
for k = violation_po
    vio = sequence;
    vio(k) = 1 - sequence(k);
    viol_seq = [viol_seq; vio];
    viol_position = [viol_position; k];
end

list_stims_violated_sequences = repmat(viol_seq,n_violations_per_position,1);

% list of the stimuli that will come after habituation

% before randomization
list_seq_and_viol = vertcat(list_stims_original_sequences,list_stims_violated_sequences);
nonviol_viol = [zeros(n_repetitions,1);repmat(viol_position,n_violations_per_position,1)];

% let's randomize controlling for the number of consecutive violations
new_indices = permutation_balanced(n_repetitions,n_violations,n_max_consecutive_violations);
list_seq_and_viol = list_seq_and_viol(new_indices,:);
nonviol_viol = nonviol_viol(new_indices');

% add at the beginning n_habituation original sequences 
list_seq_and_viol = vertcat(repmat(sequence,n_habituation,1),list_seq_and_viol);
nonviol_viol = vertcat(zeros(n_habituation,1),nonviol_viol);

% -----------------------------------------------------------------------------------------------------------------------
    function shuffed_and_balanced_indices = permutation_balanced(n_repetitions,n_violations,n_max_consecutive_violations)
        % This function shuffles the indices such that there are maximum
        % n_max_consecutive_violations consecutive violations
        
        A = ones(1,n_repetitions);
        B = 2*ones(1,n_violations);
        AB = [A,B];
        is_ok = false;
        while ~is_ok
            shuffed_and_balanced_indices = randperm(n_repetitions+n_violations);
            shuff_AB = AB(shuffed_and_balanced_indices);
            ok =0;
            for k = 1:(n_repetitions+n_violations-n_max_consecutive_violations)
                if sum(shuff_AB(k:k+n_max_consecutive_violations))==(n_max_consecutive_violations+1)*2
                    ok = 1;
                end
            end   
            if ok==0
                is_ok = true ;
            end
        end
             
        
    end



end

