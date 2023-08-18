function generate_run_sequences(list_seq_and_viol,sound_duration,ISI, SOA,path_to_save,subject_number,run_number,trigger)

% set 'trigger' to true if you want to have trigger pulses in the right ear
dv = PsychPortAudio('GetDevices');
audfreq = dv(1).DefaultSampleRate;
silence = zeros(1,audfreq*ISI);

out_l = [];
out_r = [];
for n = 1:size(list_seq_and_viol,1)
    sequence = list_seq_and_viol(n,:);
    [ output_left, output_right ]  = generate_audio_sequence(sequence,sound_duration, SOA,'false',trigger);
    out_l = [out_l,output_left];
    out_l = [out_l,silence];
    out_r = [out_r,output_right];
    out_r = [out_r,silence];

end

audiowrite([path_to_save,subject_number,'/',run_number,'.wav'],[out_l',out_r'],audfreq)
end

