function generate_sequences_sound_files(list_seq,sound_duration,audiofreq, SOA,path_to_save,trigger)

    % set 'trigger' to true if you want to have trigger pulses in the left ear
    
%     dv = PsychPortAudio('GetDevices');
%     audfreq = dv(1).DefaultSampleRate;

    out_l = [];
    out_r = [];
    for n = 1:size(list_seq,1)
        sequence = list_seq(n,:);
        [ output_left, output_right ]  = generate_audio_sequence(sequence,sound_duration, audiofreq,SOA,'false',trigger);
        fname = fullfile(path_to_save, [strrep(num2str(sequence),' ','') '.wav']);
        audiowrite(fname,[output_left',output_right'],audiofreq)
        disp(fname)
    end
   
end

