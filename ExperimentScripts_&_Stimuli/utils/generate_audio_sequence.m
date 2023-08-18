function [ output_left, output_right ] = generate_audio_sequence( sequence,sound_duration,audiofreq,SOA,path_to_save,trigger_in_right_ear )
%GENERATE_AUDIO_SEQUENCE This function creates the sound sequence
% corresponding to a given input sequence
% 'sequence' should be coded as a sequence of [0 1 0 etc]
% sound duration is the duration of the sound. default 0.070
% trigger_in_right_ear: set this to true if you want to have some triggers
% in the left ear. Code 1 for 'A' and 1 for 'B'.

% dv = PsychPortAudio('GetDevices');
trigger_duration = 0.010;

sound_A = create_sound([294, 294*2, 294*4], 0.007, sound_duration, audiofreq);
sound_B = create_sound([440 440*2 440*4],0.007,sound_duration, audiofreq);
trigger_A = 1*ones(1,audiofreq*trigger_duration);
trigger_B = 1*ones(1,audiofreq*trigger_duration);

output_left = [];
output_right = [];
trigger = [];

len = length(sequence);
ISI = SOA-sound_duration;

silence = zeros(1,audiofreq*ISI);
silence_trigger = zeros(1,audiofreq*(SOA-trigger_duration));

for i =1:len
    if sequence(i)==0
        trigger = [trigger,trigger_A,silence_trigger];
        output_left = [output_left,sound_A,silence];
        output_right = [output_right,sound_A,silence];
    elseif sequence(i)==1
        trigger = [trigger,trigger_B,silence_trigger];
        output_left = [output_left,sound_B,silence];
        output_right = [output_right,sound_B,silence];
    end
end
name = num2str(sequence);
name= name(find(~isspace(name)));

if trigger_in_right_ear
    output_right = trigger;
end

if ~strcmp(path_to_save,'false')
    save_name = [path_to_save,name,'.wav'];
    audiowrite(save_name,[output_left',output_right'],audiofreq)
end


end

