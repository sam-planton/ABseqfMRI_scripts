
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% ABSeq **************************************************** %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ======= OVERVIEW OF THE EXPERIMENT ======= %%
%
% 
%
%
%
%

%% ============= GETTING READY ============== %%
clear all;                   % Clean Workspace
clc;                         % Clean command window
AssertOpenGL;                % Make sure the script is running on Psychtoolbox-3
varscreen = 1;               % "1" for SkipSyncTests
addpath('utils')             % Folder containing functions for IO32 markers // TTL
no_screen = false;           % Set to true if fixation cross is not needed

% Root directory
rootdir = pwd;

% To continue a session...       
if exist('run_in_progress.mat')
    load('run_in_progress.mat')
    answ=questdlg(['Continue for subject ' run_in_progress{1} ', run ' num2str(run_in_progress{2}+1) '?']);
%     answ=input(['Continue for subject ' run_in_progress{1} ', run ' num2str(run_in_progress{2}+1) '?']);
    if strcmp(answ, 'Yes')
        subjectID = run_in_progress{1};
        run_number= run_in_progress{2}+1;
    else
        subjectID=inputdlg('Subject ID?','Output',1,{'sub_00'}); subjectID = subjectID{1};
        run_number=inputdlg('Run number?'); run_number = str2num(run_number{1});
    end
else
    % Subject ID
    subjectID=inputdlg('Subject ID?','Output',1,{'sub_00'}); subjectID = subjectID{1};
    % Run number
    run_number=inputdlg('Run number?'); run_number = str2num(run_number{1});
end      
 
% Choose sound device
clear PsychPortAudio;
audiodevices=PsychPortAudio('GetDevices');
audiodevicenames={audiodevices.DeviceName};
[logicalaudio,locationaudio]=ismember({'Aureon5.1MkII: USB Audio (hw:1,0)','Aureon5.1MkII: USB Audio (hw:2,0)'},audiodevicenames);% digital
% [logicalaudio,locationaudio]=ismember({'Aureon5.1MkII: USB Audio (hw:2,0)'},audiodevicenames);% analog
audiodeviceindex=audiodevices(locationaudio(logicalaudio)).DeviceIndex;

InitializePsychSound(audiodeviceindex);
sound_device = audiodeviceindex;
audiofreq = audiodevices(sound_device+1).DefaultSampleRate; % Hertz

% Directory for the log file
idcs   = strfind(rootdir,filesep);
outdir = fullfile(rootdir(1:idcs(end)-1), 'Logs', subjectID);
if ~exist(outdir, 'dir')
   mkdir(outdir)
end

% fMRI acquisition
TR      = 1.81;              % TR, in seconds, used for the fixation cross duration (after the TTL)
nb_TR   = 272;               % Number of TRs for the whole run
run_duration = TR*nb_TR;     % Run duration (s)

%% ============= LOAD THE STIMULI/TRIAL LIST ============== %%

% Delay between sound files (s)t
ITI = 0.500;

% Rest block duration (s)
rest_block_dur = 8.0;
 
% Load presentation list
presentation_list_file = fullfile(rootdir, 'stimuli_presentation_lists', subjectID, ['run' num2str(run_number, '%02.f') '_presentation_list.mat']);
load(presentation_list_file)

% Prepare list of sound files and load sound files
stim_list = string(presentation_list.Presented_sequence);
[seqfiles,ia,stimIDlist] = unique(stim_list);
audio_data = {};
for ii = 1 : numel(seqfiles)
    waveData = audioread(fullfile(rootdir, 'stimuli', [char(seqfiles(ii)) '.wav']));
    audio_data{ii} = waveData';
end

try 
%% ============= Experiment Environment (Screen) =============== %%

% Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
% they are loaded and ready when we need them - without delays:
KbCheck;
GetSecs;

if no_screen == false
    
    % /!\ WITH UBUNTU: running the helper scripts XOrgConfCreator and
    % XOrgConfSelector to automatically create and install optimized xorg.conf configuration files for the X-Server /!\ 
    
    % Get screenNumber of stimulation display
    screenNumber = 1; %max(Screen('Screens'));
    framerate = Screen(screenNumber,'FrameRate');
    fprintf('\nscreenNumber= %d and framerate=%d\n', screenNumber, framerate);
    Screen('Preference', 'SkipSyncTests', varscreen); % 1=skip! Use the value 1 with caution
    Screen('Preference', 'VisualDebugLevel', 0); % No exclamation point ??
    
    % Old GDI text renderer to avoid issues ? // not with Ubuntu
%     Screen('Preference', 'TextRenderer', 0);

    HideCursor; % hide the mouse cursor
    black=0;   white=255;   grey=80;% (white+black)/2; % screen colors
    inc=abs(white-grey); % Contrast increment range for given white and gray values

    % Open a double buffered fullscreen window on the stimulation screen 'screenNumber'. 
    % 'w' is the handle used to direct all drawing commands to that window.
    % 'screenRect' is a rectangle defining the size of the window.
    % See "help PsychRects" for help on such rectangles and useful helper functions:
    [w,screenRect] = Screen('OpenWindow', screenNumber,grey,[],[],2);
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel); % set high priority for max performance
    ifi = Screen('GetFlipInterval', w); % Query the frame duration

    [screenXpixels, screenYpixels] = Screen('WindowSize', w); % Get the size of the on screen window
    [xCenter, yCenter] = RectCenter(screenRect);              % Get the centre coordinate of the window

    % Generate Fixation cross
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % Set up alpha-blending for smooth (anti-aliased) lines
    fixCrossDimPix = 16; % Here we set the size of the arms of our fixation cross
    [xCenter, yCenter] = RectCenter(screenRect);
    xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    allCoords = [xCoords; yCoords];
    lineWidthPix = 3; % Set the line width for our fixation cross  
end

%% ========== KEYBOARD ... ========== %%
KbName('UnifyKeyNames'); 
escapeKey = KbName('ESCAPE');
spaceKey = KbName('space');

%% ============================ %
%%======= RUN EXPERIMENT ======%%
%%=============================%%
fsep = strfind(presentation_list_file,filesep);
disp(' '); disp(['#===== Running: ' presentation_list_file(fsep(end-1):end) ' =====#']);

if no_screen == false
    % Introduction message
    Screen('TextFont', w, 'Segoe UI');
    Screen('TextSize', w, 32); % Setup the text type for the window
    DrawFormattedText(w, 'Préparez-vous, l''expérience va débuter.\n\n\n\n Veillez à rester concentré(e), à conserver le regard au centre de l''écran\n\n et à ne pas bouger tout au long de l''expérience' , 'center', 'center', white);
    Screen('Flip', w);
    WaitSecs(2);
end

% Play a beep to initialize and avoid delay/distorsion on the first sound ?
myBeep = MakeBeep(300, 0.050, audiofreq); % 300 Hz, 50 ms
pahandle = PsychPortAudio('Open', sound_device, 1, 1, audiofreq, 2); % http://psychtoolbox.org/docs/PsychPortAudio-Open
PsychPortAudio('FillBuffer', pahandle, [myBeep; myBeep]); 
PsychPortAudio('Start', pahandle, 1, 0, 1); % http://psychtoolbox.org/docs/PsychPortAudio-Start

% commandwindow;
sig=waitTTL; % Calling the waitTTL script
TTL_time = GetSecs;
if sig == 1

    tstart=GetSecs;
    fprintf('\nRun_start_time = %.3f\n',tstart);

    % send event marker to Eyelink system
%         Eyelink('Message', 'StartTriggerReceived');

    %% ========== RUN TRIAL LOOP ========== %%
    waitForStart = 1; % Should we wait for the device to really start (1 = yes) http://psychtoolbox.org/docs/PsychPortAudio-Start
    startCue = 0;

    % Display fixation cross
    if no_screen == false
        Screen('DrawLines', w, allCoords, lineWidthPix, white, [xCenter yCenter], 2); % Fixation cross
        Screen('Flip', w);
    end

    % First rest block
    fprintf('\n   << %d secs rest block >>  ', rest_block_dur)
    WaitSecs(rest_block_dur);

    % All trials loop
    ttime = GetSecs;
    for trial_num = 1 : numel(stimIDlist)

        stimID = stimIDlist(trial_num);
        current_block = presentation_list.block(trial_num);

        % Load audio
        PsychPortAudio('FillBuffer', pahandle, [audio_data{stimID}(1,:); audio_data{stimID}(2,:)]);

        % Play audio
        Audio_start_time = PsychPortAudio('Start', pahandle, 1, startCue, waitForStart);       
        presentation_list.onset_TTLtime(trial_num) = Audio_start_time-tstart;

        % Check for presses on Escape key (a bit shorter than trial duration to
        % avoid delays)
        while GetSecs < Audio_start_time + (16*0.250) + ITI - 0.200
        [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                sca; PsychPortAudio('Close', pahandle); error('Escape key pressed. Screen closed.\n%s', 'Abort!')
            end   
        end

        % Compute new start time for follow-up stimulus
        [actualStartTime, ~, ~, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
        if trial_num<numel(stimIDlist)
            if presentation_list.block(trial_num+1) ~= current_block % next trial is a different block
                startCue = estStopTime + ITI + rest_block_dur;  
                tdelay = Audio_start_time - ttime; ttime = Audio_start_time;
                fprintf('\nTrial n°%d: %s, dev pos = %d, delay from previous trial = %0.2f ms', trial_num, stim_list(trial_num), presentation_list.Position_Violation(trial_num), tdelay*1000)
                fprintf('\n   << %d secs rest block >>  ', rest_block_dur)
            else % next trial is in the same block
                startCue = estStopTime + ITI;  
                tdelay = Audio_start_time - ttime; ttime = Audio_start_time;
                fprintf('\nTrial n°%d: %s, dev pos = %d, delay from previous trial = %0.2f ms', trial_num, stim_list(trial_num), presentation_list.Position_Violation(trial_num), tdelay*1000)
            end
        else
            fprintf('\nTrial n°%d: %s, dev pos = %d, delay from previous trial = %0.2f ms', trial_num, stim_list(trial_num), presentation_list.Position_Violation(trial_num), tdelay*1000)
        end
    end
end
finalTime = GetSecs;


catch  % Catch error: this is executed in case something goes wrong in the
       % 'try' part due to programming error etc...
       % if PTB crashes it will land here, allowing us to reset the screen to normal.
       
    disp(' ')
    disp(['CRITICAL ERROR: ' lasterr ])
    disp('Exiting program ...')
    
    fname = fullfile(outdir,['run' num2str(run_number, '%02.f') '_presentation_list_out_' datestr(now,'yy-mm-dd') '.mat']);
    save(fname, 'presentation_list') % save current run with timings
    save(fullfile(outdir,['run' num2str(run_number, '%02.f') '_workspace.mat'])) % save all workspace
    
    % Do same cleanup as at the end of a regular session...
    ShowCursor;
    Screen('CloseAll');
    fclose('all');
    Priority(0);
    PsychPortAudio('Close')
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
    
%     % stop the eye-tracker
%     Eyelink('Message', 'ExperimentAborted'); 
%     Eyelink('StopRecording');
%     Eyelink('CloseFile');
end

%% ======= END EXPERIMENT ======= %%

% % stop the eye-tracker
% Eyelink('Message', 'ExperimentEnd'); %#ok<*UNRCH>
% Eyelink('StopRecording');
% Eyelink('CloseFile');
% try
%     fprintf('\nAsk edf file...');
%     Eyelink('ReceiveFile');
%     if 2==exist(eyedatafile, 'file')
%         fprintf('\nEdf file received\n');
%         fprintf('Data file ''%s'' can be found in ''%s''\n', eyedatafile, pwd );
%         % convert file
% %         EyeData = Edf2Mat(eyedatafile);
% %         save(sprintf('%s%s',outdir,eyedatafile), 'EyeData');
%     end  
% catch
%     fprintf('Problem receiving data file ''%s''\n', eyedatafile );   
% end

Priority(0); % priority back to normal
PsychPortAudio('Close');% Close the audio device:
FlushEvents('keyDown'); % removed typed characters from queue.
fname = fullfile(outdir,['run' num2str(run_number, '%02.f') '_presentation_list_out_' datestr(now,'yy-mm-dd') '.mat']);
save(fname, 'presentation_list') % save current run with timings
save(fullfile(outdir,['run' num2str(run_number, '%02.f') '_workspace.mat'])) % save all workspace
fprintf('\n\n|========================================|')
fprintf('\n|            Run %s complete !            |', num2str(run_number))
fprintf('\n|            Data file saved             |')
fprintf('\n|========================================|')
fprintf('\n        --->> "%s"\n\n', fname)
%fclose(s1);

run_in_progress = {subjectID run_number};
save([rootdir filesep 'run_in_progress.mat'], 'run_in_progress');
% Wait till the end of acquisition (with a blank screen)...
fprintf('\nVolumes acquisition duration: %.2f sec\n', run_duration)
fprintf('Stimuli presentation duration: %.2f sec\n', (finalTime-TTL_time))
fprintf('Waiting the end of acquisition..............')
if no_screen == false
    while GetSecs < (TTL_time+run_duration)
        Screen('Fillrect', w, grey); % blank screen
        Screen('Flip', w);
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            sca; error('Escape key pressed. Screen closed.\n%s', 'Abort!')
        end
    end
end
fprintf('Done !')
ShowCursor; % un-hide the mouse cursor
Screen('closeall'); % this line deallocates the space pointed to by the buffers, and returns the screen to normal.
