function AB16seqfMRI_Generate_SPM_onsets(w, iS)

stim_dur         = 0.050; % item duration
SOA              = 0.250; % stimulus-onset-asynchrony, to get onset of each item
messageDuration  = 3.000; % duration of the block introduction message
responseDuration = 0.0; % onset duration for manual responses


% ======= EXTRACT AND SAVE ONSETS FOR EACH CONDITION IN EACH RUN ======= %
for nrun = 1:numel(w.sessions)

    r = extract_run_info(w.logPath, nrun);

    % =========================================================================== %
    % VERSION 1: sequence*(hab+stand+1dev) - with "overlapping"
    % regressors for dev - NO MANUAL RESPONSE

    % Set onsets list with corresponding names
    onsets_list    = [r.items;           r.items_deviants;           r.message];
    condnames_list = [r.items_ABBA_names;r.items_deviants_ABBA_names;repmat({'Message'},numel(r.message),1);];

    % List of all different item types/conditions occuring in the run
    condlist = string(unique(condnames_list));

    % Onsets for each item type/condition
    onsets = {}; names = {}; durations = {};
    for ii=1:size(condlist,1)
        names{ii}    = condlist{ii};
        onsets{ii}   = onsets_list(find(condnames_list==condlist(ii)));
        if strcmp(names{ii}, 'Message')
            durations{ii} = messageDuration;
        elseif startsWith(names{ii}, 'right') || startsWith(names{ii}, 'left') || startsWith(names{ii}, 'FA') || endsWith(names{ii}, '_resp') || startsWith(names{ii}, 'ButtonPress')
            durations{ii} = responseDuration;
        else
            durations{ii} = stim_dur;
        end
    end

    % Save SPM ONSETS MAT file
    fname = fullfile(w.stimPath,  [w.subjects{iS} '_run_' num2str(nrun, '%02.f') '_onsets_v1.mat']);
    save(fname, 'names' , 'onsets', 'durations')
    disp(['File "' fname '" saved'])


    % =========================================================================== %
    % VERSION 2: sequence*(hab+stand+1dev) - with "overlapping"regressors for DetectedDev - NO MANUAL RESPONSE

    % Set onsets list with corresponding names
    onsets_list    = [r.items;           r.items_detecteddev;           r.message];
    condnames_list = [r.items_ABBA_names;r.items_detecteddev_ABBA_names;repmat({'Message'},numel(r.message),1);];

    % List of all different item types/conditions occuring in the run
    condlist = string(unique(condnames_list));

    % Onsets for each item type/condition
    onsets = {}; names = {}; durations = {};
    for ii=1:size(condlist,1)
        names{ii}    = condlist{ii};
        onsets{ii}   = onsets_list(find(condnames_list==condlist(ii)));
        if strcmp(names{ii}, 'Message')
            durations{ii} = messageDuration;
        elseif startsWith(names{ii}, 'right') || startsWith(names{ii}, 'left') || startsWith(names{ii}, 'FA') || endsWith(names{ii}, '_resp') || startsWith(names{ii}, 'ButtonPress')
            durations{ii} = responseDuration;
        else
            durations{ii} = stim_dur;
        end
    end

    % Save SPM ONSETS MAT file
    fname = fullfile(w.stimPath,  [w.subjects{iS} '_run_' num2str(nrun, '%02.f') '_onsets_v2.mat']);
    save(fname, 'names' , 'onsets', 'durations')
    disp(['File "' fname '" saved'])
end

end

function [run_info] = extract_run_info(logPath, nrun)
% ======= EXTRACT INFO FROM LOG FILE ======= %
% Load log files
cd(logPath)
f = dir(['run' num2str(nrun, '%02.f') '_presentation_list_out*.mat']);
w.logdata = load(f.name);
data = w.logdata.presentation_list;
f = dir(['run' num2str(nrun, '%02.f') '_button_presses*.mat']);
tmp = load(f.name);
responseTime_TTL = tmp.responseTime_TTL;
f = dir(['run' num2str(nrun, '%02.f') '_messages*.mat']);
tmp = load(f.name);
messageTime_TTL = tmp.messageTime_TTL;

% Retrieve which button was pressed
semirun_onsets = []; semirun_buttons = [];
for kk=1:4
    tmp = find(data.semirun == kk);
    semirun_onsets = [semirun_onsets; data.onset_TTLtime(tmp(1))];
    semirun_buttons = [semirun_buttons; data.response_button(tmp(1))];
end
rightClick = []; leftClick = [];
for kk=1:numel(responseTime_TTL)
    semirun = find(responseTime_TTL(kk)>semirun_onsets);
    semirun = semirun(end);
    if semirun_buttons(semirun) == 1
        rightClick = [rightClick; responseTime_TTL(kk)];
    elseif semirun_buttons(semirun) == 2
        leftClick = [leftClick; responseTime_TTL(kk)];
    end
end

% Retrieve TTL-locked onsets of correct responses using "RT" field
zero_rawtime = data.onset_rawtime(1) - data.onset_TTLtime(1);
onset_of_correct_dev = (data.deviant_onset_time(~isnan(data.RT))) - zero_rawtime;
correct_resp = [];
correct_resp = onset_of_correct_dev + (data.RT(~isnan(data.RT))/1000);
RTs = data.RT(~isnan(data.RT));
% Retrieve conditions names for responses using RTs
correct_resp_names = {};
correct_resp_LoTchunkcomp = [];
if numel(RTs) ~= numel(correct_resp); error('Issue with responses...'); end
for irt=1:numel(RTs)
    idx = find(data.RT== RTs(irt));
    correct_resp_names{irt} = [char(data.condition(idx)) '_correct_resp'];
    [sequence_name, seqID, LoTcomp, LoTchunkcomp, pAlt, color] = get_seq_info(data.condition(idx));
    correct_resp_LoTchunkcomp(irt) = LoTchunkcomp;
end

% Classify each click onset as correct or FA
correct_resp_Right = []; FA_resp_Right = [];
for kk=1:numel(rightClick)
    if ismember(rightClick(kk), correct_resp)
        correct_resp_Right = [correct_resp_Right; rightClick(kk)];
    else
        FA_resp_Right = [FA_resp_Right; rightClick(kk)];
    end
end
correct_resp_Left = []; FA_resp_Left = [];
for kk=1:numel(leftClick)
    if ismember(leftClick(kk), correct_resp)
        correct_resp_Left = [correct_resp_Left; leftClick(kk)];
    else
        FA_resp_Left = [FA_resp_Left; leftClick(kk)];
    end
end
% Make sure numbers add up
nCorr = sum(~isnan(data.RT));
nFA = numel(responseTime_TTL) - nCorr;
if ((numel(FA_resp_Right) + numel(FA_resp_Left)) ~= nFA) || ((numel(correct_resp_Right) + numel(correct_resp_Left)) ~= nCorr)
    error('Issue with the classification of correct/FA button presses...')
end
FA_resp = sort([FA_resp_Right; FA_resp_Left]);

% Sound onsets
items = [];
items_names = {}; nitem = 1;
items_deviants = [];
items_deviants_names = {};
items_detecteddev = [];
items_detecteddev_names = {};
dev_LoTchunkcomp = [];
detecteddev_LoTchunkcomp = [];
n_add = 1; n_add2 = 1;
for ntrial=1:size(data,1)
    % fill a 'cond_name'cell for each individual item'
    nitem_pre=nitem;
    while nitem < nitem_pre+16
        if data.block(ntrial) == 1 || data.block(ntrial) == 2
            items_names{nitem,1} = [char(data.condition(ntrial)) '_hab'];
        else
            items_names{nitem,1} = [char(data.condition(ntrial)) '_stand'];
        end
        nitem = nitem + 1;
    end
    % onsets for all items
    stimonsetlist = data.onset_TTLtime(ntrial);
    for ii=2:16
        stimonsetlist(ii) = stimonsetlist(ii-1)+0.250;
    end
    items = [items; stimonsetlist'];
    % additional cond name & onsets for dev items
    devpos = data.Position_Violation(ntrial);
    if devpos ~= 0
        fulldevpos = (ntrial-1)*16 + devpos;
        items_deviants_names{n_add,1} = [char(data.condition(ntrial)) '_dev']; % _pos_' num2str(devpos, '%02.f')];
        items_deviants = [items_deviants; stimonsetlist(devpos)];
        [sequence_name, seqID, LoTcomp, LoTchunkcomp, pAlt, color] = get_seq_info(data.condition(ntrial));
        dev_LoTchunkcomp(n_add) = LoTchunkcomp;
        if ~isnan(data.RT(ntrial))
            items_detecteddev_names{n_add2,1} = [char(data.condition(ntrial)) '_detecteddev']; % _pos_' num2str(devpos, '%02.f')];
            items_detecteddev = [items_detecteddev; stimonsetlist(devpos)];
            detecteddev_LoTchunkcomp(n_add2) = LoTchunkcomp;
            n_add2 = n_add2 + 1;
        end
        n_add = n_add + 1;
    end
end

% Sound onsets names - separating for AB/BA
items_ABBA_names = {}; nitem = 1;
items_deviants_ABBA_names = {}; n_add = 1;
items_detecteddev_ABBA_names = {}; n_add2 = 1;
for ntrial=1:size(data,1)
    % prefix cond name
    if strcmp(data.Presented_sequence(ntrial,1),'0')
        prefix_items_ABBA_names = ['AB_' char(data.condition(ntrial))];
    elseif strcmp(data.Presented_sequence(ntrial,1),'1')
        prefix_items_ABBA_names = ['BA_' char(data.condition(ntrial))];
    end
    % fill a 'items_ABBA_names'cell for each individual item'
    nitem_pre=nitem;
    while nitem < nitem_pre+16
        if data.block(ntrial) == 1 || data.block(ntrial) == 2
            items_ABBA_names{nitem,1} = [prefix_items_ABBA_names '_hab'];
        else
            items_ABBA_names{nitem,1} = [prefix_items_ABBA_names '_stand'];
        end
        nitem = nitem + 1;
    end
    % additional cond name & onsets for dev items
    devpos = data.Position_Violation(ntrial);
    if devpos ~= 0
        fulldevpos = (ntrial-1)*16 + devpos;
        items_deviants_ABBA_names{n_add,1} = [prefix_items_ABBA_names '_dev']; % _pos_' num2str(devpos, '%02.f')];
        n_add = n_add + 1;
        if ~isnan(data.RT(ntrial))
            items_detecteddev_ABBA_names{n_add2,1} = [prefix_items_ABBA_names '_detecteddev']; % _pos_' num2str(devpos, '%02.f')];
            n_add2 = n_add2 + 1;
        end
    end
end

% Sound onsets names - separating for hab block1/block2
items_habblocks_names = {}; nitem = 1;
for ntrial=1:size(data,1)
    % fill a 'cond_name'cell for each individual item'
    nitem_pre=nitem;
    while nitem < nitem_pre+16
        if data.block(ntrial) == 1
            items_habblocks_names{nitem,1} = [char(data.condition(ntrial)) '_hab_block1'];
        elseif data.block(ntrial) == 2
            items_habblocks_names{nitem,1} = [char(data.condition(ntrial)) '_hab_block2'];
        else
            items_habblocks_names{nitem,1} = [char(data.condition(ntrial)) '_stand'];
        end
        nitem = nitem + 1;
    end
end


% ========== Output ========== %
% Sound onsets
run_info.items = items;
run_info.items_names = items_names;
run_info.items_ABBA_names = items_ABBA_names;
run_info.items_habblocks_names = items_habblocks_names;
run_info.items_deviants = items_deviants;
run_info.items_deviants_names = items_deviants_names;
run_info.items_deviants_ABBA_names = items_deviants_ABBA_names;
run_info.items_detecteddev = items_detecteddev;
run_info.items_detecteddev_names = items_detecteddev_names;
run_info.items_detecteddev_ABBA_names = items_detecteddev_ABBA_names;

% Responses
run_info.correct_resp = correct_resp;
run_info.correct_resp_names = correct_resp_names;
run_info.correct_resp_Right = correct_resp_Right;
run_info.correct_resp_Left = correct_resp_Left;
run_info.FA_resp = FA_resp;
run_info.FA_Right = FA_resp_Right;
run_info.FA_Left = FA_resp_Left;
run_info.RTs = RTs;

run_info.rightClick = sort([correct_resp_Right; FA_resp_Right]);
run_info.leftClick  = sort([correct_resp_Left;  FA_resp_Left]);
run_info.ButtonPress= sort([correct_resp; FA_resp]);

% Various
run_info.message = messageTime_TTL;
run_info.detecteddev_LoTchunkcomp = detecteddev_LoTchunkcomp;
run_info.dev_LoTchunkcomp = dev_LoTchunkcomp;

end
