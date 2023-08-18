function Sig=waitTTL(deviceIndex) 

if nargin < 1
    deviceIndex = [];
end

clear PsychHID; 
[keyboardIndices, productNames, allInfos] = GetKeyboardIndices;
[logicalTrig,locationTrig]=ismember({'Current Designs, Inc. TRIGI-USB'},productNames);% trigger device
[logicalButt,locationButt]=ismember({'Arduino LLC Arduino Leonardo'},productNames);% 2-button device
[logicalKey,locationKey]=ismember({'Dell Dell USB Keyboard'},productNames);% PC keyboard
devicenumtrigger=allInfos{locationTrig}.index;
devicenum=allInfos{locationButt}.index;
devicenumkey=allInfos{locationKey}.index;


% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems:
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
TTLKey = KbName('t');
keysOfInterest=zeros(1,256);
keysOfInterest(TTLKey)=1;  %Wait for the "S" key with KbQueueWait. 
keysOfInterest(escapeKey)=1;

% fprintf('\n\n|========================================|')
% fprintf('\n|                                        |')
% fprintf('\n|                                        |')
% fprintf('\n|                                        |')
% fprintf('\n|             Waiting TTL...             |')
% fprintf('\n|       (or press Esc to escape)         |')
% fprintf('\n|                                        |')
% fprintf('\n|                                        |')
% fprintf('\n|                                        |')
% fprintf('\n|========================================|\n\n')

fprintf('\n\n....................................................................................................................')
fprintf('\n....................................................................................................................')
fprintf('\n........;@@r...r@@:...@@r............@@@:..i@@,..@@@:...........................r@@@@@@@@@@@@@@@@@.i@@r.............')
fprintf('\n........;@@@:.,#@@H..i@@S............SSS;.:@@@;..;SSi...................5@@i....h@@@@@@@@@@@@@@@@@,#@@H.............')
fprintf('\n.........G@@G.5@@@@;,B@#,.;X@@@@@9:.,#@@iS@@@@@@:2@@#.,@@@9@@@#s..iB@@@@@H........,@@@i....:@@@;...#@@A.............')
fprintf('\n.........;@@@;#@&@@Hs@@i.,SSi,.i@@#.,@@@i.:@@@;..5@@#.,@@@5.5@@@,s@@#..#@@i.......,@@@i....:@@@;...#@@H.............')
fprintf('\n..........G@@@@@,A@@@@M,..rS@@@@@@#.,@@@i.:@@@;..2@@#.,@@@i.5@@@,.G@@@@@@S,.......,@@@i....:@@@;...#@@H.............')
fprintf('\n..........;@@@@2.;@@@@r..i@@@s.5@@#.,@@@i.:@@@;..5@@#.,#@@i.5@@#,3@@&SSSSS:.......,@@@i....:@@@;...#@@#SSSSr........')
fprintf('\n...........&@@@,..G@@M...,h@@@MX@@#,,@@@i..S@@@@:5@@#.,@@@i.2@@@,2@@@@@@@@@;......,#@@i....:@@@;...#@@@@@@@9........')
fprintf('\n................................................................:@@#SSSSM@@;........................................')
fprintf('\n..................................................................r@@@@@@r..........................................')
fprintf('\n........................................................................................(or press Esc to escape)....\n\n')



KbQueueCreate(devicenumtrigger, keysOfInterest);
KbQueueStart(devicenumtrigger);

while 1
    % Check the queue for key presses.
    [ pressed, firstPress]=KbQueueCheck(devicenumtrigger);
    if pressed
        if firstPress(TTLKey)
            %display(firstPress(FIND(firstPress)));
            %display(KbName(firstPress));
            %display('TTL')
            disp('!~~ TTL ~~!');
            KbQueueFlush(devicenumtrigger);
            KbQueueRelease(devicenumtrigger);
            Sig = 1;
            break;
        end
        if firstPress(escapeKey)
            KbQueueRelease(devicenumtrigger);
            Sig = -1;
            break
        end
    end
end
return;