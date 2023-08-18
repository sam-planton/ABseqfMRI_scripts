function sequences_and_violations = sequences_and_violations_list()

    % Sequence + list of positions for violations
    % 2 versions per sequence ("ab" & "ba")
    s1    = [{'AAAAAAAAAAAAAAAA'-'A'}     ,{[9,  12, 13, 15]}];
    s2    = [{'ABABABABABABABAB'-'A'}     ,{[9,  12, 14, 15]}];
    s3    = [{'AABBAABBAABBAABB'-'A'}     ,{[10, 11, 14, 15]}];
    s4    = [{'AAAABBBBAAAABBBB'-'A'}     ,{[9,  12, 13, 15]}];
    s5    = [{'AABBABABAABBABAB'-'A'}     ,{[10, 11, 14, 15]}];
    s6    = [{'AAAABBBBAABBABAB'-'A'}     ,{[10, 11, 14, 15]}];
    s7    = [{'ABAAABBBBABBAAAB'-'A'}     ,{[9,  12, 14, 15]}];
    s1bis = [{abs('AAAAAAAAAAAAAAAA'-'B')},{[9,  12, 13, 15]}];
    s2bis = [{abs('ABABABABABABABAB'-'B')},{[9,  12, 14, 15]}];
    s3bis = [{abs('AABBAABBAABBAABB'-'B')},{[10, 11, 14, 15]}];
    s4bis = [{abs('AAAABBBBAAAABBBB'-'B')},{[9,  12, 13, 15]}];
    s5bis = [{abs('AABBABABAABBABAB'-'B')},{[10, 11, 14, 15]}];
    s6bis = [{abs('AAAABBBBAABBABAB'-'B')},{[10, 11, 14, 15]}];
    s7bis = [{abs('ABAAABBBBABBAAAB'-'B')},{[9,  12, 14, 15]}];
    
%     sequences_and_violations = [s1;s2;s3;s4;s5;s6;s7;s1bis;s2bis;s3bis;s4bis;s5bis;s6bis;s7bis]; % Human MEG
    sequences_and_violations = [s1;s2;s3;s4;s5;s6;s7;s1bis;s2bis;s3bis;s4bis;s5bis;s6bis;s7bis]; % Human fMRI
%     sequences_and_violations = [s2;s3;s5;s2bis;s3bis;s5bis]; % Human ECOG ?
%     sequences_and_violations = [s2;s3;s5;s2bis;s3bis;s5bis]; % Marmoset

    rundur = (((16*250) + 500)*46)/60000;
    alldur = rundur*size(sequences_and_violations,1);
    disp(['With run duration of ' num2str(rundur) ' min, total duration with no pauses will be ' num2str(alldur, '%0.1f') ' min'])
end