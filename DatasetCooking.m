for subject=1:8
clearvars -except subject*
subjects=subject;
load(sprintf('%s/P300Dataset/p300-subject-%02d.mat',getdatasetpath(),subject));
subject=subjects;
clear a
clear ab
clear b
clear c
clear ans
clear channels
clear counterhits
clear counternohits
clear d
clear duration
clear durs
clear e
clear endd
clear endset
clear f
clear h
clear i
clear idx
clear idxset
clear j
clear n
clear samplingFreq
clear sls
clear ss
clear start
clear startset
clear stimulations
clear stimuls
clear stopsample
clear stoptime
clear subjectartifacts
clear subjects
clear subjectaverages
clear subjectsingletriality
clear targets
clear total
clear trial
clear validstimuls
clear vl
clear z
save(sprintf('%s/P300Dataset/P300S%02d.mat',getdatasetpath(),subject));
end
