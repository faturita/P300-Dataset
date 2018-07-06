clear all
close all
globalrepetitions=10;
KS=10:50;
for globaliterations=1:size(KS,2)
    run('OfflineProcessP300.m');
end




%%
KKS=[];
for subject=subjectRange
    performances=reshape(globalspellers(subject,:,:),[8 41]);
    figure;plot(KS,max(performances))
    
    [~,best] = max(max(performances));
    best = best(1);
    
    [ChAcc,ChNum] = max(globalspellers(subject,:,best));
    
    KKS = [KKS KS(best)];
end

KKS
