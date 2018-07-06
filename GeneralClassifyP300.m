
for subject=subjectRange   
    epochRange=1:epoch;
    trainingRange = 1:nbofclassespertrial*15;
    testRange=nbofclassespertrial*15+1:min(nbofclassespertrial*35,epoch);
    
    %trainingRange=1:nbofclasses*35;
    
    SBJ(subject).F = F;
    SBJ(subject).epochRange = epochRange;
    SBJ(subject).labelRange = labelRange;
    SBJ(subject).trainingRange = trainingRange;
    SBJ(subject).testRange = testRange;
    

    switch classifier
        case 5
            for channel=channelRange
                [DE(channel), ACC, ERR, AUC, SC(channel)] = LDAClassifier(F,labelRange,trainingRange,testRange,channel);
                globalaccij1(subject,channel)=ACC;
                globalsigmaaccij1 = globalaccij1;
                globalaccij2(subject,channel)=AUC;
            end  
        case 4
            for channel=channelRange
                [DE(channel), ACC, ERR, AUC, SC(channel)] = SVMClassifier(F,labelRange,trainingRange,testRange,channel);
                globalaccij1(subject,channel)=ACC;
                globalsigmaaccij1 = globalaccij1;
                globalaccij2(subject,channel)=AUC;
            end            
        case 1
            for channel=channelRange
                [DE(channel), ACC, ERR, AUC, SC(channel)] = NNetClassifier(F,labelRange,trainingRange,testRange,channel);
                globalaccij1(subject,channel)=ACC;
                globalsigmaaccij1 = globalaccij1;
                globalaccij2(subject,channel)=AUC;
            end
        case 2
            [AccuracyPerChannel, SigmaPerChannel] = CrossValidated(F,epochRange,labelRange,channelRange, @IterativeNBNNClassifier,1);
            globalaccij1(subject,:)=AccuracyPerChannel
            globalsigmaaccij1(subject,:)=SigmaPerChannel;
            globalaccijpernumberofsamples(globalnumberofepochs,subject,:) = globalaccij1(subject,:);
        case 3
            for channel=channelRange
                [DE(channel), ACC, ERR, AUC, SC(channel)] = IterativeNBNNClassifier(F,channel,trainingRange,labelRange,testRange,false,false);

                globalaccij1(subject,channel)=1-ERR/size(testRange,2);
                globalaccij2(subject,channel)=AUC;
                globalsigmaaccij1 = globalaccij1;
            end
        case 6
            for channel=channelRange
                DE(channel) = NBNNFeatureExtractor(F,channel,trainingRange,labelRange,[1 2],false); 
   
                %[ACC, ERR, AUC, SC(channel)] = NBMultiClass(F,DE(channel),channel,testRange,labelRange,false);
                [ACC, ERR, AUC, SC(channel)] = NBNNClassifier4(F,DE(channel),channel,testRange,labelRange,false,distancetype,k);                                                        
                
                globalaccij1(subject,channel)=1-ERR/size(testRange,2);
                globalaccij2(subject,channel)=AUC;
                globalsigmaaccij1 = globalaccij1;
            end

    end
    SBJ(subject).DE = DE;
    SBJ(subject).SC = SC;
end


%%
for subject=subjectRange
    % '2'    'B'    'A'    'C'    'I'    '5'    'R'    'O'    'S'    'E'    'Z'  'U'    'P'    'P'    'A'   
    % 'G' 'A' 'T' 'T' 'O'    'M' 'E' 'N''T' 'E'   'V''I''O''L''A'  'R''E''B''U''S'
    Speller = SpellMe(F,channelRange,16*nbofclassespertrial/12:35*nbofclassespertrial/12+(nbofclassespertrial/12-1),labelRange,trainingRange,testRange,SBJ(subject).SC);

    S = 'MANSOCINCOJUEGOQUESO';
    S = repmat(S,nbofclassespertrial/12);
    S = reshape( S, [1 size(S,1)*size(S,2)]);
    S=S(1:size(S,2)/(nbofclassespertrial/12));
    
    SpAcc = [];
    for channel=channelRange
        counter=0;
        for i=1:size(S,2)
            if Speller{channel}{i}==S(i)
                counter=counter+1;
            end
        end
        spellingacc = counter/size(S,2);
        SpAcc(end+1) = spellingacc;
        globalspeller(subject,channel) = spellingacc;
        %globalspeller(subject,channel,globalrepetitions) = spellingacc;
    
    end
    [a,b] = max(SpAcc);
end

graphics=true;
if (graphics)
    %%
    subject=1;
    channel=2;
    SC=SBJ(subject).SC(channel);
    ML=SBJ(subject).DE(channel);
    F=SBJ(subject).F;

    for i=1:30
        figure;DisplayDescriptorImageFull(F,subject,ML.C(2).IX(i,3),ML.C(2).IX(i,2),ML.C(2).IX(i,1),ML.C(2).IX(i,4),false);
    end

    %%

    figure('Name','Class 2 P300','NumberTitle','off');
    setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
    fcounter=1;
    for i=1:30
        ah=subplot_tight(6,5,fcounter,[0 0]);
        DisplayDescriptorImageFull(F,subject,ML.C(2).IX(i,3),ML.C(2).IX(i,2),ML.C(2).IX(i,1),ML.C(2).IX(i,4),true);
        fcounter=fcounter+1;
    end
    figure('Name','Class 1','NumberTitle','off');
    setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
    fcounter=1;
    for i=1:30
        ah=subplot_tight(6,5,fcounter,[0 0]);
        DisplayDescriptorImageFull(F,subject,ML.C(1).IX(i,3),ML.C(1).IX(i,2),ML.C(1).IX(i,1),ML.C(1).IX(i,4),true);
        fcounter=fcounter+1;
    end
    [TM, TIX] = BuildDescriptorMatrix(F,channel,labelRange,testRange(labelRange(testRange)==2));
    fcounter=1;
    figure('Name','P300 Query','NumberTitle','off');
    setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
    for i=1:30
        ah=subplot_tight(6,5,fcounter,[0 0]);
        DisplayDescriptorImageFull(F,subject,TIX(i,3),TIX(i,2),TIX(i,1),TIX(i,4),true);
        fcounter=fcounter+1;
    end
    figure('Name','P300 Query (resto)','NumberTitle','off');
    setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
    fcounter=1;
    for i=30:40
        ah=subplot_tight(2,5,fcounter,[0 0]);
        DisplayDescriptorImageFull(F,subject,TIX(i,3),TIX(i,2),TIX(i,1),TIX(i,4),true);
        fcounter=fcounter+1;
    end

    %%
    experiment=sprintf('Hellinger. Butter de 3 a 4, K = %d, upsampling a 16, zscore a 3,NBNN con artefactos, cosine float, unweighted without artifacts ',k);
    fid = fopen('experiment.log','a');
    fprintf(fid,'Experiment: %s \n', experiment);
    fprintf(fid,'st %f sv %f scale %f timescale %f qKS %d\n',siftscale(1),siftscale(2),imagescale,timescale,qKS);
    %totals = DisplayTotals(subjectRange,globalaccij1,globalspeller,globalaccij2,globalspeller,channels)
    %totals(:,6)
    fclose(fid)
    %%
    DisplayDescriptorImageFull(F,1,2,1,1,1,false);
    %%
    figure
    hold on
    for i=1:size(ML.C(2).M,2)
        plot(ML.C(2).M(:,i),'x');
    end
    hold off
    figure
    hold on
    pp=randperm(size(ML.C(1).M,2),size(ML.C(2).M,2));
    for i=1:size(ML.C(2).M,2)
        plot(ML.C(1).M(:,pp(i)),'x');
    end
    hold off



    DisplayDescriptorImageFull(F,1,1,1,1,-1,true)


    %%

    for l=1:20
        string=[];
        for c=1:8
            string= [string Speller{c}{l}];
        end
        string
        histogram=zeros(1,50);
        for n=1:length(string)
            currentLetter=string(n);
            histogram(currentLetter-47)=histogram(currentLetter-47)+1;
        end
        [val, ensembleletter] = max(histogram);
        Speller{9}{l} = char(ensembleletter+47);
    end
end











