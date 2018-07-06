% SVM Full Classification

clear mex
clc
clear
clearvars
close all

%run('C:/vlfeat/toolbox/vl_setup')
rng(396544);

subjectaverages= cell(0);
subjectartifacts = 0;
subjectsingletriality=119;

%for subjectsingletriality=12*[10:-3:1]-1
for subject = 1:1
clear mex;clearvars  -except subject*;close all;clc;

% Clean all the directories where the images are located.
cleanimagedirectory();

%subject = 2;
% 3 letters
%convert_ov2mat('C:/Workspace/GuessMe/signals/p300-train-[2017.05.09-15.46.13].ov','C:/Workspace/GuessMe/signals/p300-train.mat')
% full 7 letters meaningless.
%convert_ov2mat('C:/Workspace/GuessMe/signals/p300-train-[2017.09.26-14.59.43].ov','C:/Workspace/GuessMe/signals/p300-train.mat')
load('./signals/p300-train.mat');

% Clean EEG image directory
if (exist(sprintf('%s','sidgnals\\*.ov'),'dir'))
    delete(sprintf('%s%s*.ov','signals\\',filesep));
end

if (exist('p300.mat','file'))
    delete('p300.mat');
end

% NN.NNNNN
% data.X(sample, channel)
% data.y(sample)  --> 0: no, 1:nohit, 2:hit
% data.y_stim(sample) --> 1-12, 1-6 cols, 7-12 rows

%     'Fz'    'Cz'    'Pz'    'Oz'    'P3'    'P4'    'PO7'    'PO8'
channelRange=1:8;


samples(find(sampleTime==stims(1,1)),:)

a=find(stims(:,2)==hex2dec('00008205')); % 32779
b=find(stims(:,2)==hex2dec('00008206')); % 32780

%%
fprintf('%04x\n',stims(:,2))

% First, get the stimulus to different events
c=find(stims(:,2)==hex2dec('00008005')); % 32773 Trial Start
d=find(stims(:,2)==hex2dec('00008006')); % 32774 Trial Stop
a=find(stims(:,2)==hex2dec('00008205')); % 33285 Hit
b=find(stims(:,2)==hex2dec('00008206')); % 32286 Nohit
e=find(stims(:,2)==hex2dec('0000800C')); % Visual Stimulus Stop
f=find(stims(:,2)==hex2dec('0000800B')); % Visual Stimulus Start

% Find all the stimulus assosiated with row/col flashing.
stimuls = [];
for i=1:12
    stimuls = [stimuls; find(stims(:,2)==33025-1+i)];
end


%%
% Chequear si la cantidad de estimulos encontradas coincide.
total=0;
for i=1:12
    size(find(stims(:,2)==33025-1+i))
    total=total+size(find(stims(:,2)==33025-1+i))
end


% Los stimulos pueden venir despues de muchas cosas.
% Filtrar solo aquellos estimulos que estan asociados a targets.
counterhits=0;
counternohits=0;
validstimuls=[];
for i=1:size(stimuls,1)
    vl=stims(stimuls(i)-1,1:2)
    if (vl(2) == 33285) % Hit
        counterhits = counterhits + 1;
        validstimuls(end+1) = stimuls(i);
    elseif (vl(2) == 33286) % Nohit
        counternohits = counternohits + 1;
        validstimuls(end+1) = stimuls(i);
    end
    assert ( vl(2)==33285 || vl(2)==33286 || vl(2)==32777 || vl(2) == 897 || vl(2)>=33025 || vl(2)<=33036);
end

% Los que valen son los que estan precedidos por una marca de target o no
% target
% Chequear si los targets estan bien asignados a los mismos estimulos
% dentro del mismo trial.
%%
for trial=1:35
    h=[];
    for i=1:20
        vl=stims(a((trial-1)*20+i)+1,1:2);
        [(trial-1)*35+i vl(2)];
        h=[h vl(2)];
    end
    h = unique(h);
    h
    % Verificar que para cada trial, solo haya dos tipos de estimulos
    % asociados a hit (el correspondiente a las filas y el de las columnas)
    assert( size(h,2) == 2);
end


%%
ab = [a; b];

% a hits, b nohits, c and d contain where trial end and stop (five of each
% per letter).

ab = sort(ab);

% Cut the stimulus from stims, getting only the time and duration of each.
targets = [ stims(ab,1:3)];

% Remap targets, assigning 1 for NoHit and 2 for Hit.
targets(targets(:,2)==33285,2) = 2;
targets(targets(:,2)==33286,2) = 1;
targets(targets(:,2)==32773,2) = 0;
targets(targets(:,2)==32774,2) = 0;
targets(targets(:,2)==32780,2) = 0;



% Sort validstimuls based on time.
sls = sort(validstimuls);

% Pick the stimuls sorted.
stimulations = [ stims(sls,1:3) ];

% Map stimulus to 1-12 values.
stimulations(:,2) = stimulations(:,2) - 33025 + 1;
stimulations( stimulations(:,2) < 0) = 0; 

% trials
z = stims(c,1);

% Stop time is where the first invalid trial starts.
if (size(c,1)>35)
    
    stoptime=stims(c(36),1);

    stopsample=find(sampleTime>stims(c(36),1));


    sampleTime(stopsample(1):end,:) = [];
    samples(stopsample(1):end,:) = [];

    z(36) = [];

    targets(4201:end,:) = [];
    stimulations(4201:end,:) = [];
end
    
    
% Check target consistency
Word = [];
for trial=1:35
    h=[];
    for i=1:120
        if (targets((trial-1)*120+i,2)==2)
            h = [h stimulations((trial-1)*120+i,2)];
        end
    end
    % There must be only TWO targets per trial (row and col).
    h = unique(h);
    assert( size(h,2) == 2);
    Word = [Word SpellMeLetter(h(1),h(2))];
end

Word

% Data Structure
data = cell(0);

data.X = samples;
data.y = zeros(size(samples,1),1);
data.y_stim = zeros(size(samples,1),1);
data.trial=zeros(5,1);

Fs=250;

data.flash = [];

for i=1:size(targets,1)
    % Obtengo el ID del sample que esta justo despues del estimulo
    maximalsampleidx=find(sampleTime>=targets(i,1));
    maximalsampleidx=maximalsampleidx(1);
    
    % Obtengo la localizacion donde esta el marcador de fin del estimulo
    % (e)
    loc = find(stims(e,1)>targets(i,1));
    loc = loc(1); % Location on e.
    duration = stims(e(loc),1)-targets(i,1);
    
    % Marco donde inicia el flash y la duracion en sample points.
    data.flash(end+1,1) = maximalsampleidx 
    data.flash(end,2) = ceil(Fs*duration);
    
    %fakeEEG=fakeeegoutput(4,targets(i,2),channelRange,25,100,4);
    
    % Marco todos los estimulos y targets donde el flash estuvo presente.
    for j=1:ceil(Fs*duration)
        data.y(maximalsampleidx+j-1) = targets(i,2); 
        data.y_stim(maximalsampleidx+j-1) = stimulations(i,2);
        
        %fakeEEG(j,:);
    end
    
    %if (targets(i,2)==2)
    %    data.X(maximalsampleidx+1-1:maximalsampleidx+1-1+ceil(Fs*0.33),:) = zeros(ceil(Fs*0.33)+1,size(data.X,2));
    %end    
    
    
end


% Marco los inicios de los trials.
for i=1:size(z)
    n=find(sampleTime>z(i));
    data.trial(i)=n(1);
end

data.trial = data.trial';


%%
% Antes de cada uno de los inicios de los flash, los estimulos tienen que
% estar marcados con zero.
for i=1:4200
    ss=data.y_stim(data.flash(i)-5:data.flash(i)+40)'
    
    assert ( ss(5) == 0, 'Not zero');
end
%%

%data.X = data.X * 10;
save('p300.mat');

% LISTOOOOOO
 
end

%%
load('./signals/p300-subject-01.mat');

channels={ 'Fz'  ,  'Cz',    'P3' ,   'Pz'  ,  'P4'  , 'PO7'   , 'PO8',  'Oz'};
windowsize=1;
downsize=10;
imagescale=2;
timescale=4;
amplitude=3;
sqKS=[44];
siftscale=[2 2];
siftdescriptordensity=1;
minimagesize=floor(sqrt(2)*15*siftscale(2)+1);
nbofclassespertrial=12;
k=7;
adaptative=false;
subjectRange=1:1;
distancetype='cosine';
applyzscore=false;

%SVM
featuretype=2;
timescale=1;
applyzscore=false;



for subject=1:1
    
EEG = prepareEEG(Fs,windowsize,downsize,120,1:1,1:8);
Fs=ceil(Fs/downsize);

for subject=1:1
    for trial=1:35
        for i=1:12 rcounter{subject}{trial}{i} = 0; end
        for flash=1:120
            rcounter{subject}{trial}{EEG(subject,trial,flash).stim} = rcounter{subject}{trial}{EEG(subject,trial,flash).stim}+1;
        end
        % Check if all the epochs contain 10 repetitions.
        for i=1:12
            %assert( rcounter{subject}{trial}{i} == 10 );
        end
    end
end


%%
% Build routput pasting epochs toghether...
clear hit
for subject=1:1
    for trial=1:35
        for i=1:12 hit{subject}{trial}{i} = 0; end
        for i=1:12 routput{subject}{trial}{i} = []; end
        for flash=1:120
            output = EEG(subject,trial,flash).EEG;
            routput{subject}{trial}{EEG(subject,trial,flash).stim} = [routput{subject}{trial}{EEG(subject,trial,flash).stim} ;output];
            hit{subject}{trial}{EEG(subject,trial,flash).stim} = EEG(subject,trial,flash).label;
        end
    end
end

%%
h=[];
Word=[];
for subject=1:1
    for trial=1:35
        hh = [];
        for i=1:12
            rput{i} = routput{subject}{trial}{i};
            channelRange = (1:size(rput{i},2));
            channelsize = size(channelRange,2);

            assert( size(rput{i},1)/(Fs*windowsize) == rcounter{subject}{trial}{i}, 'Something wrong with PtP average. Sizes do not match.');

            rput{i}=reshape(rput{i},[(Fs*windowsize) size(rput{i},1)/(Fs*windowsize) channelsize]); 

            for channel=channelRange
                rmean{i}(:,channel) = mean(rput{i}(:,:,channel),2);
            end

            if (hit{subject}{trial}{i} == 2)
                h = [h i];
                hh = [hh i];
            end    
            routput{subject}{trial}{i} = rmean{i};
        end
        Word = [Word SpellMeLetter(hh(1),hh(2))];
    end
end

clear rsignal
for subject=1:1
    for trial=1:35
        
        for i=1:12

            rmean{i} = routput{subject}{trial}{i};
            
            for c=channelRange
                %rsignal{i}(:,c) = resample(rmean{i}(:,c),size(rmean{i},1)*timescale,size(rmean{i},1));
                rsignal{i}(:,c) = resample(rmean{i}(:,c),1:size(rmean{i},1),timescale);
            end

            if (applyzscore)
                rsignal{i} = zscore(rsignal{i})*amplitude;
            else
                rsignal{i} = rsignal{i};
            end
            
            routput{subject}{trial}{i} = rsignal{i};
        end
    end
end


%%

if (featuretype == 1)
    epoch=0;
    labelRange=[];
    epochRange=[];
    stimRange=[];
    for subject=1:1
        for trial=1:35        
            for i=1:12
            epoch=epoch+1;    
            label = hit{subject}{trial}{i};
            labelRange(epoch) = label;
            stimRange(epoch) = i;
            DS = [];
            rsignal{i}=routput{subject}{trial}{i};
            for channel=channelRange
                [eegimg, DOTS, zerolevel] = eegimage(channel,rsignal{i},imagescale,1, false,minimagesize);

                saveeegimage(subject,epoch,label,channel,eegimg);
                zerolevel = size(eegimg,1)/2;

    %             if ((size(find(trainingRange==epoch),2)==0))
    %                qKS=ceil(0.20*(Fs)*timescale):floor(0.20*(Fs)*timescale+(Fs)*timescale/4-1);
    %             else
                    qKS=sqKS(subject);
    %             end

                [frames, desc] = PlaceDescriptorsByImage(eegimg, DOTS,siftscale, siftdescriptordensity,qKS,zerolevel,false,distancetype);            
                F(channel,label,epoch).stim = i;
                F(channel,label,epoch).hit = hit{subject}{trial}{i};


                F(channel,label,epoch).descriptors = desc;
                F(channel,label,epoch).frames = frames; 
            end
            end
        end
    end
else
    epoch=0;
    labelRange=[];
    epochRange=[];
    stimRange=[];
    for subject=1:1
        for trial=1:35        
            for i=1:12
                epoch=epoch+1;    
                label = hit{subject}{trial}{i};
                labelRange(epoch) = label;
                stimRange(epoch) = i;
                DS = [];
                rsignal{i}=routput{subject}{trial}{i};

                feature = [];

                for channel=channelRange
                    feature = [feature ; rsignal{i}(:,channel)];
                end  

                for channel=channelRange
                    F(channel,label,epoch).hit = hit{subject}{trial}{i};
                    F(channel,label,epoch).descriptors = feature;
                    F(channel,label,epoch).frames = [];   
                    F(channel,label,epoch).stim = i;
                end    
            end
        end
    end
end
 
end


classifier=4;

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







