clear all;

load('p300.mat');

channels={ 'Fz'  ,  'Cz',    'P3' ,   'Pz'  ,  'P4'  , 'PO7'   , 'PO8',  'Oz'};
windowsize=1;
downsize=10;
imagescale=2;
timescale=4;
amplitude=2;

sqKS=[14*4];
siftscale=[2 2];
siftdescriptordensity=1;
minimagesize=floor(sqrt(2)*15*siftscale(2)+1);
minimagesize=minimagesize*2;

nbofclassespertrial=12;
k=7;
adaptative=false;
subjectRange=1:1;
distancetype='cosine';
applyzscore=true;
featuretype=1;
classifier=6;

%SVM
%featuretype=2;
%timescale=1;
%applyzscore=false;
%classifier=4;


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
            assert( rcounter{subject}{trial}{i} == 10 );
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
                rsignal{i} = rsignal{i}*amplitude;
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
                %minimagesize=1;
                [eegimg, DOTS, zerolevel] = eegimage2(channel,rsignal{i},imagescale,1, false,minimagesize);
                %siftscale(1) = 11.7851;
                %siftscale(2) = (height-1)/(sqrt(2)*15);
                saveeegimage(subject,epoch,label,channel,eegimg);
                zerolevel = size(eegimg,1)/2;

    %             if ((size(find(trainingRange==epoch),2)==0))
    %                qKS=ceil(0.20*(Fs)*timescale):floor(0.20*(Fs)*timescale+(Fs)*timescale/4-1);
    %             else
                    qKS=sqKS(subject);
                    %qKS=125;
    %             end
    
                qKS=[qKS+reshape(repmat(-11:11,[23 1]),[1 23*23]);zerolevel+repmat(-11:11,[1 23])]';

                [frames, desc] = PlaceDescriptorsByFrames(eegimg, siftscale,qKS,distancetype);
                
                
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



%%
epochRange=1:epoch;
trainingRange = 1:nbofclassespertrial*15;
testRange=nbofclassespertrial*15+1:min(nbofclassespertrial*35,epoch);

channel=7;
hits = find(labelRange==2);
M = [];
% Agarro el descriptor del centro.
M = [M F(channel,2,hits(1)).descriptors(:,265)];
for i=2:30
    desc = F(channel,2,hits(i)).descriptors(:,:);
    [Z,I] = pdist2(desc',M',distancetype,'Smallest',1 );
    M = [M desc(:,I(1))];
end

TM = [];
for test=testRange
    desc = F(channel,labelRange(test),test).descriptors(:,:);
    [Z,I] = pdist2(desc',M',distancetype,'Smallest',1 );
    TM = [TM desc(:,I(1))];
end
kparam=7;
DE.C(2).M = M;
fdsds
%DE = NBNNFeatureExtractor(F,channel,trainingRange,labelRange,[1 2], false);

% Este metodo toma en consideracion que la clasificacion binaria se usa
% para un speller de p300.

assert( mod(size(testRange,2),12)==0, 'This method only works for P300 spellers');

mind = 1;
maxd = 6;

SC.CLSF = {};
predicted=[];
score=[];

% W contiene los pesos de los descriptores de la bolsa de hit
% K = size(DE.C(2).M,2);
% 
% D=[];
% for i=1:2:K
%    ni=floor(i/2)+1;
%    Z= pdist2(DE.C(1).M(:,(ni-1)*10+1:(ni-1)*10+10)',DE.C(2).M(:,i:i+1)',distancetype);
%    %Z= pdist2(DE.C(1).M(:,i:i+1)',DE.C(2).M(:,i:i+1)','euclidean');
%    Di = sum(Z) ;
%    D(end+1)=Di(1);
%    D(end+1)=Di(2);
% end
% 
% DR=(D-min(D))/range(D);
% 
% Wdi = normpdf(DR,0,1);


for f=1:size(testRange,2)/12
    
    K = size(DE.C(2).M,2);

    [Z,I] = pdist2(DE.C(2).M',(TM(:,mind:maxd+6)'),distancetype,'Smallest',K );
    
    k = kparam;

    %Wi = Wdi(I(1:k,1:6)) ./ repmat( sum(Wdi(I(1:k,1:6))),k,1);     
    Wi = ones(k,6);
    if (k==1)
        sumsrow = Z(1:k,1:6).*Wi(1:k,1:6);
    else
        sumsrow = dot(Z(1:k,1:6),Wi(1:k,1:6));
    end
    
    
    %Wi = Wdi(I(1:k,7:12)) ./ repmat( sum(Wdi(I(1:k,7:12))),k,1); 
    Wi = ones(k,6);
    if (k==1)
        sumscol = Z(1:k,7:12).*Wi(1:k,1:6);
    else
        sumscol = dot(Z(1:k,7:12),Wi(1:k,1:6));
    end

    % Me quedo con aquel que la suma contra todos, dio menor.
    [c, row] = min(sumsrow);
    [c, col] = min(sumscol);
    %col=col+6;

    % I(1:3,1:6) Me da en cada columna los ids de los descriptores de M mas
    % cercaos a cada uno de los descriptores de 1 a 6.
    
    assert( sum(sumsrow)>0, 'Problem with distance function for this feature.');
    assert( sum(sumscol)>0, 'Problem with distance function for this feature.');
    

    % Las predicciones son 1 para todos excepto para row y col.
    for i=1:6
        if (i==row)
            predicted(end+1) = 2;
        else
            predicted(end+1) = 1;
        end
        score(end+1) = 1-sumsrow(i)/sum(sumsrow);
    end
    for i=1:6
        if (i==col)
            predicted(end+1) = 2;
        else
            predicted(end+1) = 1;
        end
        score(end+1) = 1-sumscol(i)/sum(sumscol);
    end

    mind=mind+12;
    maxd=maxd+12;
end
score=score';

%for channel=channelRange
fprintf ('Channel %d -------------\n', channel);

%M = MM(channel).M;
%IX = MM(channel).IX;

expected = labelRange(testRange);


%predicted=randi(unique(labelRange),size(expected))

C=confusionmat(expected, predicted)


%if (C(1,1)+C(2,2) > 65)
%    error('done');
%end

%[X,Y,T,AUC] = perfcurve(expected,single(predicted==2),2);
[X,Y,T,AUC] = perfcurve(expected,score,2);

%figure;plot(X,Y)
%xlabel('False positive rate')
%ylabel('True positive rate')
%title('ROC for Classification of P300')

ACC = (C(1,1)+C(2,2)) / size(predicted,2);
ERR = size(predicted,2) - (C(1,1)+C(2,2));

SC.FP = C(2,1);
SC.TP = C(2,2);
SC.FN = C(1,2);
SC.TN = C(1,1);

[ACC, (SC.TP/(SC.TP+SC.FP))]

SC.expected = expected;
SC.predicted = predicted;    

SD(channelRange) = SC;

% '2'    'B'    'A'    'C'    'I'    '5'    'R'    'O'    'S'    'E'    'Z'  'U'    'P'    'P'    'A'   
% 'G' 'A' 'T' 'T' 'O'    'M' 'E' 'N''T' 'E'   'V''I''O''L''A'  'R''E''B''U''S'
Speller = SpellMe(F,channelRange,16*nbofclassespertrial/12:35*nbofclassespertrial/12+(nbofclassespertrial/12-1),labelRange,trainingRange,testRange,SD);

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

