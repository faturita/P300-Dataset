%load('p300.mat');

channels={ 'Fz'  ,  'Cz',    'P3' ,   'Pz'  ,  'P4'  , 'PO7'   , 'PO8',  'Oz'};
windowsize=1;

downsize=40;
imagescale=1*6;
timescale=1*6;
amplitude=1;

sqKS=[18];

downsize=10;
imagescale=2;
timescale=4;
amplitude=3;
sqKS=[44]+zeros(1,40);

siftscale=[2 2];
siftdescriptordensity=1;
minimagesize=floor(sqrt(2)*15*siftscale(2)+1);
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

featuretype=2;
timescale=1;
applyzscore=false;

featuretype=1;
timescale=4;
applyzscore=true;

clear rcounter;
clear routput;
clear rmean;
Fs=250;



for subject=23:23
    
EEG = prepareEEG(Fs,windowsize,downsize,120,23:23,1:8);
Fs=ceil(Fs/downsize);


for subject=23:23
    for trial=1:35
        for i=1:12 rcounter{subject}{trial}{i} = 0; end
        for flash=1:120
            rcounter{subject}{trial}{EEG(subject,trial,flash).stim} = rcounter{subject}{trial}{EEG(subject,trial,flash).stim}+1;
        end
        % Check if all the epochs contain 10 repetitions.
        for i=1:12
            %assert( rcounter{subject}{trial}{i} == 10 );
        end
        rcounter{subject}{trial}
    end
end


%%
% Build routput pasting epochs toghether...
clear hit
for subject=23:23
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
for subject=23:23
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
for subject=23:23
    for trial=1:35
        
        for i=1:12

            rmean{i} = routput{subject}{trial}{i};
            
            for c=channelRange
                rsignal{i}(:,c) = resample(rmean{i}(:,c),size(rmean{i},1)*timescale,size(rmean{i},1));
                %rsignal{i}(:,c) = resample(rmean{i}(:,c),1:size(rmean{i},1),timescale);
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
    for subject=23:23
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
    for subject=23:23
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
