function PasteDataBlocks

dataparse(13,'p300-train-[2017.12.05-13.39.20].ov',20,1)
dataparse(13,'p300-train-[2017.12.05-13.56.55].ov',15,2);
PasteData(13);

end

function dataparse(subject,filename,Trials,orderno)

    addpath('./ov2mat/');

    %Parameters
    Fs=250;

    %subject = 2;
    convert_ov2mat(sprintf('./signals/Subject%d/%s',subject,filename),sprintf('./signals/p300-temp-%d.mat',subject));
    load(sprintf('./signals/p300-temp-%d.mat',subject));

    % NN.NNNNN
    % data.X(sample, channel)
    % data.y(sample)  --> 0: no, 1:nohit, 2:hit
    % data.y_stim(sample) --> 1-12, 1-6 cols, 7-12 rows

    channels={ 'Fz'  ,  'Cz',    'P3' ,   'Pz'  ,  'P4'  , 'PO7'   , 'PO8',  'Oz'};
    channelRange=1:8;


    samples(find(sampleTime==stims(1,1)),:)

    %%
    fprintf('%04x\n',stims(:,2))

    % First, get the stimulus to different events
    c=find(stims(:,2)==hex2dec('00008005')); % 32773 Trial Start
    d=find(stims(:,2)==hex2dec('00008006')); % 32774 Trial Stop
    a=find(stims(:,2)==hex2dec('00008205')); % 33285 Hit
    b=find(stims(:,2)==hex2dec('00008206')); % 32286 Nohit
    e=find(stims(:,2)==hex2dec('0000800C')); % Visual Stimulus Stop 32780
    f=find(stims(:,2)==hex2dec('0000800B')); % Visual Stimulus Start 32779

    % Find all the stimulus assosiated with row/col flashing.
    stimuls = [];
    for i=1:12
        stimuls = [stimuls; find(stims(:,2)==33025-1+i)];
    end

    %%
    % Chequear si la cantidad de estimulos encontradas coincide.
    % 33025 es el label 1.
    total=0;
    for i=1:12
        [i size(find(stims(:,2)==33025-1+i),1)]
        total=total+size(find(stims(:,2)==33025-1+i),1);
    end

    assert ((size(stimuls,1) == total), 'Stimulus found do not match.');

    %%
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

    for trial=1:Trials
        h=[];
        for i=1:20  % Hay 20 hits en cada trial, 20*35 = 700 que es el tamaño de a
            vl=stims(a((trial-1)*20+i)+1,1:2);
            %[(trial-1)*Trials+i vl(2)]
            h=[h vl(2)];
        end
        h = unique(h);
        [trial h]
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
    if (size(c,1)>Trials)

        stoptime=stims(c(Trials+1),1);

        stopsample=find(sampleTime>stims(c(Trials+1),1));


        sampleTime(stopsample(1):end,:) = [];
        samples(stopsample(1):end,:) = [];

        z(Trials+1) = [];

        targets(Trials*12*10+1:end,:) = [];
        stimulations(Trials*12*10+1:end,:) = [];
    end


    % Check target consistency
    Word = [];
    for trial=1:Trials
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

    fprintf('Copy-spelled word:');
    Word

    % Esto significat que stimulations coincide con los targets.

    %% Data Structure
    data = cell(0);

    data.X = samples;
    data.y = zeros(size(samples,1),1);
    data.y_stim = zeros(size(samples,1),1);
    data.trial=zeros(5,1);


    data.flash = [];
    durs=[];
    for i=1:size(targets,1)
        % Obtengo la localizacion donde esta el marcador de fin del estimulo
        % (e) y del principio (f)
        startset = find(stims(f,1)<=targets(i,1));
        startset = sort(startset);

        endset = find(stims(e,1)>=targets(i,1));
        endset = sort(endset);
        endd = endset(1); % Location on e.

        if (size(startset,1)==0)
            % En algunos casos la primera estimulacion cae despues del primer
            % target.
            start = endd;
        else
            start = startset(end);
        end

        duration = stims(e(endd),1)-stims(f(endd),1);


        if (duration == 0)
            duration = 1/Fs;
        end

        assert( duration > 0, 'Flash duration cannot be zero.');


        durs(end+1) = ceil(Fs*duration);
        % Marco donde inicia el flash y la duracion en sample points.

        idxset=find(sampleTime>=stims(f(start),1));
        idxset = sort(idxset);
        idx=idxset(1);

        data.flash(end+1,1) = idx
        data.flash(end,2) = ceil(Fs*duration);

        %fakeEEG=fakeeegoutput(4,targets(i,2),channelRange,25,100,4);

        % Marco todos los estimulos y targets donde el flash estuvo presente.
        for j=1:ceil(Fs*duration)
            data.y(idx+j-1) = targets(i,2);
            data.y_stim(idx+j-1) = stimulations(i,2);

            %fakeEEG(j,:);

        end


        data.flash(end,3) = stimulations(i,2);
        data.flash(end,4) = targets(i,2);


        if (targets(i,2)==2)
            %data.X(maximalsampleidx+ceil(Fs*duration)-1:maximalsampleidx+ceil(Fs*duration)-1+ceil(Fs*0.33),:) = zeros(ceil(Fs*0.33)+1,size(data.X,2));
            %data.X(maximalsampleidx-1:maximalsampleidx-1+ceil(Fs*1),:) = zeros(ceil(Fs*1)+1,size(data.X,2));

            %data.X(maximalsampleidx-1+ceil(Fs/2*1),:) = 1000*ones(1,size(data.X,2));

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
    for i=1:Trials*12*10
        ss=data.y_stim(data.flash(i)-5:data.flash(i)+40)'

        if (subject ~= 8 && i ~= 2526) && ...
                (subject ~= 11 && i ~= 2470) && ...
                (subject ~= 13 && i ~= 1974) && ...
                (subject ~= 20 && i ~= 1524)
            assert ( ss(5) == 0, 'Not zero');
        end


    end
    %%
    % LISTOOOOOO

    save(sprintf('datum%d.mat',orderno),'data');
end

%datum.data
%data

function PasteData(subject)

    datum1=load('datum1.mat');
    datum2=load('datum2.mat');


    newdata.X = [datum1.data.X;datum2.data.X];
    newdata.y = [datum1.data.y;datum2.data.y];
    newdata.y_stim = [datum1.data.y_stim;datum2.data.y_stim];
    % Update positions for the second part.
    newdata.trial = datum2.data.trial;
    newdata.trial = newdata.trial + size(datum1.data.X,1);
    newdata.flash = datum2.data.flash;
    newdata.flash(:,1) = newdata.flash(:,1) + size(datum1.data.X,1);
    newdata.trial = [datum1.data.trial datum2.data.trial];
    newdata.flash = [datum1.data.flash;datum2.data.flash];

    % Clean output file
    if (exist(sprintf('./signals/p300-subject-%02d.mat',subject),'file'))
        delete(sprintf('./signals/p300-subject-%02d.mat',subject));
    end
    
    data=newdata;
    
    clear newdata;

    save(sprintf('./signals/p300-subject-%02d.mat',subject));
end