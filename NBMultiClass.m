function [ACC, ERR, AUC, SC] = NBMultiClass(F,DE,channel,testRange,labelRange,distancetype,kparam)

fprintf('Channel %d\n', channel);
fprintf('Building Test Matrix M for Channel %d:', channel);
[TM, TIX] = BuildDescriptorMatrix(F,channel,labelRange,testRange);
fprintf('%d\n', size(TM,2));

%fprintf('Channel %d\n', channel);
%fprintf('Building Test Matrix M for Channel %d:', channel);
%[M, IX] = BuildDescriptorMatrix(F,channel,labelRange,trainingRange);
%fprintf('%d\n', size(TM,2));


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
K = size(DE.C(1).M,2);


% Obtengo las K distancias de los 30 a los 150.
[Z,I] = pdist2(DE.C(1).M',DE.C(2).M',distancetype,'Smallest',K );

k = K;

% Z es de 150 x 30.  En D guardo las sumas para cada descriptor hit de
% 1-30 (de las distancias a los K vecinos de cada uno de los descriptores
% nohits).   D me da entonces una medida de "lo diferente" de cada
% descriptor de la bolsa de descriptores.
D = sum(Z(1:k,1:size(DE.C(2).M,2)),1);

%ED = sum( Epanechnikov(D) );

%Wi = (1.-D/ED)/(30-1);
%Wi = Epanechnikov(D) / ED;


for f=1:size(testRange,2)/12


    K = size(DE.C(2).M,2);

    [Z,I] = pdist2(DE.C(2).M',(TM(:,mind:mind+12*9-1)'),distancetype,'Smallest',K );
    
    k = kparam;

    assert( k > 1, 'error');
    
    %sumsrow = dot(Z(1:k,1:6),Wi(I(1:k,1:6)));
    
    sumsrow = dot(Z(1:k,1:54),ones(k,54));
    sumsrow = reshape(sumsrow,[6 9]);
    sumsrow = sum(sumsrow,2);
    
    sumscol = dot(Z(1:k,55:108),ones(k,54));
    sumscol = reshape(sumscol,[6 9]);
    sumscol = sum(sumscol,2);

    % Me quedo con aquel que la suma contra todos, dio menor.
    [c, row] = min(sumsrow);
    [c, col] = min(sumscol);
    %col=col+6;

    % I(1:3,1:6) Me da en cada columna los ids de los descriptores de M mas
    % cercaos a cada uno de los descriptores de 1 a 6.

    %SC.CLSF{test}.predicted = DE.C(I(1)).Label;  
    %SC.CLSF{test}.IDX{clster} = IDX; 

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

    mind=mind+12*9;
    maxd=maxd+12*9;
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

end

function E = Epanechnikov(t)

E = 3/4 * (1 - (abs(t) <= 1).^2 );


end
