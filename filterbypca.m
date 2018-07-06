function dataX = filterbypca(dataX, zerocoeff)    

if (nargin<2)
    zerocoeff=[1,2]
end

[coeff, score, latent] = princomp(dataX);
    
channels = { '1','2','3','4','5','6','7','8'};
time1=275222.343262227;
time2= 279948.995187305;
channelRange=1:8;
plotthiseeg(score,channels,channelRange,time1/1,time2/1,false);

WW = inv(coeff);

for i=zerocoeff
    score(:,i) = zeros(size(score,1),1);
end

eeg = WW * score';
eeg=eeg'; 

dataX = eeg;    
    
    
end