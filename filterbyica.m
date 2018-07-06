function dataX = filterbyica(dataX, zerocoeff)

if (nargin<2)
    zerocoeff=[1,2]
end

[output, A, W] = fastica(dataX');

icas = output';

channels = { '1','2','3','4','5','6','7','8'};
time1=275222.343262227;
time2= 279948.995187305;
channelRange=1:8;
plotthiseeg(icas,channels,channelRange,time1/1,time2/1,false,10);

WW = inv(W);
for i=zerocoeff
    icas(:,i) = zeros(size(icas,1),1);
end
eeg = WW * icas';
eeg=eeg'; 

dataX = eeg;


end