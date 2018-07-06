% Frequency components analysis    
clear all
close all

downsize=40;

load('p300.mat');

%dataX = notchsignal(data.X, channelRange,Fs);
datatrial = data.trial;

figure;
drawfft(data.X(:,7)',true,Fs);
axis([0 60 0 0.6]);

dataX = data.X;
%dataX = decimateaveraging(dataX,channelRange,downsize);

dataX = notchsignal(data.X, channelRange,Fs);

figure;
drawfft(dataX(:,7)',true,Fs);
axis([0 60 0 0.6]);

dataX = bandpasseeg(dataX, channelRange,Fs,4);

figure;
drawfft(dataX(:,7)',true,Fs);
axis([0 60 0 0.6]);

dataX = decimatesignal(dataX,channelRange,downsize); 
%dataX = decimateaveraging(dataX,channelRange,downsize);


figure;
drawfft(dataX(:,7)',true,Fs/downsize);
axis([0 60 0 0.6]); 
%axis([0 12 -0.05 0.6]);
