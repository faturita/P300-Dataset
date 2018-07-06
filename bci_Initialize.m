function [box_out]=matlab_Initialize(box_in)

global  sync;
global  events;

%events = zeros([1000 2]);

sync = 0;

events = [ events [0; 0] ];

s = sprintf ( '%18.10f', now);

disp( sprintf('Value at %s : %d\n', s, sync) );

%out = myeeglogger(1,10,'RodrigoRamele','p300_init.dat');

%JUDP('SEND',7788,'127.0.0.1',int8('Start!'));

box_out = box_in;

end