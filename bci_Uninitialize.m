function [box_out]=bci_Uninitialize(box_in)

box_out = box_in;

global events;
global sync;

disp(' Unilitializing everything....');

file = strcat('C:\\Users\\User\\Desktop\\Emotiv\\events.dat');
dlmwrite ( file, events);

disp(sprintf('Saved to %s:',file));

s = sprintf ( '%18.10f', now);

fprintf('Value at %s : %d\n', s, sync);

end