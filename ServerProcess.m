disp('Waiting for plotting information on fftsignal variable.');
msg = '';

M = [];

% 7788 is the UDP port.
[MSSG,sourcehost] = judp('RECEIVE', 7788, 10024, 100000);
msg = char(MSSG');
while ( strcmpi(msg,'finishplot') == 0)
    
  
    disp(msg);
    
    %eval(msg);
    %figure(1);
    %plot(fftsignal);
    
    %M = [M ;[fftsignal]];
    

    [MSSG,sourcehost] = judp('RECEIVE', 7788, 10024, 100000);
    msg = char(MSSG');    
    disp(msg);

end