globalfullspeller = zeros(40,8);

for pZerolevel = -10:10
    run('ProcessP300.m');
    run('GeneralClassifyP300.m');
    globalfullspeller(pZerolevel+11,:) = globalspeller(:);
end

[i,v] = max(globalfullspeller)