globalfullspeller = zeros(20,8);

for pKS = 15:35
    run('ProcessP300.m');
    run('GeneralClassifyP300.m');
    globalfullspeller(pKS,:) = globalspeller(:);
end

[i,v] = max(globalfullspeller)