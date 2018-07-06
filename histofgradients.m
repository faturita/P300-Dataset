angles=[0 45 90 120 180];
values=randi(315,1,100)

edges = [0:45:315];


hold on
%[val, ord] = max(globalspellerrep(subject,:,10))
%y=globalspellerrep(subject,ord,:);
%y=reshape(y, [1 10])
%y=y*100;
%Xi = 0:0.1:size(y,2);
%Yi = pchip(1:size(y,2),y,Xi);
%plot(1:10,y,':','linestyle',linestyles{subject},'marker',mark{subject},...
%    'color',C(subject,:),'MarkerSize',3,'linewidth',2);
[n,x] = hist(values,edges)
h=bar(x,n,'BarWidth',0.3);
axis([-20 330 1 25]);
set(gca,'YTick',[0 30 70 90]);
set(gca,'XTick',edges);
h.FaceColor = [0 0.5 0.5];
set(0, 'DefaultAxesFontSize',16);

grid on
xlabel('Bin Angles');
ylabel('Number of weighted gradients')
%title({'K=3';'CSP'},'color',C(4,:),'fontsize',14,'fontweight','normal',...
    %'fontname','Courier');
%legend('1','2','3','4','5','6','7','8','location','best')