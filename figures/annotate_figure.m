%%figure
close all;

f = openfig('BER_multi.fig')
% 
% legend('DPC with equalization','HOGMT-Precoding: \epsilon {=} 10^{-2}','HOGMT-Precoding: \epsilon {=} 10^{-3}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-4}','Ideal(AWGN)','Location','Southwest','fontsize',16)
% a = get(gca,'Children');
% xdata = get(a, 'XData');
% ydata = get(a, 'YData');
% f.Children(2).Children(1).MarkerFaceColor = 'flat';
% f.Children(2).Children(2).MarkerFaceColor = 'flat';
% f.Children(2).Children(3).MarkerFaceColor = 'flat';
% f.Children(2).Children(4).MarkerFaceColor = 'flat';
% 
% f.Children(2).Children(1).SizeData = 100;
% f.Children(2).Children(2).SizeData = 100;
% f.Children(2).Children(3).SizeData = 100;
% f.Children(2).Children(4).SizeData = 100;
% grid on;
% box on;
set(gca, 'fontsize', 30)
% ylabel('Complexity')
% set(gca,'defaultMarkerSize',10);
% set(0,'defaultTextFontSize', 20)                      	% Default Font Size
% set(0,'defaultAxesFontSize', 20)  
% set(0,'defaultAxesFontSize',20)
% set(gca,'XTick',1:5:21);
% set(gca,'XTickLabel',(1:5:21));
% h = legend;
% set(h,'FontSize',20);
filename = 'BER_multi';
[h, wd, ht] = tightfig();
%%
save_file = 1;
if save_file == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end
% print -opengl -dpdf -r600 Ber_multi_mod.pdf


% legend('DPC','HOGMT-Precoding: \epsilon {=} 10^{-1}','HOGMT-Precoding: \epsilon {=} 5{\times}10^{-2}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-2}','Ideal(AWGN)','Location','Southwest','fontsize',16)