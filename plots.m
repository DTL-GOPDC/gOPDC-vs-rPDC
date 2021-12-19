clc
clear
close all;
load('Connections.mat','rPDC_mean','gOPDC_mean')
list={'1-->1','1-->2','1-->3','2-->1','2-->2','2-->3','3-->1','3-->2','3-->3'}
gOPDC_mean=gOPDC_mean';
rPDC_mean=rPDC_mean';

connectivity=[gOPDC_mean(:) rPDC_mean(:)]

figure;
bar(connectivity)
set(gca,'XTickLabel',list);
legend('GOPDC','rPDC')

grid on;
ylabel('Normalized value of DC')
xlabel('Connection')
set(gca,'fontweight','bold','FontSize',18)