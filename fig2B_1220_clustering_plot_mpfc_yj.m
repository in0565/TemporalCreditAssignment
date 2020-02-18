close all; clear all;

%%
cd('E:\Data\CRX\MeanFR\reward01');
ACC_type=load(['celltype_ACC_YJ_1220.mat'],'type');
ACC_PYR = ACC_type.type{1}(:,2:3);
ACC_INT = ACC_type.type{2}(:,2:3);

PRL_type=load(['celltype_PRL_YJ_1220.mat'],'type');
PRL_PYR = PRL_type.type{1}(:,2:3);
PRL_INT = PRL_type.type{2}(:,2:3);

IL_type=load(['celltype_IL_YJ_1220.mat'],'type');
IL_PYR = IL_type.type{1}(:,2:3);
IL_INT = IL_type.type{2}(:,2:3);

celldata=[ACC_PYR;ACC_INT;PRL_PYR;PRL_INT;IL_PYR;IL_INT];
celldata=cell2mat(celldata);

%% main
figure('PaperUnits','Centimeters','PaperPosition',[0 0 8.88 6.66]);
axes('Position',[0.1 0.4 0.5 0.5]); hold on;
  set(gca, 'FontSize', 11);
h1=scatter(celldata(:,1), celldata(:,2),10);
set(h1,'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor',[0 0 0], 'LineWidth', 1.25);
set(gcf,'color', 'w');
set(gca, 'XAxisLocation','bottom');
set(gca, 'XAxisLocation','top');

set(gca,'XLim',[0.1 0.5],'XTick',[[],0.2:0.2:0.6], 'FontSize', 11)
set(gca,'YLim',[0 45],'YTick',[0:10:40], 'FontSize', 11)
xlabel('Spike width (ms)', 'fontsize',12);
ylabel('Firing rate (Hz)', 'fontsize',12);

line([0.24 0.24],[0, 8.83],'LineWidth',1.6, 'LineStyle','--', 'Color',[.4 .4 .4]); 
line([0.24 0.8],[8.83, 8.83],'LineWidth',1.6, 'LineStyle','--', 'Color',[.4 .4 .4]); 
% Lower
axes('Position',[0.1 0.29 0.5 0.09]); hold on;
[nn, bins1]=hist(celldata(:,1),100);
bar(bins1, nn,1, 'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);

set(gca, 'XAxisLocation','top');
set(gca, 'YAxisLocation','right');
set(gca, 'YDir', 'reverse');
set(gca,'XLim',[0.1 0.5],'XTick',[], 'FontSize', 11)
set(gca,'YLim',[0 50],'YTick',[0 50], 'FontSize', 11)
line([0.24 0.24],[0, 50],'LineWidth',1.6, 'LineStyle','--', 'Color',[.4 .4 .4]);

% Right
axes('Position',[0.62 0.4 0.09 0.5]); hold on;
[nn, bins2]=hist(celldata(:,2),100);
barh(bins2, nn,1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
set(gca,'XLim',[0 150],'XTick',[[],150], 'FontSize', 11)
set(gca,'YLim',[0 20],'YTick',[], 'FontSize', 11)
line([0.24 350],[3.9244, 3.9244],'LineWidth',1.6, 'LineStyle','--', 'Color',[.4 .4 .4]);

%% save
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\fig2');
print(gcf,'-dtiff', ['fig2B_1220_cellcluster_dist']); 
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters',['fig2B_1220_cellcluster_dist' '.ai']);
close all;
