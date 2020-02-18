clear; close all;
%% 191217 discussion
%% load
cd('D:\KAIST\Grad\SNL\LDH_yj\mat\1220_reward01')
win_title = 1.0;
load(['totalCRX_ACC_yj_reward_' num2str(win_title) '.mat'])
load(['totalCRX_ACC_yj_y_100_shuffle_reward_' num2str(win_title) '.mat'])
load(['totalCRX_PRL_yj_reward_' num2str(win_title) '.mat'])
load(['totalCRX_PRL_yj_y_100_shuffle_reward_' num2str(win_title) '.mat'])
load(['totalCRX_IL_yj_reward_' num2str(win_title) '.mat'])
load(['totalCRX_IL_yj_y_100_shuffle_reward_' num2str(win_title) '.mat'])

%% color
c_map = [88, 166, 166; 66, 30, 34; 233, 128, 22]/255;

%% SRC %%%%%%%%%%%%%%%%%%%%%%
%% rearrange data
b_ACC = zeros(length(ACC_PYR),6);
b_sh_ACC = zeros(length(ACC_PYR)*100,6);
mean_cell_sh_ACC = zeros(length(ACC_PYR),6);
for icell = 1:length(ACC_PYR)
    b_ACC(icell,:) = SRC_ACC{icell};
    for iran = 1:100
        b_sh_ACC(100*icell-100+iran,:)=SRC_sh_ACC{icell,iran};
    end
    mean_cell_sh_ACC(icell,:) = nanmean(abs(b_sh_ACC(100*icell-99:100*icell,:)));
end

b_PRL = zeros(length(PRL_PYR),6);
b_sh_PRL = zeros(length(PRL_PYR)*100,6);
mean_cell_sh_PRL = zeros(length(PRL_PYR),6);
for icell = 1:length(PRL_PYR)
    b_PRL(icell,:) = SRC_PRL{icell};
    for iran = 1:100
        b_sh_PRL(100*icell-100+iran,:)=SRC_sh_PRL{icell,iran};
    end
    mean_cell_sh_PRL(icell,:) = nanmean(abs(b_sh_PRL(100*icell-99:100*icell,:)));
end

b_IL = zeros(length(IL_PYR),6);
b_sh_IL = zeros(length(IL_PYR)*100,6);
mean_cell_sh_IL = zeros(length(IL_PYR),6);
for icell = 1:length(IL_PYR)
    b_IL(icell,:) = SRC_IL{icell};
    for iran = 1:100
        b_sh_IL(100*icell-100+iran,:)=SRC_sh_IL{icell,iran};
    end
    mean_cell_sh_IL(icell,:) = nanmean(abs(b_sh_IL(100*icell-99:100*icell,:)));
end

%%

for itype = 1:3
    % no cut
    mean_b_ACC(itype) = nanmean(abs(b_ACC(:,itype+3)));
    SEM_b_ACC(itype) = nanstd(abs(b_ACC(:,itype+3)))./sqrt(sum(~isnan(b_ACC(:,itype+3))));
    mean_b_sh_ACC(itype) = nanmean(mean_cell_sh_ACC(:,itype+3));
    SEM_b_sh_ACC(itype) = nanstd(mean_cell_sh_ACC(:,itype+3))./sqrt(sum(~isnan(mean_cell_sh_ACC(:,itype+3))));
    
    mean_b_PRL(itype) = nanmean(abs(b_PRL(:,itype+3)));
    SEM_b_PRL(itype) = nanstd(abs(b_PRL(:,itype+3)))./sqrt(sum(~isnan(b_PRL(:,itype+3))));
    mean_b_sh_PRL(itype) = nanmean(mean_cell_sh_PRL(:,itype+3));
    SEM_b_sh_PRL(itype) = nanstd(mean_cell_sh_PRL(:,itype+3))./sqrt(sum(~isnan(mean_cell_sh_PRL(:,itype+3))));
    
    mean_b_IL(itype) = nanmean(abs(b_IL(:,itype+3)));
    SEM_b_IL(itype) = nanstd(abs(b_IL(:,itype+3)))./sqrt(sum(~isnan(b_IL(:,itype+3))));
    mean_b_sh_IL(itype) = nanmean(mean_cell_sh_IL(:,itype+3));
    SEM_b_sh_IL(itype) = nanstd(mean_cell_sh_IL(:,itype+3))./sqrt(sum(~isnan(mean_cell_sh_IL(:,itype+3))));
    
end

b_means = {mean_b_ACC, mean_b_PRL, mean_b_IL};
b_SEM = {SEM_b_ACC, SEM_b_PRL, SEM_b_IL};
b_sh_means = {mean_b_sh_ACC, mean_b_sh_PRL, mean_b_sh_IL};
b_sh_SEM = {SEM_b_sh_ACC, SEM_b_sh_PRL, SEM_b_sh_IL};
b_results = {abs(b_ACC(:,4:6)), abs(b_PRL(:,4:6)), abs(b_IL(:,4:6))};
b_sh_results = {mean_cell_sh_ACC(:,4:6), mean_cell_sh_PRL(:,4:6), mean_cell_sh_IL(:,4:6)};

%% draw
close all;
f3 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 9.5]);
figtitle = {'ACC','PRL','ILC'};
po = 0.17; %0.6

for iresult = 1:length(b_means)
    ifigtitle = figtitle{iresult};
    plot_src = b_means{iresult}; plot_sem = b_SEM{iresult};
    plot_src_sh = b_sh_means{iresult}; plot_sem_sh = b_sh_SEM{iresult};
    d_set = b_results{iresult};
    d_sh_set = b_sh_results{iresult};
    
    %% plot bar graph
    % QL vs QR vs QC
    subplot(length(b_means),1,iresult)
    hold on
    bar(1,plot_src(1),0.5,'FaceColor',c_map(1,:),'EdgeColor',c_map(1,:),'LineWidth',0.5)
    errorbar(1, plot_src(1), plot_sem(1),'Color','k','LineWidth',0.5)
    
    bar(2,plot_src(2),0.5,'FaceColor',c_map(2,:),'EdgeColor',c_map(2,:),'LineWidth',0.5)
    errorbar(2, plot_src(2), plot_sem(2),'Color','k','LineWidth',0.5)
    
    bar(3,plot_src(3),0.5,'FaceColor',c_map(3,:),'EdgeColor',c_map(3,:),'LineWidth',0.5)
    errorbar(3, plot_src(3), plot_sem(3),'Color','k','LineWidth',0.5)
    t=title(ifigtitle);
    set(t, 'FontSize', 10);
    % settings
    if iresult ==3
        xticks = {'' 'Q_{contra}' 'Q_{ipsi}' 'Q_C' ''};
    else
        xticks = {};
    end
    set(gca,'XLim',[0 4], 'XTick',[0:4],'XTickLabel',xticks,'XColor','w');
    set(gca,'YLim',[0 po],'YTick',[0:0.05:1],'YTickLabel',[0:0.05:1],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    y=ylabel(['|SRC|']); set(y,'FontSize',10);
    set(gca,'YGrid','on')
    % 1-ANOVA
    [~,~,st_a] = anova1([d_set(:,1); d_set(:,2); d_set(:,3)],...
        [ones(length(d_set(:,1)),1); ones(length(d_set(:,1)),1)*2; ones(length(d_set(:,1)),1)*3],'off');
    ca = multcompare(st_a,'display','off')
    plot_anova = ca(:,3).*ca(:,5)>0;

    % sig line
    src_line = 0;
    if plot_anova(1); src_line = 1; line([1 1.95], [po-0.025 po-0.025], 'LineStyle','-','Color','k','LineWidth',0.5); end
    if plot_anova(2); src_line = 1; line([1 3], [po-0.02 po-0.02], 'LineStyle','-','Color','k','LineWidth',0.5); end
    if plot_anova(3); src_line = 1; line([2.05 3], [po-0.025 po-0.025], 'LineStyle','-','Color','k','LineWidth',0.5); end
    if src_line==0; text(1,po*0.8,'n.s.','FontSize',12); end

end

%% save
% cd(['D:\KAIST\Grad\SNL\LDH_yj\Figures\fig5'])
% print(gcf, '-dtiff', ['fig5_1215_SRC_bar_V_reward_' num2str(win_title) '.tif']);
% close all
% cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
% print(gcf, '-depsc','-painters',['fig5_1220_SRC_bar_V_reward_' num2str(win_title) '.ai'])
