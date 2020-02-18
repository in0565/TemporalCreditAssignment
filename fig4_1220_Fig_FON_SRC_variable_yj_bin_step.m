clear all; close all;

%% load
cd('D:\KAIST\Grad\SNL\LDH_yj\mat\1220\Q')
file_list = ls('*.mat');
for ifile = 1:size(file_list,1)
    load(file_list(ifile,:))
end

load('E:\Data\Neural analysis\mean_memory_lengths_1220.mat');
mean_memory_length = nanmean(memory_length_1220);

ran_iter = 100;
%% FON
for icell=1:length(ACC_PYR)
    P_value_ACC{icell,1}=[P_value_ACC{icell,1}'; P_value_ACC{icell,2}'; P_value_ACC{icell,3}'; P_value_ACC{icell,4}'; P_value_ACC{icell,5}']';
    for iran = 1:ran_iter
        P_value_sh_ACC{icell,1,iran}=[P_value_sh_ACC{icell,1,iran}'; P_value_sh_ACC{icell,2,iran}'; P_value_sh_ACC{icell,3,iran}'; P_value_sh_ACC{icell,4,iran}'; P_value_sh_ACC{icell,5,iran}']';
    end
    for itype=1:6
        sig_ACC{itype,1}(icell,:)=P_value_ACC{icell}(itype+1,:);
        for iran = 1:ran_iter
            sig_sh_ACC{itype,1}(ran_iter*icell-ran_iter+iran,:)=P_value_sh_ACC{icell,1,iran}(itype+1,:);
        end
    end
end

for icell=1:length(PRL_PYR)
    P_value_PRL{icell,1}=[P_value_PRL{icell,1}'; P_value_PRL{icell,2}'; P_value_PRL{icell,3}'; P_value_PRL{icell,4}'; P_value_PRL{icell,5}']';
    for iran = 1:ran_iter
        P_value_sh_PRL{icell,1,iran}=[P_value_sh_PRL{icell,1,iran}'; P_value_sh_PRL{icell,2,iran}'; P_value_sh_PRL{icell,3,iran}'; P_value_sh_PRL{icell,4,iran}'; P_value_sh_PRL{icell,5,iran}']';
    end
    for itype=1:6
        sig_PRL{itype,1}(icell,:)=P_value_PRL{icell}(itype+1,:);
        for iran = 1:ran_iter
            sig_sh_PRL{itype,1}(ran_iter*icell-ran_iter+iran,:)=P_value_sh_PRL{icell,1,iran}(itype+1,:);
        end
    end
end

for icell=1:length(IL_PYR)
    P_value_IL{icell,1}=[P_value_IL{icell,1}'; P_value_IL{icell,2}'; P_value_IL{icell,3}'; P_value_IL{icell,4}'; P_value_IL{icell,5}']';
    for iran = 1:ran_iter
        P_value_sh_IL{icell,1,iran}=[P_value_sh_IL{icell,1,iran}'; P_value_sh_IL{icell,2,iran}'; P_value_sh_IL{icell,3,iran}'; P_value_sh_IL{icell,4,iran}'; P_value_sh_IL{icell,5,iran}']';
    end
    for itype=1:6
        sig_IL{itype,1}(icell,:)=P_value_IL{icell}(itype+1,:);
        for iran = 1:ran_iter
            sig_sh_IL{itype,1}(ran_iter*icell-ran_iter+iran,:)=P_value_sh_IL{icell,1,iran}(itype+1,:);
        end
    end
end

%% Finding significant cells (p<0.05)
for itype=1:3
    sig_ACC_05{itype,1}=sig_ACC{itype+3,1} < 0.05;
    sig_PRL_05{itype,1}=sig_PRL{itype+3,1} < 0.05;
    sig_IL_05{itype,1}=sig_IL{itype+3,1} < 0.05;

    sig_sh_ACC_05{itype,1}=sig_sh_ACC{itype+3,1} < 0.05;
    sig_sh_PRL_05{itype,1}=sig_sh_PRL{itype+3,1} < 0.05;
    sig_sh_IL_05{itype,1}=sig_sh_IL{itype+3,1} < 0.05;
    %% FON: Y shuffle (whole)
    FON_sh_ACC(itype,:)=sum(sig_sh_ACC_05{itype,1})/length(sig_sh_ACC_05{itype,1}(:,1))*100;
    FON_sh_PRL(itype,:)=sum(sig_sh_PRL_05{itype,1})/length(sig_sh_PRL_05{itype,1}(:,1))*100;
    FON_sh_IL(itype,:)=sum(sig_sh_IL_05{itype,1})/length(sig_sh_IL_05{itype,1}(:,1))*100;
    %% Fraction of Neurons (whole)
    FON_ACC(itype,:)=sum(sig_ACC_05{itype,1})/length(sig_ACC_05{itype,1}(:,1))*100;
    FON_PRL(itype,:)=sum(sig_PRL_05{itype,1})/length(sig_PRL_05{itype,1}(:,1))*100;
    FON_IL(itype,:)=sum(sig_IL_05{itype,1})/length(sig_IL_05{itype,1}(:,1))*100;
end

FON_n = {sig_ACC_05, sig_PRL_05, sig_IL_05};
FON_n_shuffle = {sig_sh_ACC_05, sig_sh_PRL_05, sig_sh_IL_05};
FON_results = {FON_ACC, FON_PRL, FON_IL};
FON_sh_results = {FON_sh_ACC, FON_sh_PRL, FON_sh_IL};

%% chance level: binomial test
bino_chance_p_t = binoinv([0.05 0.95], length(PRL_PYR), 0.05);
bino_chance_PRL = bino_chance_p_t(2)/length(PRL_PYR);
bino_chance_a_t = binoinv([0.05 0.95], length(ACC_PYR), 0.05);
bino_chance_ACC = bino_chance_a_t(2)/length(ACC_PYR);
bino_chance_i_t = binoinv([0.05 0.95], length(IL_PYR), 0.05);
bino_chance_IL = bino_chance_i_t(2)/length(IL_PYR);
N_nocut = [length(ACC_PYR) length(PRL_PYR) length(IL_PYR)];
plot_bino_chance = (bino_chance_PRL+bino_chance_ACC+bino_chance_IL)/3;

%% SRC
%% Rearranging data
%%% ACC
b_ACC = cell(6,1);
b_sh_ACC = cell(6,1);
mean_cell_sh_ACC = cell(6,1);
b_reward_ACC = zeros(length(ACC_PYR),3); b_reward_sh_ACC = zeros(length(ACC_PYR),3);
for icell=1:length(ACC_PYR)
    SRC_ACC{icell,1}=[SRC_ACC{icell,1}'; SRC_ACC{icell,2}'; SRC_ACC{icell,3}'; SRC_ACC{icell,4}'; SRC_ACC{icell,5}']';
    for iran = 1:ran_iter
        SRC_sh_ACC{icell,1,iran}=[SRC_sh_ACC{icell,1,iran}'; SRC_sh_ACC{icell,2,iran}'; SRC_sh_ACC{icell,3,iran}'; SRC_sh_ACC{icell,4,iran}'; SRC_sh_ACC{icell,5,iran}']';
    end
    for itype=1:6
        b_ACC{itype,1}(icell,:)=SRC_ACC{icell}(itype,:);
        for iran = 1:ran_iter
            b_sh_ACC{itype,1}(ran_iter*icell-ran_iter+iran,:)=SRC_sh_ACC{icell,1,iran}(itype,:);
        end
        mean_cell_sh_ACC{itype,1}(icell,:) = nanmean(abs(b_sh_ACC{itype,1}(ran_iter*icell-ran_iter+1:ran_iter*icell,:)));
    end
end

%%% PRL
b_PRL = cell(6,1);
b_sh_PRL = cell(6,1);
mean_cell_sh_PRL = cell(6,1);
b_reward_PRL = zeros(length(PRL_PYR),3); b_reward_sh_PRL = zeros(length(PRL_PYR),3);
for icell=1:length(PRL_PYR)
    SRC_PRL{icell,1}=[SRC_PRL{icell,1}'; SRC_PRL{icell,2}'; SRC_PRL{icell,3}'; SRC_PRL{icell,4}'; SRC_PRL{icell,5}']';
    for iran = 1:ran_iter
        SRC_sh_PRL{icell,1,iran}=[SRC_sh_PRL{icell,1,iran}'; SRC_sh_PRL{icell,2,iran}'; SRC_sh_PRL{icell,3,iran}'; SRC_sh_PRL{icell,4,iran}'; SRC_sh_PRL{icell,5,iran}']';
    end
    for itype=1:6
        b_PRL{itype,1}(icell,:)=SRC_PRL{icell}(itype,:);
        for iran = 1:ran_iter
            b_sh_PRL{itype,1}(ran_iter*icell-ran_iter+iran,:)=SRC_sh_PRL{icell,1,iran}(itype,:);
        end
        mean_cell_sh_PRL{itype,1}(icell,:) = nanmean(abs(b_sh_PRL{itype,1}(ran_iter*icell-ran_iter+1:ran_iter*icell,:)));
    end
end

%%% IL
b_IL = cell(6,1);
b_sh_IL = cell(6,1);
mean_cell_sh_IL = cell(6,1);
b_reward_IL = zeros(length(IL_PYR),3); b_reward_sh_IL = zeros(length(IL_PYR),3);
for icell=1:length(IL_PYR)
    SRC_IL{icell,1}=[SRC_IL{icell,1}'; SRC_IL{icell,2}'; SRC_IL{icell,3}'; SRC_IL{icell,4}'; SRC_IL{icell,5}']';
    for iran = 1:ran_iter
        SRC_sh_IL{icell,1,iran}=[SRC_sh_IL{icell,1,iran}'; SRC_sh_IL{icell,2,iran}'; SRC_sh_IL{icell,3,iran}'; SRC_sh_IL{icell,4,iran}'; SRC_sh_IL{icell,5,iran}']';
    end
    for itype=1:6
        b_IL{itype,1}(icell,:)=SRC_IL{icell}(itype,:);
        for iran = 1:ran_iter
            b_sh_IL{itype,1}(ran_iter*icell-ran_iter+iran,:)=SRC_sh_IL{icell,1,iran}(itype,:);
        end
        mean_cell_sh_IL{itype,1}(icell,:) = nanmean(abs(b_sh_IL{itype,1}(ran_iter*icell-ran_iter+1:ran_iter*icell,:)));
    end
end

%% mean and SEM to plot
% for sh: 1st mean,std inside neuron and then among neurons
for itype=1:6
    %% no Hz cut
    mean_b_ACC{itype} = nanmean(abs(b_ACC{itype,1}));
    SEM_b_ACC{itype} = nanstd(abs(b_ACC{itype,1}))./sqrt(sum(~isnan(b_ACC{itype,1})));
    mean_b_sh_ACC{itype} = nanmean(mean_cell_sh_ACC{itype,1});
    SEM_b_sh_ACC{itype} = nanstd(mean_cell_sh_ACC{itype,1})./sqrt(sum(~isnan(mean_cell_sh_ACC{itype,1})));
    
    mean_b_PRL{itype} = nanmean(abs(b_PRL{itype,1}));
    SEM_b_PRL{itype} = nanstd(abs(b_PRL{itype,1}))./sqrt(sum(~isnan(b_PRL{itype,1})));
    mean_b_sh_PRL{itype} = nanmean(mean_cell_sh_PRL{itype,1});
    SEM_b_sh_PRL{itype} = nanstd(mean_cell_sh_PRL{itype,1})./sqrt(sum(~isnan(mean_cell_sh_PRL{itype,1})));
    
    mean_b_IL{itype} = nanmean(abs(b_IL{itype,1}));
    SEM_b_IL{itype} = nanstd(abs(b_IL{itype,1}))./sqrt(sum(~isnan(b_IL{itype,1})));
    mean_b_sh_IL{itype} = nanmean(mean_cell_sh_IL{itype,1});
    SEM_b_sh_IL{itype} = nanstd(mean_cell_sh_IL{itype,1})./sqrt(sum(~isnan(mean_cell_sh_IL{itype,1})));
end

b_means = {mean_b_ACC, mean_b_PRL, mean_b_IL};
b_SEM = {SEM_b_ACC, SEM_b_PRL, SEM_b_IL};
b_sh_means = {mean_b_sh_ACC, mean_b_sh_PRL, mean_b_sh_IL};
b_sh_SEM = {SEM_b_sh_ACC, SEM_b_sh_PRL, SEM_b_sh_IL};
b_results = {(b_ACC{4:6,:}), (b_PRL{4:6,:}), (b_IL{4:6,:})}; % no abs
b_sh_results = {mean_cell_sh_ACC{4:6,:}, mean_cell_sh_PRL{4:6,:}, mean_cell_sh_IL{4:6,:}};

%% draw settings
x_win = [-1 3;-1 1;-1 3;-2 2;-2 3];  sliding2 = 0.05; bin2 = 0.5;
y_win=[1 4/sliding2+1;4/sliding2+2 6/sliding2+2;6/sliding2+3 10/sliding2+3;10/sliding2+4 14/sliding2+4;14/sliding2+5 19/sliding2+5]; % 500ms with 100ms
x_lim = [-1 3;-1 1;-1 3;-1 2;-2 2];
typetitle = {'Delay','Approach','Choice','Memory','Reward'};
regiontitle = {'ACC','PRL','ILC'};
maintitle = {'C','R','X'; 'Q_{contra}','Q_{ipsi}','Q_C'};
c_map = [0 0 255; 255 0 0; 0 204 0]/255;
c_map_sem = [153 153 255; 255 102 102; 153 255 153]/255;
y_sh_color = [64 64 191; 191 64 64; 179 230 179]/255;
po_fon = 30; po_src = 0.17;% po_cpd = 2;
data_linewidth = 0.5;
data_markersize = 7;
sh_linewidth = 0.35;
sh_markersize = 2.5;
%% draw FON
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 9]);
i=5;
for iregion=1:3
    % QL
    subplot(3,3,iregion)
    hold on
    line([-0.5 -0.5],[0, 99],'LineStyle','--', 'Color',[.6 .6 .6]);
    line([0 0],[0, 99],'LineStyle','-', 'Color',[.6 .6 .6]);
    line([-mean_memory_length/1000000 -mean_memory_length/1000000],[0 99],'LineStyle','--','Color',[120 177 255]/255);
    % shuffled(y)
    plot(x_win(i,1):sliding2:x_win(i,2), FON_sh_results{iregion}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','s', 'MarkerSize',sh_markersize,'MarkerFaceColor','none','LineWidth',sh_linewidth,'Color',y_sh_color(iregion,:));
    %%% main plot
    L2=plot(x_win(i,1):sliding2:x_win(i,2), FON_results{iregion}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','.', 'MarkerSize',data_markersize,'MarkerFaceColor',c_map(iregion,:),'LineWidth',data_linewidth,'Color',c_map(iregion,:));
    t=title(regiontitle{iregion});
    set(t, 'FontSize', 9);
    set(gca,'XLim',[x_lim(i,1) x_lim(i,2)],'XTick',[[], x_lim(i,1)+1:1:x_lim(i,2)-1,[]],'XColor','w');
    set(gca,'YLim',[0 po_fon],'YTick',[0:po_fon/2:po_fon],'YTickLabel',[0:po_fon/2:po_fon],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    if iregion==1
        text(-1.8, po_fon*0.8, maintitle(2,1), 'FontSize', 10);
    end
    
    % QR
    subplot(3,3,iregion+3)
    hold on
    line([-0.5 -0.5],[0, 99],'LineStyle','--', 'Color',[.6 .6 .6]);
    line([0 0],[0, 99],'LineStyle','-', 'Color',[.6 .6 .6]);
    line([-mean_memory_length/1000000 -mean_memory_length/1000000],[0 99],'LineStyle','--','Color',[120 177 255]/255);
    % shuffled(y)
    plot(x_win(i,1):sliding2:x_win(i,2), FON_sh_results{iregion}(2,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','s', 'MarkerSize',sh_markersize ,'MarkerFaceColor','none','LineWidth',sh_linewidth ,'Color',y_sh_color(iregion,:));
    %%% main plot
    plot(x_win(i,1):sliding2:x_win(i,2), FON_results{iregion}(2,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','.', 'MarkerSize',data_markersize ,'MarkerFaceColor',c_map(iregion,:),'LineWidth',data_linewidth,'Color',c_map(iregion,:));
    set(gca,'XLim',[x_lim(i,1) x_lim(i,2)],'XTick',[[], x_lim(i,1)+1:1:x_lim(i,2)-1,[]],'XColor','w');
    set(gca,'YLim',[0 po_fon],'YTick',[0:po_fon/2:po_fon],'YTickLabel',[0:po_fon/2:po_fon],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    if iregion==1
        y = ylabel(['Fraction of neurons(%)']); set(y,'FontSize',11);
        text(-1.8, po_fon*0.8, maintitle(2,2), 'FontSize', 10);
    end
    
    % QC
    subplot(3,3,iregion+6)
    hold on
    line([-0.5 -0.5],[0, 99],'LineStyle','--', 'Color',[.6 .6 .6]);
    line([0 0],[0, 99],'LineStyle','-', 'Color',[.6 .6 .6]);
    line([-mean_memory_length/1000000 -mean_memory_length/1000000],[0 99],'LineStyle','--','Color',[120 177 255]/255);
    % shuffled(y)
    plot(x_win(i,1):sliding2:x_win(i,2), FON_sh_results{iregion}(3,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','s', 'MarkerSize',sh_markersize ,'MarkerFaceColor','none','LineWidth',sh_linewidth ,'Color',y_sh_color(iregion,:));
    %%% main plot
    plot(x_win(i,1):sliding2:x_win(i,2), FON_results{iregion}(3,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','.', 'MarkerSize',data_markersize ,'MarkerFaceColor',c_map(iregion,:),'LineWidth',data_linewidth,'Color',c_map(iregion,:));
    set(gca,'XLim',[x_lim(i,1) x_lim(i,2)],'XTick',[[], x_lim(i,1)+1:1:x_lim(i,2)-1,[]],'XColor','w');
    set(gca,'YLim',[0 po_fon],'YTick',[0:po_fon/2:po_fon],'YTickLabel',[0:po_fon/2:po_fon],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    if iregion == 2
        x=xlabel(['Time from reward onset (s)']); set(x, 'FontSize', 11);
    end
    if iregion==1
%         ylabel(['Fraction of neurons(%)']);
        text(-1.8, po_fon*0.8, maintitle(2,3), 'FontSize', 10);
    end
end

%% save FON
% cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\fig5')
% print(gcf, '-dtiff', ['fig5_1217_V_FON_variable_sliding_' num2str(sliding2) '_win_' num2str(bin2) '_newev_newq_sh_y_' num2str(ran_iter) '.tif'])
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters',['fig5_1220_V_FON_variable_sliding_' num2str(sliding2) '_win_' num2str(bin2) '_newev_newq_sh_y_' num2str(ran_iter) '.ai']);
% close all;

%% draw SRC
f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 9]);
i=5;
for iregion=1:3
    % QL
    subplot(3,3,iregion)
    hold on
    line([-0.5 -0.5],[0, 99],'LineStyle','--', 'Color',[.6 .6 .6]);
    line([0 0],[0, 99],'LineStyle','-', 'Color',[.6 .6 .6]);
    line([-mean_memory_length/1000000 -mean_memory_length/1000000],[0 99],'LineStyle','--','Color',[120 177 255]/255);
    % shading
    fill_x= [x_win(i,1):sliding2:x_win(i,2), fliplr(x_win(i,1):sliding2:x_win(i,2))];
    % shading (y shuffle)
    fill_a_y1 = [(b_sh_means{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))-b_sh_SEM{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))),...
        fliplr((b_sh_means{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))+b_sh_SEM{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))))];
    fill_a_y = fill(fill_x, fill_a_y1, [0.8 0.8 0.8]); set(fill_a_y,'EdgeColor','none'); alpha(fill_a_y,0.5)
    % shuffled(y)
    plot(x_win(i,1):sliding2:x_win(i,2), b_sh_means{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','s', 'MarkerSize',sh_markersize ,'MarkerFaceColor','none','LineWidth',sh_linewidth ,'Color',y_sh_color(iregion,:));
    % shading (main)
    fill_a1 = [(b_means{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))-b_SEM{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))),...
        fliplr((b_means{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))+b_SEM{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))))];
    fill_a = fill(fill_x, fill_a1, c_map_sem(iregion,:)); set(fill_a,'EdgeColor','none'); alpha(fill_a,0.5)
    %%% main plot
    L2=plot(x_win(i,1):sliding2:x_win(i,2), b_means{iregion}{4}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','.', 'MarkerSize',data_markersize ,'MarkerFaceColor',c_map(iregion,:),'LineWidth',data_linewidth,'Color',c_map(iregion,:));
    t=title(regiontitle{iregion});
    set(t, 'FontSize', 9);
    set(gca,'XLim',[x_lim(i,1) x_lim(i,2)],'XTick',[[], x_lim(i,1)+1:1:x_lim(i,2)-1,[]],'XColor','w');
    set(gca,'YLim',[0.045 po_src],'YTick',[0:0.05:1],'YTickLabel',[0:0.05:1],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    if iregion==1
        text(-1.8, po_src*0.8, maintitle(2,1), 'FontSize', 10);
    end
    
    % QR
    subplot(3,3,iregion+3)
    hold on
    line([-0.5 -0.5],[0, 99],'LineStyle','--', 'Color',[.6 .6 .6]);
    line([0 0],[0, 99],'LineStyle','-', 'Color',[.6 .6 .6]);
    line([-mean_memory_length/1000000 -mean_memory_length/1000000],[0 99],'LineStyle','--','Color',[120 177 255]/255);
    % shading (y shuffle)
    fill_a_y1 = [(b_sh_means{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))-b_sh_SEM{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))),...
        fliplr((b_sh_means{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))+b_sh_SEM{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))))];
    fill_a_y = fill(fill_x, fill_a_y1, [0.8 0.8 0.8]); set(fill_a_y,'EdgeColor','none'); alpha(fill_a_y,0.5)
    % shuffled(y)
    plot(x_win(i,1):sliding2:x_win(i,2), b_sh_means{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','s', 'MarkerSize',sh_markersize ,'MarkerFaceColor','none','LineWidth',sh_linewidth ,'Color',y_sh_color(iregion,:));
    % shading (main)
    fill_a1 = [(b_means{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))-b_SEM{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))),...
        fliplr((b_means{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))+b_SEM{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))))];
    fill_a = fill(fill_x, fill_a1, c_map_sem(iregion,:)); set(fill_a,'EdgeColor','none'); alpha(fill_a,0.5)
    %%% main plot
    plot(x_win(i,1):sliding2:x_win(i,2), b_means{iregion}{5}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','.', 'MarkerSize',data_markersize ,'MarkerFaceColor',c_map(iregion,:),'LineWidth',data_linewidth,'Color',c_map(iregion,:));
    set(gca,'XLim',[x_lim(i,1) x_lim(i,2)],'XTick',[[], x_lim(i,1)+1:1:x_lim(i,2)-1,[]],'XColor','w');
    set(gca,'YLim',[0.045 po_src],'YTick',[0:0.05:1],'YTickLabel',[0:0.05:1],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    if iregion==1
        y2 = ylabel(['|SRC|']); set(y2, 'FontSize',11);
        text(-1.8, po_src*0.8, maintitle(2,2), 'FontSize', 10);
    end
    
    % QC
    subplot(3,3,iregion+6)
    hold on
    line([-0.5 -0.5],[0, 99],'LineStyle','--', 'Color',[.6 .6 .6]);
    line([0 0],[0, 99],'LineStyle','-', 'Color',[.6 .6 .6]);
    line([-mean_memory_length/1000000 -mean_memory_length/1000000],[0 99],'LineStyle','--','Color',[120 177 255]/255);
    % shading (y shuffle)
    fill_a_y1 = [(b_sh_means{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))-b_sh_SEM{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))),...
        fliplr((b_sh_means{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))+b_sh_SEM{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))))];
    fill_a_y = fill(fill_x, fill_a_y1, [0.8 0.8 0.8]); set(fill_a_y,'EdgeColor','none'); alpha(fill_a_y,0.5)
    % shuffled(y)
    plot(x_win(i,1):sliding2:x_win(i,2), b_sh_means{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','s', 'MarkerSize',sh_markersize ,'MarkerFaceColor','none','LineWidth',sh_linewidth ,'Color',y_sh_color(iregion,:));
    % shading (main)
    fill_a1 = [(b_means{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))-b_SEM{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))),...
        fliplr((b_means{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))+b_SEM{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2))))];
    fill_a = fill(fill_x, fill_a1, c_map_sem(iregion,:)); set(fill_a,'EdgeColor','none'); alpha(fill_a,0.5)
    %%% main plot
    plot(x_win(i,1):sliding2:x_win(i,2), b_means{iregion}{6}(1,y_win(i,1)-y_win(3,2):y_win(i,2)-y_win(3,2)),...
        'Marker','.', 'MarkerSize',data_markersize ,'MarkerFaceColor',c_map(iregion,:),'LineWidth',data_linewidth,'Color',c_map(iregion,:));
    set(gca,'XLim',[x_lim(i,1) x_lim(i,2)],'XTick',[[], x_lim(i,1)+1:1:x_lim(i,2)-1,[]],'XColor','w');
    set(gca,'YLim',[0.045 po_src],'YTick',[0:0.05:1],'YTickLabel',[0:0.05:1],'YColor','k'); hold on;
    set(gca, 'FontSize', 9);
    if iregion == 2
        x2=xlabel(['Time from reward onset (s)']); set(x2, 'FontSize', 11);
    end
    if iregion==1
        text(-1.8, po_src*0.8, maintitle(2,3), 'FontSize', 10);
    end
end

%% save SRC
% cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\fig5')
% print(gcf, '-dtiff', ['fig5_1217_V_SRC_variable_sliding_' num2str(sliding2) '_win_' num2str(bin2) '_newev_newq_sh_y_' num2str(ran_iter) '.tif'])
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters',['fig5_1220_V_SRC_variable_sliding_' num2str(sliding2) '_win_' num2str(bin2) '_newev_newq_sh_y_' num2str(ran_iter) '.ai']);
close all;
