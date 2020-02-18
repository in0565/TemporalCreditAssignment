%%
clear all; close all;

% Pyr cell data
cd('E:\Data\CRX\MeanFR\reward01');
ACC_type = load('celltype_ACC_YJ_1220.mat');       ACC_PYR = ACC_type.type{1}(:,1);
PRL_type = load('celltype_PRL_YJ_1220.mat');       PRL_PYR = PRL_type.type{1}(:,1);
IL_type = load('celltype_IL_YJ_1220.mat');         IL_PYR = IL_type.type{1}(:,1);

clear ACC_type, clear PRL_type, clear IL_type,

cd('D:\KAIST\Grad\SNL\LDH_yj\mat\1220\CRX')
load('ACC_C_t_R_t_1220_idx_bin0.5_sliding0.05.mat','C_sig_tf_ACC');
load('PRL_C_t_R_t_1220_idx_bin0.5_sliding0.05.mat','C_sig_tf_PRL');
load('IL_C_t_R_t_1220_idx_bin0.5_sliding0.05.mat','C_sig_tf_IL');

% settings
saving_path = 'D:\KAIST\Grad\SNL\LDH_yj\Figures\Final';
raster_line_width = 1;
raster_marker_size = 0.1;
%% choice cell - memory period
win = [-3 3]*1000000;
xlims = [-1 2];     bin = 0.01; % unit: second
wintick=(xlims(1,1)-0.25):bin:(xlims(1,2)+0.25); % for fig4A
% ACC
spktime_ACC = cell(length(PRL_PYR),1);
spk_data_ACC = cell(length(PRL_PYR),1);
fns_neuron_ACC = cell(length(PRL_PYR),1);
spk_time_ACC = cell(length(PRL_PYR),1);
memory_length_ACC = zeros(length(PRL_PYR),1);

cell_id = 'iL_60';
icell = 60;
choice_cell = IL_PYR{icell};
[path, name]=fileparts(choice_cell);
cd(path);
load('Event_AC_YJ_1215.mat');
ntrial = size(Event_Time,2);
eventtime = Event_Time;
mean_memory_length = mean(eventtime(5,:)-eventtime(4,:));
memory_length_ACC(icell,1)=mean_memory_length;

load('beh_.mat')
index{1}=(beh_{1,2}(1:end-1,1)==-1);
index{2}=(beh_{1,2}(1:end-1,1)==1);

name=[name '.t'];
spk_data_ACC{icell} = get_Spiketime(name)*100; %s->ms
date=strsplit(path,'\');
date = date{end};
date = strsplit(date,'_');
date = date{1};
fns_neuron_ACC{icell} = [date(end-4:end) '_' name];

% PSTH
for itrial=1:ntrial
    ievent=1;
    [~, time_idx_ACC] = histc(spk_data_ACC{icell,1},eventtime(ievent+3,itrial)+win(ievent,:));
    if isempty(time_idx_ACC); continue; end;
    spk_time_ACC{icell,1}{itrial,ievent}=spk_data_ACC{icell,1}(logical(time_idx_ACC))-eventtime(ievent+3,itrial);
end
spktime_ACC{icell,1}= {spk_time_ACC{icell,1}(index{1}(:,1),:) spk_time_ACC{icell,1}(index{2}(:,1),:)}; %index{1}=Lt, index{2}=Rt

% Rearrange data
tmp_spk=spktime_ACC{icell};
n_trial_choice = [size(tmp_spk{1},1) size(tmp_spk{2},1)];
n_trial = sum(n_trial_choice);
yy = [0:n_trial-1; 1:n_trial; NaN(1,n_trial)];

for ichoice=1:2
    n_spikes_each = cellfun(@length,tmp_spk{ichoice}(:,1));
    n_spikes = sum(n_spikes_each);
    if n_spikes==0
        x_pts_ACC{icell}{ichoice} = [];
        y_pts_ACC{icell}{ichoice} = [];
        continue;
    end
    tmp_x = [[cell2mat(tmp_spk{ichoice}(:,1))'; cell2mat(tmp_spk{ichoice}(:,1))']; NaN(1,n_spikes)];
    tmp_x = tmp_x(:); tmp_x(end) = [];
    tmp_y = [];
    for iy=1:n_trial_choice(ichoice)
        if n_spikes_each(iy)==0; continue; end;
        if ichoice==1
            tmp_y = [tmp_y repmat(yy(:,iy),1,n_spikes_each(iy))];
        else
            tmp_y = [tmp_y repmat(yy(:,iy+sum(n_trial_choice(1:ichoice-1))),1,n_spikes_each(iy))];
        end
    end
    tmp_y = tmp_y(:); tmp_y(end) = [];
    x_pts_ACC{icell}{ichoice} = tmp_x/1000000;
    y_pts_ACC{icell}{ichoice} = tmp_y;
end

%% Raster PSTH그리기
lineclr={[0.0 0.0 0.4] [0.7 0.0 0.0]}; %choice
linewth=[0.8 0.8];

% Axes Point 생성
npp = 1; % numbers of cells per page
PSTH=cell(1,2);
f=figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 6]);
ylims = [1 n_trial];

% spike train
subplot(2,11,1:3)
set(gca, 'FontSize', 8);
hold on;

line([0 0],[0, 300],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([mean_memory_length/1000000 mean_memory_length/1000000], [0 300],'LineWidth',0.5,'LineStyle','-', 'Color',[1.0 0.6 0.4]);
t=title('C(t)');  %%%%%%%%%%%%%%%%%%% brain region 각각 지정하기. %%%%%%%%%
set(t, 'FontSize', 9,'FontAngle', 'italic');

for iichoice=1:2
    plot(x_pts_ACC{1,icell}{1,iichoice}, y_pts_ACC{1,icell}{1,iichoice},...
        'LineStyle','-','LineWidth',raster_line_width,'Color',lineclr{iichoice},'Marker','o','MarkerSize',raster_marker_size);
end
set(gca,'XLim',xlims(1,:),'XTick',[],'XColor','w');
set(gca,'YLim',ylims,'YTick',ylims,'YTickLabel',{[],ylims(2)},'YColor','k');
set(gca, 'FontSize', 9);
set(gca,'box','off','TickDir','out','FontSize',8);
L1=ylabel({'Trial','number'});set(L1, 'FontSize',9);

% sdf
for ichoice=1:2
    tmp_x=[x_pts_ACC{icell}{ichoice};NaN];
    if isnan(tmp_x); continue; end;
    tmp_x2=reshape(tmp_x,3,length(tmp_x)/3);
    spike_hist = histc(tmp_x2(1,:),wintick);
    PSTH{ichoice} = conv(spike_hist,fspecial('Gaussian',[1 100],10),'same')/(bin*n_trial_choice(ichoice));
end

% PSTH with SDF
subplot(2,11,12:14)
hold on;

lim=max([max(PSTH{1}), max(PSTH{2})]);

for ichoice=1:2
    plot(wintick,PSTH{ichoice},...
        'LineStyle','-','LineWidth',linewth(ichoice),'Color',lineclr{ichoice});
end

ylims2 = ceil(lim);

set(gca,'XLim',xlims(1,:), 'XTick',[[],xlims(1,1)+1:1:xlims(1,2)-1,[]]);
set(gca,'YLim',[0 ylims2],'YTick',[0 ylims2],'YTickLabel',[0,ylims2],'YColor','k'); %'
set(gca,'box','off','TickDir','out','FontSize',11);

line([0 0],[0, 99],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([mean_memory_length/1000000 mean_memory_length/1000000], [0 300],'LineWidth',0.5,'LineStyle','-', 'Color',[1.0 0.6 0.4]);

L1=ylabel({'Firing', 'rate (Hz)'}); set(L1, 'FontSize',9)

%% choice cell - reward period
win = [-3 3]*1000000;
xlims = [-2 2];     bin = 0.01; % unit: second
wintick=(xlims(1,1)-0.25):bin:(xlims(1,2)+0.25);
% ACC
spktime_ACC = cell(length(PRL_PYR),1);
spk_data_ACC = cell(length(PRL_PYR),1);
fns_neuron_ACC = cell(length(PRL_PYR),1);
spk_time_ACC = cell(length(PRL_PYR),1);
memory_length_ACC = zeros(length(PRL_PYR),1);

[path, name]=fileparts(choice_cell);
cd(path);
load('Event_AC_YJ_1215.mat');
ntrial = size(Event_Time,2);
eventtime = Event_Time;
mean_memory_length = mean(eventtime(5,:)-eventtime(4,:));
memory_length_ACC(icell,1)=mean_memory_length;

load('beh_.mat')
index{1}=(beh_{1,2}(1:end-1,1)==-1);
index{2}=(beh_{1,2}(1:end-1,1)==1);

name=[name '.t'];
spk_data_ACC{icell} = get_Spiketime(name)*100; %s->ms
date=strsplit(path,'\');
date = date{end};
date = strsplit(date,'_');
date = date{1};
fns_neuron_ACC{icell} = [date(end-4:end) '_' name];

% PSTH
for itrial=1:ntrial
    ievent=1;
    [~, time_idx_ACC] = histc(spk_data_ACC{icell,1},eventtime(ievent+4,itrial)+win(ievent,:));
    if isempty(time_idx_ACC); continue; end;
    spk_time_ACC{icell,1}{itrial,ievent}=spk_data_ACC{icell,1}(logical(time_idx_ACC))-eventtime(ievent+4,itrial);
end
spktime_ACC{icell,1}= {spk_time_ACC{icell,1}(index{1}(:,1),:) spk_time_ACC{icell,1}(index{2}(:,1),:)}; %index{1}=Lt, index{2}=Rt

% Rearrange data
tmp_spk=spktime_ACC{icell};
n_trial_choice = [size(tmp_spk{1},1) size(tmp_spk{2},1)];
n_trial = sum(n_trial_choice);
yy = [0:n_trial-1; 1:n_trial; NaN(1,n_trial)];

for ichoice=1:2
    n_spikes_each = cellfun(@length,tmp_spk{ichoice}(:,1));
    n_spikes = sum(n_spikes_each);
    if n_spikes==0
        x_pts_ACC{icell}{ichoice} = [];
        y_pts_ACC{icell}{ichoice} = [];
        continue;
    end
    tmp_x = [[cell2mat(tmp_spk{ichoice}(:,1))'; cell2mat(tmp_spk{ichoice}(:,1))']; NaN(1,n_spikes)];
    tmp_x = tmp_x(:); tmp_x(end) = [];
    tmp_y = [];
    for iy=1:n_trial_choice(ichoice)
        if n_spikes_each(iy)==0; continue; end;
        if ichoice==1
            tmp_y = [tmp_y repmat(yy(:,iy),1,n_spikes_each(iy))];
        else
            tmp_y = [tmp_y repmat(yy(:,iy+sum(n_trial_choice(1:ichoice-1))),1,n_spikes_each(iy))];
        end
    end
    tmp_y = tmp_y(:); tmp_y(end) = [];
    x_pts_ACC{icell}{ichoice} = tmp_x/1000000;
    y_pts_ACC{icell}{ichoice} = tmp_y;
end

%% Raster PSTH그리기 2
lineclr={[0.0 0.0 0.4] [0.7 0.0 0.0]}; %choice
linewth=[0.8 0.8];


PSTH=cell(1,2);
ylims = [1 n_trial];

% spike train
subplot(2,11,4:7)
set(gca, 'FontSize', 8);
hold on;

line([0 0],[0, 300],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-0.5 -0.5],[0, 300],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 300],'LineWidth',0.5,'LineStyle','-', 'Color',[1.0 0.6 0.4]);

for iichoice=1:2
    plot(x_pts_ACC{1,icell}{1,iichoice}, y_pts_ACC{1,icell}{1,iichoice},...
        'LineStyle','-','LineWidth',raster_line_width,'Color',lineclr{iichoice},'Marker','o','MarkerSize',raster_marker_size);
end
set(gca,'XLim',xlims(1,:),'XTick',[],'XColor','w');
set(gca,'YLim',ylims,'YTick',ylims,'YTickLabel',{},'YColor','k');
set(gca, 'FontSize', 9);
set(gca,'box','off','TickDir','out','FontSize',8);

% sdf
for ichoice=1:2
    tmp_x=[x_pts_ACC{icell}{ichoice};NaN];
    if isnan(tmp_x); continue; end;
    tmp_x2=reshape(tmp_x,3,length(tmp_x)/3);
    spike_hist = histc(tmp_x2(1,:),wintick);
    PSTH{ichoice} = conv(spike_hist,fspecial('Gaussian',[1 100],10),'same')/(bin*n_trial_choice(ichoice));
end


% PSTH with SDF
subplot(2,11,15:18)
hold on;

lim=max([max(PSTH{1}), max(PSTH{2})]);

for ichoice=1:2
    plot(wintick,PSTH{ichoice},...
        'LineStyle','-','LineWidth',linewth(ichoice),'Color',lineclr{ichoice});
end

ylims2 = ceil(lim);

set(gca,'XLim',xlims(1,:), 'XTick',[[],xlims(1,1)+1:1:xlims(1,2)-1,[]]);
set(gca,'YLim',[0 ylims2],'YTick',[0 ylims2],'YTickLabel',[],'YColor','k'); %'
set(gca,'box','off','TickDir','out','FontSize',11);

line([0 0],[0, 99],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-0.5 -0.5],[0, 99],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 300],'LineWidth',0.5,'LineStyle','-', 'Color',[1.0 0.6 0.4]);

%% reward cell - reward period
win = [-3 3]*1000000;
xlims = [-2 2];     bin = 0.01; % unit: second
wintick=(xlims(1,1)-0.25):bin:(xlims(1,2)+0.25);

% ACC
spktime_ACC = cell(length(ACC_PYR),1);
spk_data_ACC = cell(length(ACC_PYR),1);
fns_neuron_ACC = cell(length(ACC_PYR),1);
ttest_ACC = cell(length(ACC_PYR),1); plot_ACC = cell(length(ACC_PYR),1);
time_idx_ACC = [];   spk_time_ACC = cell(length(ACC_PYR),1);

icell = 53;    %1220

[path, name]=fileparts(ACC_PYR{icell});
cd(path);
load('Event_AC_YJ_1220.mat');
ntrial = size(Event_Time,2);
eventtime = Event_Time;
mean_memory_length = mean(eventtime(5,:)-eventtime(4,:));

load('beh_.mat')
index{1}=(beh_{1,2}(1:end-1,2)==-1); % reward
index{2}=(beh_{1,2}(1:end-1,2)==1);

name=[name '.t'];
spk_data_ACC{icell} = get_Spiketime(name)*100; %s->ms
date=strsplit(path,'\');
date = date{end};
date = strsplit(date,'_');
date = date{1};
fns_neuron_ACC{icell} = [date(end-4:end) '_' name];

% PSTH
for itrial=1:ntrial
    ievent=1;
    [~, time_idx_ACC] = histc(spk_data_ACC{icell,1},eventtime(ievent+4,itrial)+win(ievent,:));
    if isempty(time_idx_ACC); continue; end;
    spk_time_ACC{icell,1}{itrial,ievent}=spk_data_ACC{icell,1}(logical(time_idx_ACC))-eventtime(ievent+4,itrial);
end
spktime_ACC{icell,1}= {spk_time_ACC{icell,1}(index{1}(:,1),:) spk_time_ACC{icell,1}(index{2}(:,1),:)}; %index{1}=Lt, index{2}=Rt

% Rearrange data
tmp_spk=spktime_ACC{icell};
n_trial_choice = [size(tmp_spk{1},1) size(tmp_spk{2},1)];
n_trial = sum(n_trial_choice);
yy = [0:n_trial-1; 1:n_trial; NaN(1,n_trial)];

for ichoice=1:2
    n_spikes_each = cellfun(@length,tmp_spk{ichoice}(:,1));
    n_spikes = sum(n_spikes_each);
    if n_spikes==0
        x_pts_ACC{icell}{ichoice} = [];
        y_pts_ACC{icell}{ichoice} = [];
        continue;
    end
    tmp_x = [[cell2mat(tmp_spk{ichoice}(:,1))'; cell2mat(tmp_spk{ichoice}(:,1))']; NaN(1,n_spikes)];
    tmp_x = tmp_x(:); tmp_x(end) = [];
    tmp_y = [];
    for iy=1:n_trial_choice(ichoice)
        if n_spikes_each(iy)==0; continue; end;
        if ichoice==1
            tmp_y = [tmp_y repmat(yy(:,iy),1,n_spikes_each(iy))];
        else
            tmp_y = [tmp_y repmat(yy(:,iy+sum(n_trial_choice(1:ichoice-1))),1,n_spikes_each(iy))];
        end
    end
    tmp_y = tmp_y(:); tmp_y(end) = [];
    x_pts_ACC{icell}{ichoice} = tmp_x/1000000;
    y_pts_ACC{icell}{ichoice} = tmp_y;
end


%% Raster PSTH그리기 3
lineclr={[0.6 0.6 0.6] [0.0 0.4 0.8]}; %reward
linewth=[0.8 0.8];

% Axes Point 생성
npp = 1; % numbers of cells per page
PSTH=cell(1,2);
ylims = [1 n_trial];

% spike train
a1 = subplot(2,11,8:11);
set(gca, 'FontSize', 8);
set(a1,'YAxisLocation','Right')
hold on;

line([-0.5 -0.5],[0, 900],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([0 0],[0, 900],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 900],'LineWidth',0.5,'LineStyle','-', 'Color',[120 177 255]/255);
t=title('R(t)');  %%%%%%%%%%%%%%%%%%% brain region 각각 지정하기. %%%%%%%%%
set(t, 'FontSize', 9,'FontAngle', 'italic');

for iichoice=1:2
    plot(x_pts_ACC{1,icell}{1,iichoice}, y_pts_ACC{1,icell}{1,iichoice},...
        'LineStyle','-','LineWidth',raster_line_width,'Color',lineclr{iichoice},'Marker','o','MarkerSize',raster_marker_size);
end
set(gca,'XLim',xlims(1,:),'XTick',[],'XColor','w');
set(gca,'YLim',ylims,'YTick',ylims,'YTickLabel',{},'YColor','k');
set(gca,'box','off','TickDir','out','FontSize',8);

% sdf
for ichoice=1:2
    tmp_x=[x_pts_ACC{icell}{ichoice};NaN];
    if isnan(tmp_x); continue; end;
    tmp_x2=reshape(tmp_x,3,length(tmp_x)/3);
    spike_hist = histc(tmp_x2(1,:),wintick);
    PSTH{ichoice} = conv(spike_hist,fspecial('Gaussian',[1 100],10),'same')/(bin*n_trial_choice(ichoice));
end

% PSTH with SDF
a2 = subplot(2,11,19:22);
set(a2,'YAxisLocation','Right')
hold on;

lim=max([max(PSTH{1}), max(PSTH{2})]);

for ichoice=1:2
    plot(wintick,PSTH{ichoice},...
        'LineStyle','-','LineWidth',linewth(ichoice),'Color',lineclr{ichoice});
end

ylims2 = ceil(lim);

set(gca,'XLim',xlims(1,:), 'XTick',[[],xlims(1,1)+1:1:xlims(1,2)-1,[]]);
set(gca,'YLim',[0 ylims2],'YTick',[0 ylims2],'YTickLabel',[],'YColor','k'); %'
set(gca,'box','off','TickDir','out','FontSize',11);
line([-0.5 -0.5],[0, 900],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([0 0],[0, 900],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 900],'LineWidth',0.5,'LineStyle','-', 'Color',[120 177 255]/255);

%% save
cd(saving_path)
% print(gcf,'-dtiff',['choice_cell_' cell_id])
print(gcf, '-depsc','-painters',['choice_cell_' cell_id '.ai']);
close all
