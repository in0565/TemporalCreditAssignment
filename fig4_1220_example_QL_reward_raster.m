clear all; close all;

%% Pyr cell data
cd('E:\Data\CRX\MeanFR\reward01');
PRL_type = load('celltype_PRL_YJ_1215.mat');       PRL_PYR = PRL_type.type{1}(:,1);
IL_type = load('celltype_IL_YJ_1215.mat');         IL_PYR = IL_type.type{1}(:,1);

clear PRL_type, clear IL_type,

%% settings
saving_path = 'D:\KAIST\Grad\SNL\LDH_yj\Figures\Final';
win = [-3 3]*1000000;
xlims = [-2 2];     bin = 0.01; % unit: second
wintick=(xlims(1,1)-0.25):bin:(xlims(1,2)+0.25);

%% PRL
spktime_PRL = cell(length(PRL_PYR),1);
spk_data_PRL = cell(length(PRL_PYR),1);
fns_neuron_PRL = cell(length(PRL_PYR),1);
plot_PRL = cell(length(PRL_PYR),2);
time_idx_PRL = [];   spk_time_PRL = cell(length(PRL_PYR),1);

icell = 94;    %%

[path, name]=fileparts(PRL_PYR{icell});
cd(path);
load('Event_AC_YJ_1220.mat');
ntrial = size(Event_Time,2);
eventtime = Event_Time;
mean_memory_length = mean(eventtime(5,:)-eventtime(4,:));

load('new_V_value_yj.mat')
QL=value_v_final(1:end-2,1);
index{1}=(QL<0.25);
index{2}=(0.25<QL).*(QL<0.5);
index{3}=(0.5<QL).*(QL<0.75);
index{4}=(0.75<QL).*(QL<1);

name=[name '.t'];
spk_data_PRL{icell} = get_Spiketime(name)*100; %s->ms
date=strsplit(path,'\');
date = date{end};
date = strsplit(date,'_');
date = date{1};
fns_neuron_PRL{icell} = [date(end-4:end) '_' name];

% PSTH
for itrial=1:ntrial
    for ievent=1:1
        [~, time_idx_PRL] = histc(spk_data_PRL{icell,1},eventtime(ievent+4,itrial)+win(ievent,:));
        if isempty(time_idx_PRL); continue; end;
        spk_time{icell,1}{itrial,ievent}=spk_data_PRL{icell,1}(logical(time_idx_PRL))-eventtime(ievent+4,itrial);
    end
end
spktime_PRL{icell,1}= {spk_time{icell,1}(index{1}(:,1)==1,:) spk_time{icell,1}(index{2}(:,1)==1,:) spk_time{icell,1}(index{3}(:,1)==1,:) spk_time{icell,1}(index{4}(:,1)==1,:)};

% Rearrange data
tmp_spk=spktime_PRL{icell};
n_trial_choice = [size(tmp_spk{1},1) size(tmp_spk{2},1) size(tmp_spk{3},1) size(tmp_spk{4},1)];
n_trial = sum(n_trial_choice);
yy = [0:n_trial-1; 1:n_trial; NaN(1,n_trial)];

for ichoice=1:4
    n_spikes_each = cellfun(@length,tmp_spk{ichoice}(:,1));
    n_spikes = sum(n_spikes_each);
    if n_spikes==0
        x_pts_PRL{icell}{ichoice} = [];
        y_pts_PRL{icell}{ichoice} = [];
        plot_PRL{icell}{ichoice} = [];
        continue;
    end
    tmp_x = [[cell2mat(tmp_spk{ichoice}(:,1))'; cell2mat(tmp_spk{ichoice}(:,1))']; NaN(1,n_spikes)];
    tmp_x = tmp_x(:); tmp_x(end) = [];
    tmp_y = [];
    tmp_x2 = zeros(n_trial_choice(ichoice),length(wintick));
    for iy=1:n_trial_choice(ichoice)
        if n_spikes_each(iy)==0; continue; end;
        tmp_x3 = tmp_spk{ichoice}{iy,1}/1000000;
        tmp_x2(itrial,:) = histc(tmp_x3(:,1),wintick);
        if ichoice==1
            tmp_y = [tmp_y repmat(yy(:,iy),1,n_spikes_each(iy))];
        else
            tmp_y = [tmp_y repmat(yy(:,iy+sum(n_trial_choice(1:ichoice-1))),1,n_spikes_each(iy))];
        end
    end
    tmp_y = tmp_y(:); tmp_y(end) = [];
    x_pts_PRL{icell}{ichoice} = tmp_x/1000000;
    y_pts_PRL{icell}{ichoice} = tmp_y;
    plot_PRL{icell}{ichoice} = {tmp_x2};
end


%% Raster PSTH그리기
lineclr={[1.0 0.8 0.8], [0.8 0.6 0.8],[0.6 0.4 0.6], [0.6 0.0 0.4]}; %choice
linewth=[0.8 0.8 0.8 0.8];

% Axes Point
npp = 1; % numbers of cells per page
axpt = AxesPoint_1raster_new(1,npp);
PSTH=cell(1,4); lim=[];
f2=figure('PaperUnits','Centimeters','PaperPosition',[0 0 7 8]);
ylims = [1 ntrial];

% spike train
axes('Position',axpt(1,mod((icell-1),npp)+1,1,:));
hold on;
set(gca, 'FontSize', 8);

line([-0.5 -0.5],[0, 900],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([0 0],[0, 900],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 900],'LineWidth',0.5,'LineStyle','-', 'Color',[120 177 255]/255);

for iichoice=1:4
    plot(x_pts_PRL{1,icell}{1,iichoice}, y_pts_PRL{1,icell}{1,iichoice},...
        'LineStyle','-','LineWidth',1,'Color',lineclr{iichoice},'Marker','o','MarkerSize',0.1);
end
set(gca,'XLim',xlims(1,:),'XTick',[],'XColor','w');
set(gca,'YLim',ylims,'YTick',ylims,'YTickLabel',{[],ylims(2)},'YColor','k');
set(gca,'box','off','TickDir','out','FontSize',8);
L1=ylabel({'Trial', 'number'}); set(L1, 'FontSize',8);
t=title('Q_{contra}');
set(t, 'FontSize', 9, 'FontAngle','italic');

% sdf
for ichoice=1:4
    tmp_x=[x_pts_PRL{icell}{ichoice};NaN];
    if isnan(tmp_x); PSTH{ichoice}=zeros(size(wintick)); continue; end;
    tmp_x2=reshape(tmp_x,3,length(tmp_x)/3);
    spike_hist = histc(tmp_x2(1,:),wintick);
    PSTH{ichoice} = conv(spike_hist,fspecial('Gaussian',[1 100],10),'same')/(bin*n_trial_choice(ichoice));
end


% PSTH with SDF
axes('Position',axpt(1,mod((icell-1),npp)+1,2,:));
hold on;

for ichoice=1:4
    plot(wintick,PSTH{ichoice},...
        'LineStyle','-','LineWidth',linewth(ichoice),'Color',lineclr{ichoice});
end
lim=max([max((PSTH{1,1}(1,:))), max((PSTH{1,2}(1,:))), max((PSTH{1,3}(1,:))) max((PSTH{1,4}(1,:)))]);
ylims2 = ceil(lim);

set(gca,'XLim',xlims(1,:), 'XTick',[xlims(1,1):1:xlims(1,2)],'XTickLabel',{});
set(gca,'YLim',[0 ylims2],'YTick',[0 ylims2],'YColor','k'); %'
set(gca,'box','off','TickDir','out','FontSize',8);

line([-0.5 -0.5],[0, 99],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([0 0],[0, 99],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 300],'LineWidth',0.5,'LineStyle','-', 'Color',[120 177 255]/255);

L1=ylabel({'Firing', 'rate (Hz)'}); set(L1, 'FontSize',8)


%% save
cd([saving_path])
%     print(gcf, '-dpdf', ['QL_Raster_' fns_neuron_PRL{icell}(1:end-2)]);
print(gcf, '-depsc','-painters', ['QL_Raster_PRL_' num2str(icell) '_' fns_neuron_PRL{icell}(1:end-2) '.ai']);

%% IL
spktime_IL = cell(length(IL_PYR),1);
spk_data_IL = cell(length(IL_PYR),1);
fns_neuron_IL = cell(length(IL_PYR),1);
plot_IL = cell(length(IL_PYR),2);
time_idx_IL = [];   spk_time_IL = cell(length(IL_PYR),1);

icell = 1;
%

[path, name]=fileparts(IL_PYR{icell});
cd(path);
load('Event_AC_YJ_1220.mat');
ntrial = size(Event_Time,2);
eventtime = Event_Time;
mean_memory_length = mean(eventtime(5,:)-eventtime(4,:));

load('new_V_value_yj.mat')
QL=value_v_final(1:end-2,1);
index{1}=(QL<0.25);
index{2}=(0.25<QL).*(QL<0.5);
index{3}=(0.5<QL).*(QL<0.75);
index{4}=(0.75<QL).*(QL<1);

name=[name '.t'];
spk_data_IL{icell} = get_Spiketime(name)*100; %s->ms
date=strsplit(path,'\');
date = date{end};
date = strsplit(date,'_');
date = date{1};
fns_neuron_IL{icell} = [date(end-4:end) '_' name];

% PSTH
for itrial=1:ntrial
    for ievent=1:1
        [~, time_idx_IL] = histc(spk_data_IL{icell,1},eventtime(ievent+4,itrial)+win(ievent,:));
        if isempty(time_idx_IL); continue; end;
        spk_time{icell,1}{itrial,ievent}=spk_data_IL{icell,1}(logical(time_idx_IL))-eventtime(ievent+4,itrial);
    end
end
spktime_IL{icell,1}= {spk_time{icell,1}(index{1}(:,1)==1,:) spk_time{icell,1}(index{2}(:,1)==1,:) spk_time{icell,1}(index{3}(:,1)==1,:) spk_time{icell,1}(index{4}(:,1)==1,:)};

% Rearrange data
tmp_spk=spktime_IL{icell};
n_trial_choice = [size(tmp_spk{1},1) size(tmp_spk{2},1) size(tmp_spk{3},1) size(tmp_spk{4},1)];
n_trial = sum(n_trial_choice);
yy = [0:n_trial-1; 1:n_trial; NaN(1,n_trial)];

for ichoice=1:4
    n_spikes_each = cellfun(@length,tmp_spk{ichoice}(:,1));
    n_spikes = sum(n_spikes_each);
    if n_spikes==0
        x_pts_IL{icell}{ichoice} = [];
        y_pts_IL{icell}{ichoice} = [];
        plot_IL{icell}{ichoice} = [];
        continue;
    end
    tmp_x = [[cell2mat(tmp_spk{ichoice}(:,1))'; cell2mat(tmp_spk{ichoice}(:,1))']; NaN(1,n_spikes)];
    tmp_x = tmp_x(:); tmp_x(end) = [];
    tmp_y = [];
    tmp_x2 = zeros(n_trial_choice(ichoice),length(wintick));
    for iy=1:n_trial_choice(ichoice)
        if n_spikes_each(iy)==0; continue; end;
        tmp_x3 = tmp_spk{ichoice}{iy,1}/1000000;
        tmp_x2(itrial,:) = histc(tmp_x3(:,1),wintick);
        if ichoice==1
            tmp_y = [tmp_y repmat(yy(:,iy),1,n_spikes_each(iy))];
        else
            tmp_y = [tmp_y repmat(yy(:,iy+sum(n_trial_choice(1:ichoice-1))),1,n_spikes_each(iy))];
        end
    end
    tmp_y = tmp_y(:); tmp_y(end) = [];
    x_pts_IL{icell}{ichoice} = tmp_x/1000000;
    y_pts_IL{icell}{ichoice} = tmp_y;
    plot_IL{icell}{ichoice} = {tmp_x2};
end

%% Raster PSTH
lineclr={[1.0 0.8 0.8], [0.8 0.6 0.8],[0.6 0.4 0.6], [0.6 0.0 0.4]}; %choice
linewth=[0.8 0.8 0.8 0.8];

% Axes Point
npp = 1; % numbers of cells per page
axpt = AxesPoint_1raster_new(1,npp);
PSTH=cell(1,4); lim=[];
f2=figure('PaperUnits','Centimeters','PaperPosition',[0 0 7 8]);
ylims = [1 ntrial];

% spike train
axes('Position',axpt(1,mod((icell-1),npp)+1,1,:));
hold on;
set(gca, 'FontSize', 8);

line([-0.5 -0.5],[0, 900],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([0 0],[0, 900],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 900],'LineWidth',0.5,'LineStyle','-', 'Color',[120 177 255]/255);

for iichoice=1:4
    plot(x_pts_IL{1,icell}{1,iichoice}, y_pts_IL{1,icell}{1,iichoice},...
        'LineStyle','-','LineWidth',1,'Color',lineclr{iichoice},'Marker','o','MarkerSize',0.1);
end
set(gca,'XLim',xlims(1,:),'XTick',[],'XColor','w');
set(gca,'YLim',ylims,'YTick',ylims,'YTickLabel',{[],ylims(2)},'YColor','k');
set(gca,'box','off','TickDir','out','FontSize',8);

L1=ylabel({'Trial', 'number'}); set(L1, 'FontSize',8);
t=title('Q_{contra}');
set(t, 'FontSize', 9,'FontAngle','italic');

% sdf
for ichoice=1:4
    tmp_x=[x_pts_IL{icell}{ichoice};NaN];
    if isnan(tmp_x); PSTH{ichoice}=zeros(size(wintick)); continue; end;
    tmp_x2=reshape(tmp_x,3,length(tmp_x)/3);
    spike_hist = histc(tmp_x2(1,:),wintick);
    PSTH{ichoice} = conv(spike_hist,fspecial('Gaussian',[1 100],10),'same')/(bin*n_trial_choice(ichoice));
end

% PSTH with SDF
axes('Position',axpt(1,mod((icell-1),npp)+1,2,:));
hold on;

for ichoice=1:4
    plot(wintick,PSTH{ichoice},...
        'LineStyle','-','LineWidth',linewth(ichoice),'Color',lineclr{ichoice});
end
lim=max([max((PSTH{1,1}(1,:))), max((PSTH{1,2}(1,:))), max((PSTH{1,3}(1,:))) max((PSTH{1,4}(1,:)))]);
ylims2 = ceil(lim);

set(gca,'XLim',xlims(1,:), 'XTick',[xlims(1,1):1:xlims(1,2)],'XTickLabel',{'' -1 0 1 ''});
set(gca,'YLim',[0 ylims2],'YTick',[0 ylims2],'YColor','k'); %'
set(gca,'box','off','TickDir','out','FontSize',8);

line([-0.5 -0.5],[0, 99],'LineWidth',0.5,'LineStyle','--', 'Color',[.6 .6 .6]);
line([0 0],[0, 99],'LineWidth',0.5,'LineStyle','-', 'Color',[.6 .6 .6]);
line([-mean_memory_length/1000000 -mean_memory_length/1000000], [0 300],'LineWidth',0.5,'LineStyle','-', 'Color',[120 177 255]/255);

L1=ylabel({'Firing', 'rate (Hz)'}); set(L1, 'FontSize',8)
L2=xlabel('Time from reward onset (s)'); set(L2, 'FontSize',9)

%% save
cd([saving_path])
%     print(gcf, '-dpdf', ['QL_Raster_' fns_neuron_IL{icell}(1:end-2)]);
print(gcf, '-depsc','-painters', ['QL_Raster_IL_' num2str(icell) '_' fns_neuron_IL{icell}(1:end-2) '.ai']);
close all;
