close all; clear all;
%%
load('E:\Data\Neural analysis\mean_memory_lengths_1220.mat');
mean_memory_length = nanmean(memory_length_1220);

%% Loading video files
sub_path = {'ACC','PRL','IL'};
session_mean_R = {};
session_mean_L = {};
isession = 0;
for iregion = 1:3
    cd_path = ['E:\Data\Neural analysis\' sub_path{iregion}];
    cd(cd_path);
    videofiles = FindFiles('*.nvt'); % mcluster가 directory에 있어야함
    
    %% Converging point for each session
    for ifile=1:length(videofiles)
        cd(fileparts(videofiles{ifile}));
        [T,P,angle] = nvt2mat('VT1.nvt'); %add by JWlee
        isession = isession+1;
        % beh data
        load('beh_.mat');
        beh=beh_{1,2}(:,1:3); clear beh_,
        %% Event data
        load('Event.mat');     clear basebin, clear index,
        
        if Event_f(end,2) ~= 5   ;
            Event_f(end,:) = [];
        end
        
        if Event_f(end,2) ~= 5;
            Event_f(end,:) = [];
        end
        
        if Event_f(end,2) ~= 5;
            Event_f(end,:) = [];
        end
        
        ntrial=length(Event_f)/4;
        eventtime = [reshape(Event_f(1:end,1), 4, ntrial)  [Event_f(end,1);0;0;0]];
        eventflag = [reshape(Event_f(1:end,2), 4, ntrial) [Event_f(end,2);0;0;0]];
        RT_win=[];
        RT_win(2,:)=eventtime(4,1:end-1)+1500000;
        RT_win(1,:)=eventtime(4,1:end-1)-2500000;
        p=99;
        %% VT data
        x_posi=P(:,1); y_posi=P(:,2);
        clear Target, clear dwPoints,
        %% Converging time
        position={}; con_time={}; total_L=[]; total_R=[];
        for itrial=1:size(RT_win,2)
            [~, VT_hist]=histc(T,RT_win(:,itrial));
            if isempty(find(VT_hist==1,1)); continue; end;
            position{itrial,1}= [find(VT_hist==1), x_posi(VT_hist==1),y_posi(VT_hist==1)];
            position{itrial,1}(position{itrial,1}==0)=NaN;
            con_time{itrial,1}=T(position{itrial,1}(:,1))/1000 - T(min(position{itrial,1}(:,1)))/1000;
            if beh(itrial,1)<0
                total_L= [total_L; position{itrial,1}(1:p,2)'];
            else
                total_R= [total_R; position{itrial,1}(1:p,2)'];
            end
        end
        
        %% mean ttest
        total_mean_L=nanmean(total_L);
        total_mean_R=nanmean(total_R);
        session_mean_L{isession,1}=total_mean_L(1:p);
        session_mean_R{isession,1}=total_mean_R(1:p);
    end
end

%% rearrange
session_mean_L_mat = cell2mat(session_mean_L);
session_mean_R_mat = cell2mat(session_mean_R);

%% plot
figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 5.25]);
axes('Position',[0.2 0.3 0.6 0.6]);    hold on;
fill_x = [con_time{itrial,1}(1:p,1)' fliplr(con_time{itrial,1}(1:p,1)')];
fill_Ly = [mean(session_mean_L_mat,1)+std(session_mean_L_mat,0,1),fliplr(mean(session_mean_L_mat,1)-std(session_mean_L_mat,0,1))];
fill_Ry = [mean(session_mean_R_mat,1)+std(session_mean_R_mat,0,1),fliplr(mean(session_mean_R_mat,1)-std(session_mean_R_mat,0,1))];
fill_L = fill(fill_x,fill_Ly,[153 153 255]/256); alpha(fill_L,0.5); set(fill_L,'EdgeColor','none');
fill_R = fill(fill_x,fill_Ry,[255 102 102]/256); alpha(fill_R,0.5); set(fill_R,'EdgeColor','none');
plot(con_time{itrial,1}(1:p,1),mean(session_mean_L_mat,1),'Marker','none', 'MarkerFaceColor',[0.0 0.0 0.4],'MarkerEdgeColor',[0.0 0.0 0.4],'LineWidth',0.8,'Color',[0.0 0.0 0.4]); hold on;
plot(con_time{itrial,1}(1:p,1),mean(session_mean_R_mat,1),'Marker','none', 'MarkerFaceColor',[0.7 0.0 0.0],'MarkerEdgeColor',[0.7 0.0 0.0],'LineWidth',0.8,'Color',[0.7 0.0 0.0]); hold on;
set(gca,'XLim',[0 3000],'XTick',[0:1000:3000], 'XTickLabel',[-3:1:0],'FontSize',11);
set(gca,'YLim', [200 600], 'YTick',[200 600], 'YTickLabel',[200 600]);
line([3000-mean_memory_length/1000 3000-mean_memory_length/1000],[0, 1000],'LineStyle','-','LineWidth',0.5, 'Color',[1.0 0.0 0.2]);
t=xlabel('Time from reward onset (s)');
set(t, 'FontSize', 12);

%% save
% cd('D:\KAIST\Grad\SNL\LDH_yj\Figures')
% print(gcf,'-dtiff', ['fig1C_conversing_sessions_YJ_1215']);
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters',['fig1C_conversing_sessions_YJ_1220' '.ai']);
