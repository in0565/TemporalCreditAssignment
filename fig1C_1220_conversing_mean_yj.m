close all; clear all;

%% Loading video files
cd_path = 'E:\Data\Neural analysis';
cd(cd_path);
videofiles = FindFiles('*.nvt'); % mcluster가 directory에 있어야함

%% Converging point for each session
memory_length_1220 = zeros(1,length(videofiles));
invalid_sessions_1220 = {};
for ifile=1:length(videofiles)
    cd(fileparts(videofiles{ifile}));
    [T,P,angle] = nvt2mat('VT1.nvt'); %add by JWlee
    
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
    RT_win(2,:)=eventtime(4,1:end-1)+1500000; % rw 도달 후 1.5초
    RT_win(1,:)=eventtime(4,1:end-1)-2500000; % rw 도달 전 2.5초 = rw onset전 3초
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
    total_std_L=nanstd(total_L);
    total_std_R=nanstd(total_R);
    
    ttest_result = zeros(1,p);
    for ii=1:p
        ttest_result(1,ii)=ttest2(total_L(:,ii),total_R(:,ii));
    end
    ttest_result_fig = ttest_result;
    ttest_result_fig(ttest_result_fig==0)=nan;
    
    for i2 = 1:p-2
        if isequal(ttest_result(1,i2:i2+2),[0 0 0])
            point_t = i2;
            no_conv = false;
            break
        elseif i2==p-2
            no_conv = true;
            disp('no convergence')
            disp(ifile)
            disp(fileparts(videofiles{ifile}))
            memory_length_1220(ifile) = nan;
            invalid_sessions_1220{end+1,1} = fileparts(videofiles{ifile});
        end
    end
    if no_conv
        continue
    end
    %% plot
    figure();
    axes('Position',[0.2 0.2 0.6 0.6]);    hold on;
    line([1000 1000],[0, 550],'LineStyle','-','LineWidth',1.3, 'Color',[.8 .8 .8]);
    line([2000 2000],[0, 550],'LineStyle','-','LineWidth',1.3, 'Color',[.8 .8 .8]);
    errorbar(con_time{itrial,1}(1:p,1),total_mean_L(1:p),total_std_L(1:p),'Marker','o', 'MarkerFaceColor',[0.0 0.0 0.4],'MarkerEdgeColor',[0.0 0.0 0.4],'LineWidth',1.3,'Color',[0.0 0.0 0.4]); hold on;
    errorbar(con_time{itrial,1}(1:p,1),total_mean_R(1:p),total_std_R(1:p),'Marker','o', 'MarkerFaceColor',[0.7 0.0 0.0],'MarkerEdgeColor',[0.7 0.0 0.0],'LineWidth',1.3,'Color',[0.7 0.0 0.0]); hold on;

    set(gca,'XLim',[0 3000],'XTick',[0:1000:3000], 'XTickLabel',[-3:1:0]);     ylim([200 550]);
    t=xlabel('Time from reward onset (s)'); t1=ylabel('X-position');
    set(t, 'FontSize', 15); set(t1, 'FontSize', 17);
    line([con_time{itrial,1}(point_t,1) con_time{itrial,1}(point_t,1)],[0, 550],'LineStyle','-','LineWidth',1.3, 'Color',[1.0 0.0 0.2]);
    plot(con_time{itrial,1}(1:p,1),550*ttest_result_fig,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
    
    
    print(gcf,'-dtiff', ['conversing_ttest_YJ_1215']);
    
    %% x y T point
    cross2=[];
    for itrial=1:size(RT_win,2)
        cross2(itrial,1:3)=position{itrial,1}(point_t,:);
        cross2(itrial,4)=T(cross2(itrial,1));
    end
    %% Saving new Event time
    New_Event_Time(1:3,:)=eventtime(1:3,1:end-1); % delay onset, delay offset, selection onset
    New_Event_Time(4,:)=cross2(:,4)'; % memory onset
    New_Event_Time(5,:)=eventtime(4,1:end-1)+500000; % reward onset
    if (mean(New_Event_Time(5,:)-New_Event_Time(4,:)))/1000000<0.5
        invalid_sessions_1220{end+1,1} = fileparts(videofiles{ifile});
        memory_length_1220(ifile) = nan;
    else
        memory_length_1220(ifile) = mean(New_Event_Time(5,:)-New_Event_Time(4,:));
    end
    save(['Event_New_YJ_1220'], 'New_Event_Time');
    clear New_Event_Time point_t
    close all
end

%% save memory stage length
cd(cd_path);
save('mean_memory_lengths_1220','memory_length_1220');
save('invalid_sessions_1220','invalid_sessions_1220');
