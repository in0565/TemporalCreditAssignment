close all; clear all;

%% Loading video files
cd_path = 'E:\Data\Neural analysis';
cd(cd_path);
videofiles = FindFiles('*.nvt');

%%
for ifile=1:length(videofiles)
    cd(fileparts(videofiles{ifile}));
    [T,P,angle] = nvt2mat('VT1.nvt'); %add by JWlee
    % beh data
    load('beh_.mat');
    beh=beh_{1,2}(:,1:3); clear beh_,
    %% Event data
    load('Event_New_YJ_1220.mat');
    RT_win=New_Event_Time(2,1:end-1); % definition: delay offset
    RT_win(1,:)=RT_win(1,:)-1000000;
    RT_win(2,:)=RT_win(1,:)+4000000; % 4s window: -1 ~ 3 of delay offset
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
        con_time{itrial,1}=T(position{itrial,1}(:,1))/1000 - T(min(position{itrial,1}(:,1)))/1000;
        if beh(itrial,1)<0 % position{itrial,1}(p,2)>330
            total_L= [total_L; position{itrial,1}(1:p,2)'];
        else
            total_R= [total_R; position{itrial,1}(1:p,2)'];
        end
    end
    %% mean ttest
    total_mean_L=mean(total_L);
    total_mean_R=mean(total_R);
    total_std_L=nanstd(total_L);
    total_std_R=nanstd(total_R);
    for ii=1:p
        ttest_result(1,ii)=ttest2(total_L(:,ii), total_R(:,ii));
    end
    ttest_result_fig = ttest_result;
    ttest_result_fig(ttest_result_fig==0)=nan;
    
    for i2 = 31:p-2
        if isequal(ttest_result(1,i2:i2+2),[1 1 1]) && mean(total_L(:,i2))<mean(total_R(:,i2))
            point_t = i2;
            no_conv = false;
            break
        elseif i2==p-2
            no_conv = true;
            disp('no divergence')
            disp(ifile)
        end
    end
    if no_conv
        Event_Time=New_Event_Time(:,1:end-1);
        save(['Event_AC_YJ_1220'], 'Event_Time');
        clear Event_Time
        continue
    end
    %% plot
    figure();
    axes('Position',[0.2 0.2 0.6 0.6]);    hold on;
    line([1000 1000],[0, 550],'LineStyle','-','LineWidth',1.3, 'Color',[.8 .8 .8]);     hold on;
    line([2000 2000],[0, 550],'LineStyle','-','LineWidth',1.3, 'Color',[.8 .8 .8]);
    errorbar(con_time{itrial,1}(1:p,1),total_mean_L,total_std_L,'Marker','o', 'MarkerFaceColor',[0.0 0.0 0.4],'MarkerEdgeColor',[0.0 0.0 0.4],'LineWidth',1.3,'Color',[0.0 0.0 0.4]); hold on;
    errorbar(con_time{itrial,1}(1:p,1),total_mean_R,total_std_R,'Marker','o', 'MarkerFaceColor',[0.7 0.0 0.0],'MarkerEdgeColor',[0.7 0.0 0.0],'LineWidth',1.3,'Color',[0.7 0.0 0.0]); hold on;
    set(gca,'XLim',[500 3000],'XTick',[500:500:3000], 'XTickLabel',{''; 0; ''; 1; ''; 2});     ylim([150 500]);
    t=xlabel('Time from delay offset (s)'); t1=ylabel('X-position');
    set(t, 'FontSize', 15); set(t1, 'FontSize', 17);
    line([con_time{itrial,1}(point_t,1) con_time{itrial,1}(point_t,1)],[0, 550],'LineStyle','-','LineWidth',1.3, 'Color',[1.0 0.0 0.2]);
    plot(con_time{itrial,1}(1:p,1),500*ttest_result_fig,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')

    print(gcf,'-dtiff', ['diversing_ttest_yj_1220']);
    %% x y T point
    cross1=[];
    for itrial=1:size(RT_win,2)
        cross1(itrial,1:3)=position{itrial,1}(point_t,:);
        cross1(itrial,4)=T(cross1(itrial,1));
    end
    %% Saving new Event time
    Event_Time(1,:)=New_Event_Time(1,1:end-1); % delay onset
    Event_Time(2,:)=cross1(:,4)'; % approach onset
    Event_Time(3:5,:)=New_Event_Time(3:5,1:end-1); % selection, memory, reward onset
    save(['Event_AC_YJ_1220'], 'Event_Time');
    clear Event_Time
    close all
end
