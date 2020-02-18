close all; clear all;

%% Loading video files

cd('E:\Data\Neural analysis\ACC\2016-03-29_12-31-51'); %LE9_5

filename = 'fig1C_conversing_ttest_all_1220_YJ_2016-0329';

[T,P,angle] = nvt2mat('VT1.nvt'); %add by JWlee


%% beh data
load('beh_total.mat');
beh=beh_{1}(:,1:3); clear beh_,
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
ttest_result = zeros(1,p);
for ii=1:p
    ttest_result(1,ii)=ttest2(total_L(:,ii),total_R(:,ii));
end
ttest_result_fig = ttest_result;
ttest_result_fig(ttest_result_fig==0)=nan;

for i2 = 1:p-2
    if isequal(ttest_result(1,i2:i2+2),[0 0 0])
        point_t = i2;
        break
    elseif i2==p-2
        disp('no convergence')
    end
end

%% plot
figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 5.25]);
axes('Position',[0.2 0.2 0.6 0.6]);    hold on;
for itrial=1:size(RT_win,2)
    if beh(itrial,1)<0
        plot(con_time{itrial,1}(:,1)+20,position{itrial,1}(:,2),'Marker','.', 'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',0.35,'Color',[0.0 0.0 0.4]); hold on;
    else
        plot(con_time{itrial,1}(:,1)+20,position{itrial,1}(:,2),'Marker','.', 'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',0.35,'Color',[0.7 0.0 0.0]); hold on;

    end
end
set(gca,'XLim',[0 3000],'XTick',[0:1000:3000], 'XTickLabel',[-3:1:0], 'FontSize', 11);
set(gca,'YLim', [200 600], 'YTick',[200 600], 'YTickLabel',[200 600]);
line([con_time{itrial,1}(point_t,1) con_time{itrial,1}(point_t,1)],[0, 1000],'LineStyle','-','LineWidth',0.5, 'Color',[1.0 0.0 0.2]);

%% save
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters',[filename '.ai']);
% cd('D:\KAIST\Grad\SNL\LDH_yj\Figures')
% print('-dpdf', filename);
% print(gcf,'-dtiff', filename);
