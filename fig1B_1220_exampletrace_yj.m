close all; clear all;

%% Loading video files
cd('E:\Data\Neural analysis\ACC\2016-03-29_12-31-51'); %LE9_5
[T,P,angle] = nvt2mat('VT1.nvt');

% beh data
load('beh_total.mat');
beh=beh_{1}(:,1:3); clear beh_,

% Event data
load('Event.mat');     clear basebin, clear index,

if Event_f(end,2) ~= 5
    Event_f(end,:) = [];
end
if Event_f(end,2) ~= 5
    Event_f(end,:) = [];
end
if Event_f(end,2) ~= 5
    Event_f(end,:) = [];
end

ntrial=length(Event_f)/4;
eventtime = [reshape(Event_f(1:end,1), 4, ntrial) [Event_f(end,1);0;0;0]];
eventflag = [reshape(Event_f(1:end,2), 4, ntrial) [Event_f(end,2);0;0;0]];
RT_win=[];
RT_win(2,:)=eventtime(1,2:end); 
RT_win(1,:)=eventtime(1,1:end-1); 
p=200;

% VT data
x_posi = P(:,1);
y_posi = P(:,2);
clear Target, clear dwPoints,

%% Converging time
position={}; con_time={}; total_xL=[]; total_xR=[]; total_yL=[]; total_yR=[];
Rt=0; Lt=0;
for itrial=1:size(RT_win,2)
    [~, VT_hist]=histc(T,RT_win(:,itrial));
    if isempty(find(VT_hist==1,1)); continue; end;
    position{itrial,1}= [find(VT_hist==1), x_posi(VT_hist==1),y_posi(VT_hist==1)];
    con_time{itrial,1}=T(position{itrial,1}(:,1))/1000 - T(min(position{itrial,1}(:,1)))/1000;
    if beh(itrial,1)<0
        Lt=Lt+1;
        total_xL{Lt,1}= position{itrial,1}(:,2);
        total_yL{Lt,1}= position{itrial,1}(:,3);
    else
        Rt=Rt+1;
        total_xR{Rt,1}= position{itrial,1}(:,2);
        total_yR{Rt,1}= position{itrial,1}(:,3);
    end
end
clear Lt, clear Rt,

%% plot
figure()
axes('Position',[0.2 0.1 0.2 0.8]);
hold on;
for ii=1:size(total_xR,1)
    h1=plot(total_xR{ii,1},total_yR{ii,1}, 'LineStyle','none');
    set(h1, 'Marker','.','MarkerSize', 3,'MarkerEdgeColor',[0.7 0.0 0.0]);
    set(gca, 'YDir', 'reverse');
    hold on
end
for ii=1:size(total_xL,1)
    h1=plot(total_xL{ii,1},total_yL{ii,1}, 'LineStyle','none');
    set(h1,'Marker','.', 'MarkerSize', 3, 'MarkerEdgeColor',[0.0 0.0 0.4]);
    set(gca, 'YDir', 'reverse');
    hold on
end

set(gca,'XLim',[150 600], 'XTick',[150:50:600]);
set(gca, 'ylim',[0 480],'YTick',[0:10:480]);
axis off;

%% save
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters','Fig1_B_tracking_LE9_5.ai')
% cd('D:\KAIST\Grad\SNL\LDH_yj\Figures')
% print(gcf,'-dtiff','-r300', ['Fig1_B_tracking_LE9_23']);
