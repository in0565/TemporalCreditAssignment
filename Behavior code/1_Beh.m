clc; clear all; close all;

%% Make Beh-file each day
cd('E:\Data\Neural analysis\M2\2016-07-29_14-52-23')
behfiles = FindFiles('*_LE16.txt');

beh_shuff=[]; 
for ifile=1:length(behfiles)
 beh_={};   
    [path{ifile}, name{ifile}, ext]=fileparts(behfiles{ifile});
    cd(path{ifile});
    beh_{1,1}=name{ifile};  %%index1: beh_name
%% Load txt. files wo error    
    beh_ori=load(behfiles{ifile}); 
    beh=[];
    beh=beh_ori(1:end,3:4);
    Err=find(beh_ori(:,2)==1);
    beh(Err,:)=[]; 
%% Lt: 1 -> -1, Rt: 2 -> 1, non-Rw -> -1, Rw ->1 
   beh(:,1)=beh(:,1)*2-3;
   beh(:,2)=beh(:,2)*2-1;
   beh(:,3)=beh(:,1).*beh(:,2);
   
      % load('Event_New_5.mat');  RT=[];
      load('Event_New_YJ_1220.mat');  RT=[];
      
      RT=[RT;(New_Event_Time(5,:)-New_Event_Time(3,:))'];
   %RT¿Í behÀÇ Â÷ÀÌ º¸Á¤   
    if length(RT) > length(beh)  
        RT=RT(1:end-(length(RT)-length(beh)),:);
    elseif length(RT) < length(beh)  
        beh=beh(1:end-(length(beh)-length(RT)),:);
    end    
   beh(:,4)=RT./1000000;   
   beh(:,5)=log(RT./1000000);    
      
%% index2: t=0(C,Rw,X), index3: t=-1(C,Rw,X), index4:t=-2
    beh_{1,2}=beh; % C,R,X
    beh_{1,3}=[[0,0,0,0,0]; beh(1:end-1,:)];
    beh_{1,4}=[[0,0,0,0,0]; [0,0,0,0,0]; beh(1:end-2,:)];    
%% index5: shuffle
    trials=length(beh);
    shuffle=randperm(trials)';
    cho=beh(:,1); rew=beh(:,2); int=beh(:,3);
    beh_shuff(:,1)=cho(shuffle); beh_shuff(:,2)=rew(shuffle); beh_shuff(:,3)=int(shuffle); 
    
    beh_{1,5}=beh_shuff;
    beh_{1,6}=[[0,0,0]; beh_shuff(1:end-1,:)];
    beh_{1,7}=[[0,0,0]; [0,0,0]; beh_shuff(1:end-2,:)];  
    
%cd('F:\ES_data_forH\beh\beh_');
 cd(path{ifile})
save(['beh_.mat'], 'beh_')
clear beh_shuff, 
end
