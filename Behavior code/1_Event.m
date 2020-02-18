clc; close all; clear all;

%% Loading event files
cd('E:\Data\Neural analysis\LE17\IL');
eventfiles = FindFiles('*.nev');

event_time = []; event_flag = []; basebin=[];
for ifile=1:length(eventfiles)
    cd(fileparts(eventfiles{ifile}));
    [event_time_full, event_flag_full] = EventRead;
     
    basebin=[event_time_full(1) event_time_full(2); event_time_full(end-1) event_time_full(end)];
    
    event_time = event_time_full(4:end-3); event_time = event_time(1:2:end);
    event_flag = event_flag_full(4:end-3); event_flag = event_flag(1:2:end);
    
    HH=length(event_flag); event_flag1 = zeros(size(event_flag,1),1);

  for num_str = 1:HH     
    if strcmp(event_flag{num_str,1}(end-7:end-2),num2str('0x0008')) == 1  
            event_flag1(num_str,1)  = 1;
      elseif strcmp(event_flag{num_str,1}(end-7:end-2),num2str('0x0002')) == 1 
             event_flag1(num_str,1)  = 2;
      elseif strcmp(event_flag{num_str,1}(end-7:end-2),num2str('0x0040')) == 1  %(2)
             event_flag1(num_str,1)  = 3;
      elseif strcmp(event_flag{num_str,1}(end-7:end-2),num2str('0x0020')) == 1 %(2)
             event_flag1(num_str,1)  = 4;
      elseif strcmp(event_flag{num_str,1}(end-7:end-2),num2str('0x0004')) == 1 %(1)
             event_flag1(num_str,1)  = 5; 
    end
  end
  
   if event_flag1(1,1) == 2    
             event_flag1=[1; event_flag1]; 
             event_time=[event_time_full(3,1); event_time];
             HH=length(event_flag1);
   end
    all_flag1 = [event_time, event_flag1];

    E=find(event_flag1==1);     Err=find(diff(E)==2);   E(Err)=[];
    EE=sort([E;E+1;E+2;E+3]);  EE=EE(1:end-4);
    Event_f=all_flag1(EE,:);
    
   
    A=Event_f(3:4:end,:);

  index{1}=find(A(:,2)==3);   index{2}=find(A(:,2)==4);
  %% Saving event files
     cd(fileparts(eventfiles{ifile}));
     save(['Event.mat'], 'Event_f','basebin','index')
end
