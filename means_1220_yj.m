clear; close all;
%%
main_path = 'E:\Data\Neural analysis';
sub_path = {'\ACC','\IL','\PRL'};
cd(main_path)

%% trial number per session: event file
%% reward state time: event file
trial_number = []; rewarded_reward_duration = []; unrewarded_reward_duration = [];
for iregion = 1:3
    cd([main_path sub_path{iregion}])
    videofiles = FindFiles('*.nvt');
    for ifile=1:length(videofiles)
        % trial number
        cd(fileparts(videofiles{ifile}));
        load('Event_AC_YJ_1220.mat')
        trial_number = [trial_number; size(Event_Time,2)];
        
        % reward state
        load('Event_New_YJ_1220.mat');
        reward_duration = New_Event_Time(1,2:end)-New_Event_Time(5,1:end-1);
        % rewarded/unrewarded
        load('beh_.mat');
        rewarded = beh_{1,2}(1:end-1,2);
        rewarded_reward_duration = [rewarded_reward_duration reward_duration(rewarded==1)];
        unrewarded_reward_duration = [unrewarded_reward_duration reward_duration(rewarded==-1)];
    end
end
rewarded_reward_duration=rewarded_reward_duration/1000000;
unrewarded_reward_duration=unrewarded_reward_duration/1000000;
fprintf('Trial number per session: %.3f +- %.3f (SD)\n',mean(trial_number),std(trial_number))
fprintf('Rewarded trial reward stage duration: %.3f +- %.3f (SD) sec\n',mean(rewarded_reward_duration),std(rewarded_reward_duration))
fprintf('Unrewarded trial reward stage duration: %.3f +- %.3f (SD) sec\n',mean(unrewarded_reward_duration),std(unrewarded_reward_duration))
fprintf('SEM (s): %.3f rewarded / %.3f unrewarded\n',std(rewarded_reward_duration)/sqrt(length(rewarded_reward_duration)),std(unrewarded_reward_duration)/sqrt(length(unrewarded_reward_duration)))

%% firing rate during task per region and neuron
ACC_type={}; PRL_type={}; IL_type={};
ACC_PYR=[]; PRL_PYR=[]; IL_PYR=[];
ACC_INT=[]; PRL_INT=[]; IL_INT=[];

cd('E:\Data\CRX\MeanFR\reward01');
ACC_type = load('celltype_ACC_YJ_1220.mat');       ACC_PYR = ACC_type.type{1}(:,1);
PRL_type = load('celltype_PRL_YJ_1220.mat');       PRL_PYR = PRL_type.type{1}(:,1);
IL_type = load('celltype_IL_YJ_1220.mat');         IL_PYR = IL_type.type{1}(:,1);

clear ACC_type, clear PRL_type, clear IL_type,
%% ACC
spk_data_ACC=cell(length(ACC_PYR),1); mean_FR_ACC = zeros(length(ACC_PYR),1);
for icell=1:length(ACC_PYR)
    [path, name]=fileparts(ACC_PYR{icell});
    cd(path);
    load('Event_New_YJ_1220.mat');
    % between New_Event_Time(1,1) New_Event_Time(1,end)
    session_onset = New_Event_Time(1,1);
    session_offset = New_Event_Time(1,end);
    %% Loading Spike datad
    name=[name '.t'];
    spk_data_ACC{icell} = get_Spiketime(name)*100;
    temp_spk = spk_data_ACC{icell};
    %% firing rate
    session_spk = temp_spk(temp_spk>=session_onset & temp_spk<=session_offset);
    mean_FR_ACC(icell,1) = length(session_spk)/(session_offset-session_onset)*1000000;
end
fprintf('Mean firing rate in ACC (%d neurons): %.3f +- %.3f Hz\n',length(mean_FR_ACC),mean(mean_FR_ACC),std(mean_FR_ACC));

%% PRL
spk_data_PRL=cell(length(PRL_PYR),1); mean_FR_PRL = zeros(length(PRL_PYR),1);
for icell=1:length(PRL_PYR)
    [path, name]=fileparts(PRL_PYR{icell});
    cd(path);
    load('Event_New_YJ_1220.mat');
    % between New_Event_Time(1,1) New_Event_Time(1,end)
    session_onset = New_Event_Time(1,1);
    session_offset = New_Event_Time(1,end);
    %% Loading Spike datad
    name=[name '.t'];
    spk_data_PRL{icell} = get_Spiketime(name)*100;
    temp_spk = spk_data_PRL{icell};
    %% firing rate
    session_spk = temp_spk(temp_spk>=session_onset & temp_spk<=session_offset);
    mean_FR_PRL(icell,1) = length(session_spk)/(session_offset-session_onset)*1000000;
end
fprintf('Mean firing rate in PRL (%d neurons): %.3f +- %.3f Hz\n',length(mean_FR_PRL),mean(mean_FR_PRL),std(mean_FR_PRL));

%% IL
spk_data_IL=cell(length(IL_PYR),1); mean_FR_IL = zeros(length(IL_PYR),1);
for icell=1:length(IL_PYR)
    [path, name]=fileparts(IL_PYR{icell});
    cd(path);
    load('Event_New_YJ_1220.mat');
    % between New_Event_Time(1,1) New_Event_Time(1,end)
    session_onset = New_Event_Time(1,1);
    session_offset = New_Event_Time(1,end);
    %% Loading Spike datad
    name=[name '.t'];
    spk_data_IL{icell} = get_Spiketime(name)*100;
    temp_spk = spk_data_IL{icell};
    %% firing rate
    session_spk = temp_spk(temp_spk>=session_onset & temp_spk<=session_offset);
    mean_FR_IL(icell,1) = length(session_spk)/(session_offset-session_onset)*1000000;
end
fprintf('Mean firing rate in IL (%d neurons): %.3f +- %.3f Hz\n',length(mean_FR_IL),mean(mean_FR_IL),std(mean_FR_IL));

%% all neurons
mean_FR_all = [mean_FR_ACC; mean_FR_PRL; mean_FR_IL];
fprintf('Mean firing rate all PYR (%d neurons): %.3f +- %.3f Hz\n',length(mean_FR_all),mean(mean_FR_all),std(mean_FR_all));
