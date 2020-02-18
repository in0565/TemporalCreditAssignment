clear all; close all;
tic

cd('E:\Data\CRX\MeanFR\reward01');
ACC_type = load('celltype_ACC_YJ_1220.mat');       ACC_PYR = ACC_type.type{1}(:,1);
PRL_type = load('celltype_PRL_YJ_1220.mat');       PRL_PYR = PRL_type.type{1}(:,1);
IL_type = load('celltype_IL_YJ_1220.mat');         IL_PYR = IL_type.type{1}(:,1);

clear ACC_type, clear PRL_type, clear IL_type,
saving_path = 'D:\KAIST\Grad\SNL\LDH_yj\mat\1220_reward01';

%% window& bin for Fig.
win_title = 2.0;
win = [0 win_title]*1000000;
iter = 100;
%% ACC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Loading Beh(for reward) & Event
spk_data_ACC=cell(length(ACC_PYR),1);  variable_ACC=cell(length(ACC_PYR),1);  spike_hist_ACC_reward={};
variable_sh_ACC=cell(length(ACC_PYR),iter); ran_y_ACC = cell(length(ACC_PYR),iter);
for icell=1:length(ACC_PYR)
    [path, name]=fileparts(ACC_PYR{icell});
    cd(path);
    load('Event_AC_YJ.mat');
    eventtime =[];
    eventtime(1:5,:) = Event_Time(1:5,:);
    ntrial=size(eventtime,2);

    valid_trial = (1:ntrial);
    %% Variable 1:L, 2:R, 3:Qc (delta_Q)  4,5,6: CRX
    load('beh_.mat');
    ran_trial = size(beh_{1,2},1);
    variable_ACC{icell,1}(:,1:3)=beh_{1,2}([valid_trial ran_trial],1:3);
    load('new_V_value_yj.mat');
    variable_ACC{icell,1}(:,4:6)=value_v_final([valid_trial ran_trial],1:3);
    
    %% shuffle Y
    for iran = 1:iter
        ran_y = randperm(ntrial);
        ran_y_ACC{icell,iran}=ran_y;
        variable_sh_ACC{icell,iran}(:,4:6)=value_v_final([valid_trial(ran_y) ran_trial],1:3);
        variable_sh_ACC{icell,iran}(:,1:3)=beh_{1,2}([valid_trial(ran_y) ran_trial],1:3);
    end
   %% Loading Spike datad
    name=[name '.t'];
    spk_data_ACC{icell} = get_Spiketime(name)*100;
    %% Rearrange data
    for itrial=1:ntrial
        ievent = 5;
        [spk_number, time_idx] = histc(spk_data_ACC{icell,1},eventtime(ievent,itrial)+win(1,:));
        spike_hist_ACC_reward{icell,1}(itrial,1)=spk_number(1);
    end
end
%% checking error
disp('ACC_PYR_checking error....');
for icell=1:length(ACC_PYR)
    itype=5;
        if length(spike_hist_ACC_reward{icell}(:,1))+1 == length(variable_ACC{icell}(:,1));
        else
            delta= length(spike_hist_ACC_reward{icell}(:,1))+1 - length(variable_ACC{icell}(:,1));
            spike_hist_ACC_reward{icell}= spike_hist_ACC_reward{icell}(:,1:end-delta);
            fprintf('Delta of %d in cell %d',delta,icell);
        end
end

%% GLM fitting
disp('ACC_PYR_GLM fitting....');
P_value_ACC=cell(length(ACC_PYR),1); P_value_sh_ACC=cell(length(ACC_PYR),iter);
SRC_ACC=cell(length(ACC_PYR),1); SRC_sh_ACC=cell(length(ACC_PYR),iter);
for icell=1:length(ACC_PYR)
    %% for each neuron
    fprintf('Fitting cell no %d\n',icell)
    std_x = std(variable_ACC{icell}(1:end-1,:));
    
    variable_ACC{icell}(:,7)=[[0];spike_hist_ACC_reward{icell}(1:end-1,1);[0]];
    variable_ACC{icell}(:,8)=[[0;0];spike_hist_ACC_reward{icell}(1:end-2,1);[0]];
    variable_ACC{icell}(:,9)=[[0;0;0];spike_hist_ACC_reward{icell}(1:end-3,1);[0]];
    [coeff, dev, stats] = glmfit(variable_ACC{icell}(1:end-1,:), spike_hist_ACC_reward{icell}, 'normal');
    P_value_ACC{icell,1}(1:6,1) = stats.p(2:7,:);
    std_y = std(spike_hist_ACC_reward{icell});
    SRC_ACC{icell,1}(1:6,1) = coeff(2:7,1).*std_x(1:6)'./std_y;
    
    for iran=1:iter
        variable_sh_ACC{icell,iran}(:,7)=[[0];spike_hist_ACC_reward{icell}(ran_y_ACC{icell,iran}(1,1:end-1),1);[0]];
        variable_sh_ACC{icell,iran}(:,8)=[[0;0];spike_hist_ACC_reward{icell}(ran_y_ACC{icell,iran}(1,1:end-2),1);[0]];
        variable_sh_ACC{icell,iran}(:,9)=[[0;0;0];spike_hist_ACC_reward{icell}(ran_y_ACC{icell,iran}(1,1:end-3),1);[0]];
        [coeff2, dev2, stats2] = glmfit(variable_sh_ACC{icell,iran}(1:end-1,:), spike_hist_ACC_reward{icell}, 'normal');
        P_value_sh_ACC{icell,iran}(1:6,1)=stats2.p(2:7,:);
        std_x_sh = std(variable_sh_ACC{icell,iran}(1:end-1,:));
        SRC_sh_ACC{icell,iran}(1:6,1) = coeff2(2:7,1).*std_x_sh(1:6)'./std_y;
                
    end
end

%% Saving
cd(saving_path)
save(['totalCRX_ACC_yj_reward_' num2str(win_title) '.mat'],...
    'SRC_ACC','P_value_ACC','spike_hist_ACC_reward','variable_ACC','ACC_PYR');
save(['totalCRX_ACC_yj_y_' num2str(iter) '_shuffle_reward_' num2str(win_title) '.mat'],...
    'P_value_sh_ACC', 'SRC_sh_ACC','ACC_PYR', 'spike_hist_ACC_reward','variable_sh_ACC','ran_y_ACC');
disp('ACC done')
toc

%% PRL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Loading Beh(for reward) & Event
spk_data_PRL=cell(length(PRL_PYR),1);  variable_PRL=cell(length(PRL_PYR),1);  spike_hist_PRL_reward={};
variable_sh_PRL=cell(length(PRL_PYR),iter); ran_y_PRL = cell(length(PRL_PYR),iter);
for icell=1:length(PRL_PYR)
    [path, name]=fileparts(PRL_PYR{icell});
    cd(path);
    load('Event_AC_YJ_1220.mat');
    eventtime =[];
    eventtime(1:5,:) = Event_Time(1:5,:);
    ntrial=size(eventtime,2);
    valid_trial = (1:ntrial);
    %% Variable 1:L, 2:R, 3:Qc (delta_Q)  4,5,6: CRX
    load('beh_.mat');
    ran_trial = size(beh_{1,2},1);
    variable_PRL{icell,1}(:,1:3)=beh_{1,2}([valid_trial ran_trial],1:3);
    load('new_V_value_yj.mat');
    variable_PRL{icell,1}(:,4:6)=value_v_final([valid_trial ran_trial],1:3);
    
    %% shuffle Y
    for iran = 1:iter
        ran_y = randperm(ntrial);
        ran_y_PRL{icell,iran}=ran_y;
        variable_sh_PRL{icell,iran}(:,4:6)=value_v_final([valid_trial(ran_y) ran_trial],1:3);
        variable_sh_PRL{icell,iran}(:,1:3)=beh_{1,2}([valid_trial(ran_y) ran_trial],1:3);
    end
   %% Loading Spike datad
    name=[name '.t'];
    spk_data_PRL{icell} = get_Spiketime(name)*100;
    %% Rearrange data
    for itrial=1:ntrial
        ievent = 5;
        [spk_number, time_idx] = histc(spk_data_PRL{icell,1},eventtime(ievent,itrial)+win(1,:));
        spike_hist_PRL_reward{icell,1}(itrial,1)=spk_number(1);
    end
end
%% checking error
disp('PRL_PYR_checking error....');
for icell=1:length(PRL_PYR)
    itype=5;
        if length(spike_hist_PRL_reward{icell}(:,1))+1 == length(variable_PRL{icell}(:,1));
        else
            delta= length(spike_hist_PRL_reward{icell}(:,1))+1 - length(variable_PRL{icell}(:,1));
            spike_hist_PRL_reward{icell}= spike_hist_PRL_reward{icell}(:,1:end-delta);
            fprintf('Delta of %d in cell %d',delta,icell);
        end
end

%% GLM fitting
disp('PRL_PYR_GLM fitting....');
P_value_PRL=cell(length(PRL_PYR),1); P_value_sh_PRL=cell(length(PRL_PYR),iter);
SRC_PRL=cell(length(PRL_PYR),1); SRC_sh_PRL=cell(length(PRL_PYR),iter);
for icell=1:length(PRL_PYR)
    %% for each neuron
    fprintf('Fitting cell no %d\n',icell)
    std_x = std(variable_PRL{icell}(1:end-1,:));
    
    variable_PRL{icell}(:,7)=[[0];spike_hist_PRL_reward{icell}(1:end-1,1);[0]];
    variable_PRL{icell}(:,8)=[[0;0];spike_hist_PRL_reward{icell}(1:end-2,1);[0]];
    variable_PRL{icell}(:,9)=[[0;0;0];spike_hist_PRL_reward{icell}(1:end-3,1);[0]];
    [coeff, dev, stats] = glmfit(variable_PRL{icell}(1:end-1,:), spike_hist_PRL_reward{icell}, 'normal');
    P_value_PRL{icell,1}(1:6,1) = stats.p(2:7,:);
    std_y = std(spike_hist_PRL_reward{icell});
    SRC_PRL{icell,1}(1:6,1) = coeff(2:7,1).*std_x(1:6)'./std_y;

    for iran=1:iter
        variable_sh_PRL{icell,iran}(:,7)=[[0];spike_hist_PRL_reward{icell}(ran_y_PRL{icell,iran}(1,1:end-1),1);[0]];
        variable_sh_PRL{icell,iran}(:,8)=[[0;0];spike_hist_PRL_reward{icell}(ran_y_PRL{icell,iran}(1,1:end-2),1);[0]];
        variable_sh_PRL{icell,iran}(:,9)=[[0;0;0];spike_hist_PRL_reward{icell}(ran_y_PRL{icell,iran}(1,1:end-3),1);[0]];
        [coeff2, dev2, stats2] = glmfit(variable_sh_PRL{icell,iran}(1:end-1,:), spike_hist_PRL_reward{icell}, 'normal');
        P_value_sh_PRL{icell,iran}(1:6,1)=stats2.p(2:7,:);
        std_x_sh = std(variable_sh_PRL{icell,iran}(1:end-1,:));
        SRC_sh_PRL{icell,iran}(1:6,1) = coeff2(2:7,1).*std_x_sh(1:6)'./std_y;        
    end
end

%% Saving
cd(saving_path)
save(['totalCRX_PRL_yj_reward_' num2str(win_title) '.mat'],...
    'SRC_PRL','P_value_PRL','spike_hist_PRL_reward','variable_PRL','PRL_PYR');
save(['totalCRX_PRL_yj_y_' num2str(iter) '_shuffle_reward_' num2str(win_title) '.mat'],...
    'P_value_sh_PRL', 'SRC_sh_PRL','PRL_PYR', 'spike_hist_PRL_reward','variable_sh_PRL','ran_y_PRL');
disp('PRL done')
toc

%% IL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Loading Beh(for reward) & Event
spk_data_IL=cell(length(IL_PYR),1);  variable_IL=cell(length(IL_PYR),1);  spike_hist_IL_reward={};
variable_sh_IL=cell(length(IL_PYR),iter); ran_y_IL = cell(length(IL_PYR),iter);
for icell=1:length(IL_PYR)
    [path, name]=fileparts(IL_PYR{icell});
    cd(path);
    load('Event_AC_YJ_1220.mat');
    eventtime =[];
    eventtime(1:5,:) = Event_Time(1:5,:);
    ntrial=size(eventtime,2);
    valid_trial = (1:ntrial);
    %% Variable 1:L, 2:R, 3:Qc (delta_Q)  4,5,6: CRX
    load('beh_.mat');
    ran_trial = size(beh_{1,2},1);
    variable_IL{icell,1}(:,1:3)=beh_{1,2}([valid_trial ran_trial],1:3);
    load('new_V_value_yj.mat');
    variable_IL{icell,1}(:,4:6)=value_v_final([valid_trial ran_trial],1:3);
    
    %% shuffle Y
    for iran = 1:iter
        ran_y = randperm(ntrial);
        ran_y_IL{icell,iran}=ran_y;
        variable_sh_IL{icell,iran}(:,4:6)=value_v_final([valid_trial(ran_y) ran_trial],1:3);
        variable_sh_IL{icell,iran}(:,1:3)=beh_{1,2}([valid_trial(ran_y) ran_trial],1:3);
    end
   %% Loading Spike datad
    name=[name '.t'];
    spk_data_IL{icell} = get_Spiketime(name)*100;
    %% Rearrange data
    for itrial=1:ntrial
        ievent = 5;
        [spk_number, time_idx] = histc(spk_data_IL{icell,1},eventtime(ievent,itrial)+win(1,:));
        spike_hist_IL_reward{icell,1}(itrial,1)=spk_number(1);
    end
end
%% checking error
disp('IL_PYR_checking error....');
for icell=1:length(IL_PYR)
    itype=5;
        if length(spike_hist_IL_reward{icell}(:,1))+1 == length(variable_IL{icell}(:,1));
        else
            delta= length(spike_hist_IL_reward{icell}(:,1))+1 - length(variable_IL{icell}(:,1));
            spike_hist_IL_reward{icell}= spike_hist_IL_reward{icell}(:,1:end-delta);
            fprintf('Delta of %d in cell %d',delta,icell);
        end
end

%% GLM fitting
disp('IL_PYR_GLM fitting....');
P_value_IL=cell(length(IL_PYR),1); P_value_sh_IL=cell(length(IL_PYR),iter);
SRC_IL=cell(length(IL_PYR),1); SRC_sh_IL=cell(length(IL_PYR),iter);
for icell=1:length(IL_PYR)
    %% for each neuron
    fprintf('Fitting cell no %d\n',icell)
    std_x = std(variable_IL{icell}(1:end-1,:));
    
    variable_IL{icell}(:,7)=[[0];spike_hist_IL_reward{icell}(1:end-1,1);[0]];
    variable_IL{icell}(:,8)=[[0;0];spike_hist_IL_reward{icell}(1:end-2,1);[0]];
    variable_IL{icell}(:,9)=[[0;0;0];spike_hist_IL_reward{icell}(1:end-3,1);[0]];
    [coeff, dev, stats] = glmfit(variable_IL{icell}(1:end-1,:), spike_hist_IL_reward{icell}, 'normal');
    P_value_IL{icell,1}(1:6,1) = stats.p(2:7,:);
    std_y = std(spike_hist_IL_reward{icell});
    SRC_IL{icell,1}(1:6,1) = coeff(2:7,1).*std_x(1:6)'./std_y;

    for iran=1:iter
        variable_sh_IL{icell,iran}(:,7)=[[0];spike_hist_IL_reward{icell}(ran_y_IL{icell,iran}(1,1:end-1),1);[0]];
        variable_sh_IL{icell,iran}(:,8)=[[0;0];spike_hist_IL_reward{icell}(ran_y_IL{icell,iran}(1,1:end-2),1);[0]];
        variable_sh_IL{icell,iran}(:,9)=[[0;0;0];spike_hist_IL_reward{icell}(ran_y_IL{icell,iran}(1,1:end-3),1);[0]];
        [coeff2, dev2, stats2] = glmfit(variable_sh_IL{icell,iran}(1:end-1,:), spike_hist_IL_reward{icell}, 'normal');
        P_value_sh_IL{icell,iran}(1:6,1)=stats2.p(2:7,:);
        std_x_sh = std(variable_sh_IL{icell,iran}(1:end-1,:));
        SRC_sh_IL{icell,iran}(1:6,1) = coeff2(2:7,1).*std_x_sh(1:6)'./std_y;

    end
end

%% Saving
cd(saving_path)
save(['totalCRX_IL_yj_reward_' num2str(win_title) '.mat'],...
    'SRC_IL','P_value_IL','spike_hist_IL_reward','variable_IL','IL_PYR');
save(['totalCRX_IL_yj_y_' num2str(iter) '_shuffle_reward_' num2str(win_title) '.mat'],...
    'P_value_sh_IL', 'SRC_sh_IL','IL_PYR', 'spike_hist_IL_reward','variable_sh_IL','ran_y_IL');
disp('IL done')
toc
