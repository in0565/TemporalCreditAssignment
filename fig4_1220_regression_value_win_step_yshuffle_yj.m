clear all; close all;
tic

cd('E:\Data\CRX\MeanFR\reward01');
ACC_type = load('celltype_ACC_YJ_1220.mat');       ACC_PYR = ACC_type.type{1}(:,1);
PRL_type = load('celltype_PRL_YJ_1220.mat');       PRL_PYR = PRL_type.type{1}(:,1);
IL_type = load('celltype_IL_YJ_1220.mat');         IL_PYR = IL_type.type{1}(:,1);

clear ACC_type, clear PRL_type, clear IL_type,
saving_path = ['D:\KAIST\Grad\SNL\LDH_yj\mat\1220\Q'];

ran_iter = 100;
%% window& bin for Fig.
win = [-1.5 3.5;-1.5 1.5;-1.5 3.5;-2.5 2.5;-2.5 3.5]*1000000;
bin = 0.5*1000000; % unit: second
sliding = 0.05*1000000;

%% ACC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Loading Beh(for reward) & Event
spk_data_ACC=cell(length(ACC_PYR),1);  variable_ACC=cell(length(ACC_PYR),1);  spike_hist_ACC={}; tem_hist_ACC={}; tem_hist2_ACC={};
variable_sh_ACC=cell(length(ACC_PYR),ran_iter); ran_y_ACC = cell(length(ACC_PYR),ran_iter);
for icell=1:length(ACC_PYR)
    [path, name]=fileparts(ACC_PYR{icell});
    cd(path);
    load('Event_AC_YJ_1220.mat');
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
    for iran = 1:ran_iter
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
        for ievent=1:5
            [~, time_idx] = histc(spk_data_ACC{icell,1},eventtime(ievent,itrial)+win(ievent,:));
            if isempty(time_idx); continue; end;
            spk_time_ACC{icell,1}{itrial,ievent}=spk_data_ACC{icell,1}(logical(time_idx))-eventtime(ievent,itrial);
        end
    end
    for itype=1:5
        winbin{itype}=win(itype,1):sliding:win(itype,2);
        for itrial=1:ntrial
            tem_hist_ACC{icell,1}{itype,1}(:,itrial) = histc(spk_data_ACC{icell,1}, eventtime(itype, itrial)+winbin{itype});
            tem_hist2_ACC{icell,1}{itype,1}(:,itrial) = conv(tem_hist_ACC{icell,1}{itype,1}(:,itrial),ones(bin/sliding,1),'same');
            spike_hist_ACC{icell,1}{itype,1}(:,itrial) = tem_hist2_ACC{icell,1}{itype,1}(0.5/sliding*1000000+1:end-0.5/sliding*1000000,itrial);
        end
    end
end
%% checking error
disp('ACC_PYR_checking error....');
for icell=1:length(ACC_PYR)
    for itype=1:5
        if length(spk_time_ACC{icell}(:,1))+1 == length(variable_ACC{icell}(:,1));
        else
            delta= length(spk_time_ACC{icell}(:,1))+1 - length(variable_ACC{icell}(:,1));
            spike_hist_ACC{icell}{itype,1}= spike_hist_ACC{icell}{itype,1}(:,1:end-delta);
            fprintf('delta %d in %d cell %d type\n',delta,icell,itype);
        end
    end
end

%% GLM fitting
disp('ACC_PYR_GLM fitting....');

P_value_ACC=cell(length(ACC_PYR),5); P_value_sh_ACC=cell(length(ACC_PYR),5,ran_iter);
QL_sig_tf_ACC = zeros(length(ACC_PYR),1); QR_sig_tf_ACC = zeros(length(ACC_PYR),1); QC_sig_tf_ACC = zeros(length(ACC_PYR),1); 

SRC_ACC=cell(length(ACC_PYR),5); SRC_sh_ACC=cell(length(ACC_PYR),5,ran_iter);

for icell=1:length(ACC_PYR)
    %% for each neuron
    fprintf('Fitting cell no %d\n',icell)
    std_x = std(variable_ACC{icell}(1:end-1,:));
    
    % regression only memory and reward stage
    for itype=4:5
        for i=1:(((win(itype,2)-win(itype,1))/1000000-1)/(sliding/1000000)+1) %125 for 1s,200ms
            variable_ACC{icell}(:,7)=[[0];spike_hist_ACC{icell}{itype,1}(i,1:end-1)';[0]];
            variable_ACC{icell}(:,8)=[[0;0];spike_hist_ACC{icell}{itype,1}(i,1:end-2)';[0]];
            variable_ACC{icell}(:,9)=[[0;0;0];spike_hist_ACC{icell}{itype,1}(i,1:end-3)';[0]];
            [coeff, dev, stats] = glmfit(variable_ACC{icell}(1:end-1,:), spike_hist_ACC{icell}{itype,1}(i,:)', 'normal');
            % p_value for fon
            P_value_ACC{icell,itype}(1:7,i) = stats.p(1:7,:);
            % SRC
            std_y = std(spike_hist_ACC{icell}{itype,1}(i,:));
            SRC_ACC{icell,itype}(1:6,i) = coeff(2:7,1).*std_x(1:6)'./std_y;
            
            for iran=1:ran_iter
                variable_sh_ACC{icell,iran}(:,7)=[[0];spike_hist_ACC{icell}{itype,1}(i,ran_y_ACC{icell,iran}(1:end-1))';[0]];
                variable_sh_ACC{icell,iran}(:,8)=[[0;0];spike_hist_ACC{icell}{itype,1}(i,ran_y_ACC{icell,iran}(1:end-2))';[0]];
                variable_sh_ACC{icell,iran}(:,9)=[[0;0;0];spike_hist_ACC{icell}{itype,1}(i,ran_y_ACC{icell,iran}(1:end-3))';[0]];
                [coeff2, dev2, stats2] = glmfit(variable_sh_ACC{icell,iran}(1:end-1,:), spike_hist_ACC{icell}{itype,1}(i,:)', 'normal');
                % p_value for fon
                P_value_sh_ACC{icell,itype,iran}(1:7,i)=stats2.p(1:7,:);
                % SRC
                std_x_sh = std(variable_sh_ACC{icell,iran}(1:end-1,:));
                SRC_sh_ACC{icell,itype,iran}(1:6,i) = coeff2(2:7,1).*std_x_sh(1:6)'./std_y;
            end
        end
    end
    % Q in reward stage (type=5)
    QL_sig_tf_ACC(icell,1) = sum(P_value_ACC{icell,5}(5,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    QR_sig_tf_ACC(icell,1) = sum(P_value_ACC{icell,5}(6,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    QC_sig_tf_ACC(icell,1) = sum(P_value_ACC{icell,5}(7,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    toc
end
%% Saving
cd(saving_path)
save(['totalCRX_ACC_yj_1220_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'P_value_ACC','SRC_ACC','spike_hist_ACC', 'ACC_PYR', 'bin','sliding','variable_ACC');
save(['totalCRX_ACC_yj_1220_y_',num2str(ran_iter),'_shuffle_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'P_value_sh_ACC','SRC_sh_ACC','ACC_PYR', 'spike_hist_ACC', 'bin', 'sliding','ran_y_ACC');
save(['ACC_Q_1220_idx_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'QL_sig_tf_ACC','QR_sig_tf_ACC','QC_sig_tf_ACC', 'bin','sliding');
disp('ACC done')
toc

%% PRL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Loading Beh(for reward) & Event
spk_data_PRL=cell(length(PRL_PYR),1);  variable_PRL=cell(length(PRL_PYR),1);  spike_hist_PRL={}; tem_hist_PRL={}; tem_hist2_PRL={};
variable_sh_PRL=cell(length(PRL_PYR),ran_iter); ran_y_PRL = cell(length(PRL_PYR),ran_iter);
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
    for iran = 1:ran_iter
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
        for ievent=1:5
            [~, time_idx] = histc(spk_data_PRL{icell,1},eventtime(ievent,itrial)+win(ievent,:));
            if isempty(time_idx); continue; end;
            spk_time_PRL{icell,1}{itrial,ievent}=spk_data_PRL{icell,1}(logical(time_idx))-eventtime(ievent,itrial);
        end
    end
    for itype=1:5
        winbin{itype}=win(itype,1):sliding:win(itype,2);
        for itrial=1:ntrial
            tem_hist_PRL{icell,1}{itype,1}(:,itrial) = histc(spk_data_PRL{icell,1}, eventtime(itype, itrial)+winbin{itype});
            tem_hist2_PRL{icell,1}{itype,1}(:,itrial) = conv(tem_hist_PRL{icell,1}{itype,1}(:,itrial),ones(bin/sliding,1),'same');
            spike_hist_PRL{icell,1}{itype,1}(:,itrial) = tem_hist2_PRL{icell,1}{itype,1}(0.5/sliding*1000000+1:end-0.5/sliding*1000000,itrial);
        end
    end
end
%% checking error
disp('PRL_PYR_checking error....');
for icell=1:length(PRL_PYR)
    for itype=1:5
        if length(spk_time_PRL{icell}(:,1))+1 == length(variable_PRL{icell}(:,1));
        else
            delta= length(spk_time_PRL{icell}(:,1))+1 - length(variable_PRL{icell}(:,1));
            spike_hist_PRL{icell}{itype,1}= spike_hist_PRL{icell}{itype,1}(:,1:end-delta);
            fprintf('delta %d in %d cell %d type\n',delta,icell,itype);
        end
    end
end

%% GLM fitting
disp('PRL_PYR_GLM fitting....');

P_value_PRL=cell(length(PRL_PYR),5); P_value_sh_PRL=cell(length(PRL_PYR),5,ran_iter);
QL_sig_tf_PRL = zeros(length(PRL_PYR),1); QR_sig_tf_PRL = zeros(length(PRL_PYR),1); QC_sig_tf_PRL = zeros(length(PRL_PYR),1); 

SRC_PRL=cell(length(PRL_PYR),5); SRC_sh_PRL=cell(length(PRL_PYR),5,ran_iter);

CPD_QL_PRL=cell(length(PRL_PYR),5); CPD_QR_PRL=cell(length(PRL_PYR),5); CPD_QC_PRL=cell(length(PRL_PYR),5);
CPD_sh_QL_PRL=cell(length(PRL_PYR),5,ran_iter); CPD_sh_QR_PRL=cell(length(PRL_PYR),5,ran_iter); CPD_sh_QC_PRL=cell(length(PRL_PYR),5,ran_iter);

for icell=1:length(PRL_PYR)
    %% for each neuron
    fprintf('Fitting cell no %d\n',icell)
    std_x = std(variable_PRL{icell}(1:end-1,:));
    
    % regression only memory and reward stage
    for itype=4:5
        for i=1:(((win(itype,2)-win(itype,1))/1000000-1)/(sliding/1000000)+1) %125 for 1s,200ms
            variable_PRL{icell}(:,7)=[[0];spike_hist_PRL{icell}{itype,1}(i,1:end-1)';[0]];
            variable_PRL{icell}(:,8)=[[0;0];spike_hist_PRL{icell}{itype,1}(i,1:end-2)';[0]];
            variable_PRL{icell}(:,9)=[[0;0;0];spike_hist_PRL{icell}{itype,1}(i,1:end-3)';[0]];
            [coeff, dev, stats] = glmfit(variable_PRL{icell}(1:end-1,:), spike_hist_PRL{icell}{itype,1}(i,:)', 'normal');
            % p_value for fon
            P_value_PRL{icell,itype}(1:7,i) = stats.p(1:7,:);
            % SRC
            std_y = std(spike_hist_PRL{icell}{itype,1}(i,:));
            SRC_PRL{icell,itype}(1:6,i) = coeff(2:7,1).*std_x(1:6)'./std_y;
            
            for iran=1:ran_iter
                variable_sh_PRL{icell,iran}(:,7)=[[0];spike_hist_PRL{icell}{itype,1}(i,ran_y_PRL{icell,iran}(1:end-1))';[0]];
                variable_sh_PRL{icell,iran}(:,8)=[[0;0];spike_hist_PRL{icell}{itype,1}(i,ran_y_PRL{icell,iran}(1:end-2))';[0]];
                variable_sh_PRL{icell,iran}(:,9)=[[0;0;0];spike_hist_PRL{icell}{itype,1}(i,ran_y_PRL{icell,iran}(1:end-3))';[0]];
                [coeff2, dev2, stats2] = glmfit(variable_sh_PRL{icell,iran}(1:end-1,:), spike_hist_PRL{icell}{itype,1}(i,:)', 'normal');
                % p_value for fon
                P_value_sh_PRL{icell,itype,iran}(1:7,i)=stats2.p(1:7,:);
                % SRC
                std_x_sh = std(variable_sh_PRL{icell,iran}(1:end-1,:));
                SRC_sh_PRL{icell,itype,iran}(1:6,i) = coeff2(2:7,1).*std_x_sh(1:6)'./std_y;
                
            end
        end
    end
    % Q in reward stage (type=5)
    QL_sig_tf_PRL(icell,1) = sum(P_value_PRL{icell,5}(5,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    QR_sig_tf_PRL(icell,1) = sum(P_value_PRL{icell,5}(6,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    QC_sig_tf_PRL(icell,1) = sum(P_value_PRL{icell,5}(7,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    toc
end
%% Saving
cd(saving_path)
save(['totalCRX_PRL_yj_1220_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'P_value_PRL','SRC_PRL','spike_hist_PRL', 'PRL_PYR', 'bin','sliding','variable_PRL');
save(['totalCRX_PRL_yj_1220_y_',num2str(ran_iter),'_shuffle_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'P_value_sh_PRL','SRC_sh_PRL','PRL_PYR', 'spike_hist_PRL', 'bin', 'sliding','ran_y_PRL');
save(['PRL_Q_1220_idx_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'QL_sig_tf_PRL','QR_sig_tf_PRL','QC_sig_tf_PRL', 'bin','sliding');
disp('PRL done')
toc

%% IL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Loading Beh(for reward) & Event
spk_data_IL=cell(length(IL_PYR),1);  variable_IL=cell(length(IL_PYR),1);  spike_hist_IL={}; tem_hist_IL={}; tem_hist2_IL={};
variable_sh_IL=cell(length(IL_PYR),ran_iter); ran_y_IL = cell(length(IL_PYR),ran_iter);
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
    for iran = 1:ran_iter
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
        for ievent=1:5
            [~, time_idx] = histc(spk_data_IL{icell,1},eventtime(ievent,itrial)+win(ievent,:));
            if isempty(time_idx); continue; end;
            spk_time_IL{icell,1}{itrial,ievent}=spk_data_IL{icell,1}(logical(time_idx))-eventtime(ievent,itrial);
        end
    end
    for itype=1:5
        winbin{itype}=win(itype,1):sliding:win(itype,2);
        for itrial=1:ntrial
            tem_hist_IL{icell,1}{itype,1}(:,itrial) = histc(spk_data_IL{icell,1}, eventtime(itype, itrial)+winbin{itype});
            tem_hist2_IL{icell,1}{itype,1}(:,itrial) = conv(tem_hist_IL{icell,1}{itype,1}(:,itrial),ones(bin/sliding,1),'same');
            spike_hist_IL{icell,1}{itype,1}(:,itrial) = tem_hist2_IL{icell,1}{itype,1}(0.5/sliding*1000000+1:end-0.5/sliding*1000000,itrial);
        end
    end
end
%% checking error
disp('IL_PYR_checking error....');
for icell=1:length(IL_PYR)
    for itype=1:5
        if length(spk_time_IL{icell}(:,1))+1 == length(variable_IL{icell}(:,1));
        else
            delta= length(spk_time_IL{icell}(:,1))+1 - length(variable_IL{icell}(:,1));
            spike_hist_IL{icell}{itype,1}= spike_hist_IL{icell}{itype,1}(:,1:end-delta);
            fprintf('delta %d in %d cell %d type\n',delta,icell,itype);
        end
    end
end

%% GLM fitting
disp('IL_PYR_GLM fitting....');

P_value_IL=cell(length(IL_PYR),5); P_value_sh_IL=cell(length(IL_PYR),5,ran_iter);
QL_sig_tf_IL = zeros(length(IL_PYR),1); QR_sig_tf_IL = zeros(length(IL_PYR),1); QC_sig_tf_IL = zeros(length(IL_PYR),1); 

SRC_IL=cell(length(IL_PYR),5); SRC_sh_IL=cell(length(IL_PYR),5,ran_iter);

CPD_QL_IL=cell(length(IL_PYR),5); CPD_QR_IL=cell(length(IL_PYR),5); CPD_QC_IL=cell(length(IL_PYR),5);
CPD_sh_QL_IL=cell(length(IL_PYR),5,ran_iter); CPD_sh_QR_IL=cell(length(IL_PYR),5,ran_iter); CPD_sh_QC_IL=cell(length(IL_PYR),5,ran_iter);

for icell=1:length(IL_PYR)
    %% for each neuron
    fprintf('Fitting cell no %d\n',icell)
    std_x = std(variable_IL{icell}(1:end-1,:));
    
    % regression only memory and reward stage
    for itype=4:5
        for i=1:(((win(itype,2)-win(itype,1))/1000000-1)/(sliding/1000000)+1) %125 for 1s,200ms
            variable_IL{icell}(:,7)=[[0];spike_hist_IL{icell}{itype,1}(i,1:end-1)';[0]];
            variable_IL{icell}(:,8)=[[0;0];spike_hist_IL{icell}{itype,1}(i,1:end-2)';[0]];
            variable_IL{icell}(:,9)=[[0;0;0];spike_hist_IL{icell}{itype,1}(i,1:end-3)';[0]];
            [coeff, dev, stats] = glmfit(variable_IL{icell}(1:end-1,:), spike_hist_IL{icell}{itype,1}(i,:)', 'normal');
            % p_value for fon
            P_value_IL{icell,itype}(1:7,i) = stats.p(1:7,:);
            % SRC
            std_y = std(spike_hist_IL{icell}{itype,1}(i,:));
            SRC_IL{icell,itype}(1:6,i) = coeff(2:7,1).*std_x(1:6)'./std_y;
           
            for iran=1:ran_iter
                variable_sh_IL{icell,iran}(:,7)=[[0];spike_hist_IL{icell}{itype,1}(i,ran_y_IL{icell,iran}(1:end-1))';[0]];
                variable_sh_IL{icell,iran}(:,8)=[[0;0];spike_hist_IL{icell}{itype,1}(i,ran_y_IL{icell,iran}(1:end-2))';[0]];
                variable_sh_IL{icell,iran}(:,9)=[[0;0;0];spike_hist_IL{icell}{itype,1}(i,ran_y_IL{icell,iran}(1:end-3))';[0]];
                [coeff2, dev2, stats2] = glmfit(variable_sh_IL{icell,iran}(1:end-1,:), spike_hist_IL{icell}{itype,1}(i,:)', 'normal');
                % p_value for fon
                P_value_sh_IL{icell,itype,iran}(1:7,i)=stats2.p(1:7,:);
                % SRC
                std_x_sh = std(variable_sh_IL{icell,iran}(1:end-1,:));
                SRC_sh_IL{icell,itype,iran}(1:6,i) = coeff2(2:7,1).*std_x_sh(1:6)'./std_y;
                
            end
        end
    end
    % Q in reward stage (type=5)
    QL_sig_tf_IL(icell,1) = sum(P_value_IL{icell,5}(5,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    QR_sig_tf_IL(icell,1) = sum(P_value_IL{icell,5}(6,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    QC_sig_tf_IL(icell,1) = sum(P_value_IL{icell,5}(7,(2/(sliding/1000000)+1):(5/(sliding/1000000)+1))<0.05)>0;
    toc
end
%% Saving
cd(saving_path)
save(['totalCRX_IL_yj_1220_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'P_value_IL','SRC_IL','spike_hist_IL', 'IL_PYR', 'bin','sliding','variable_IL');
save(['totalCRX_IL_yj_1220_y_',num2str(ran_iter),'_shuffle_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'P_value_sh_IL','SRC_sh_IL','IL_PYR', 'spike_hist_IL', 'bin', 'sliding','ran_y_IL');
save(['IL_Q_1220_idx_bin',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'],...
    'QL_sig_tf_IL','QR_sig_tf_IL','QC_sig_tf_IL', 'bin','sliding');
disp('IL done')
toc
save(['whole_workspace_Q_1220_bin_',num2str(bin/1000000),'_sliding',num2str(sliding/1000000),'.mat'])
toc
