clear; close all;

%% Load beh.file
cd('E:\Data\Neural analysis')
behfiles = FindFiles('beh_.mat');
%% for total session
sess={};  beh_dat={}; beh_rat = {};
for ifile=1:size(behfiles,1);
    sess=struct2cell(load(behfiles{ifile}));
    beh_dat{ifile,1}=sess{1}{2};
    sess_name = strsplit(sess{1}{1},'_');
    beh_rat{ifile,1}=sess_name{2};
    clear sess,
end

%% Q-value C,R -1, 1 --> 0, 1
no_day = length(beh_dat); value_t = cell(no_day,1); value_v = cell(no_day,1);
for i = 1:no_day
    %% total alpha, beta
    cd('E:\Data\DH_Results\alpha_beta_gamma_1');
    if isequal(beh_rat{i,1},'LE9')
        load('LE9_like');
    elseif isequal(beh_rat{i,1},'LE16')
        load('LE16_like');
    elseif isequal(beh_rat{i,1},'LE17')
        load('LE17_like');
    else
        disp('error in session animal name')
        break
    end
    alpha_t=min_tot(1,3); beta_t=min_tot(1,2);    gamma=min_tot(1,4);  S_=min_tot(1,5); S0_=min_tot(1,6);

    %%
    choice=(beh_dat{i}(:,1)+1)*0.5;
    reward=(beh_dat{i}(:,2)+1)*0.5;
    no_trial=length(choice);
    v_right_t = 0.5; v_left_t = 0.5; v_left_sr=0;  v_right_sr=0;
    value_t{i} = zeros(no_trial,2);      value_v{i} = zeros(no_trial,2);
    
    for j = 1:no_trial
        pright(j,1)=exp(beta_t.*(v_right_t+gamma))/(exp(beta_t.*(v_right_t+gamma))+exp(beta_t.*v_left_t));
        pleft(j,1)=1-pright(j,1);
        value_t{i}(j+1,1) = v_left_t;         value_t{i}(j+1,2) = v_right_t;
        value_v{i}(j+1,1) = v_left_t + v_left_sr;         value_v{i}(j+1,2) = v_right_t + v_right_sr;
        if choice(j,1) > 0
            v_right_t = v_right_t + alpha_t*(reward(j,1)-v_right_t);
            value_t{i}(j+1,2) = v_right_t;
            v_left_sr = 0;
            v_right_sr = S_*(reward(j)==1)+ S0_*(reward(j)==0) ;
%             v_right_sr = S_*(reward(j)==1)+ S_*(reward(j)==0) ;
        else
            v_left_t = v_left_t + alpha_t*(reward(j,1)-v_left_t);
            value_t{i}(j+1,1) = v_left_t;
            v_right_sr = 0;
            v_left_sr = S_*(reward(j)==1)+ S0_*(reward(j)==0);
%             v_left_sr = S_*(reward(j)==1)+ S_*(reward(j)==0);
        end
    end
    choice1=(choice-1)*(-1); % Lt:1, Rt:0
    %% Qc Qu dQ
    Q_c=value_t{i}(1:end-1,1).*choice1 + value_t{i}(1:end-1,2).*choice;
    value_t{i}(:,3)=[Q_c;[0]];
    UQ_c=value_t{i}(1:end-1,1).*choice + value_t{i}(1:end-1,2).*choice1;
    value_t{i}(:,4)=[UQ_c;[0]];
    value_t{i}(:,5)=abs(value_t{i}(:,1)-value_t{i}(:,2));
    
    Q_c2=value_v{i}(1:end-1,1).*choice1 + value_v{i}(1:end-1,2).*choice;
    value_v{i}(:,3)=[Q_c2;[0]];
    UQ_c2=value_v{i}(1:end-1,1).*choice + value_v{i}(1:end-1,2).*choice1;
    value_v{i}(:,4)=[UQ_c2;[0]];
    value_v{i}(:,5)=abs(value_v{i}(:,1)-value_v{i}(:,2));
    %% state value
    value_t{i}(1:end-1,6)=pleft(:,1).*value_t{i}(1:end-1,1) +pright(:,1).*value_t{i}(1:end-1,2);
    value_v{i}(1:end-1,6)=pleft(:,1).*value_v{i}(1:end-1,1) +pright(:,1).*value_v{i}(1:end-1,2);
    clear pleft, clear pright,
    %% alpha
    value_t{i}(:,7)=alpha_t;
    value_v{i}(:,7)=alpha_t;
    %% saving each recording folder
    value_LRCD=value_t{i};
    value_v_final=value_v{i};
    [path{i}, name, ext]=fileparts(behfiles{i});
    cd(path{i});
    save(['new_V_value_yj.mat'], 'value_v_final');
end
