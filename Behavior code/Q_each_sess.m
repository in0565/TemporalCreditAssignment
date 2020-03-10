clc; clear all; close all;

%% Load beh.file
cd('E:\systemsneuro\recording\DLS\recording_Y2')
    behfiles = FindFiles('beh_.mat');    
%% for total session
  tot_t=0; sess={}; beh_dat={};
for ifile=1:size(behfiles,1);
    sess=struct2cell(load(behfiles{ifile}));
    beh_dat{1}(:,1)=(sess{1}{2}(:,1)+1)/2; % CR: -1,1 -->0,1 
    beh_dat{1}(:,2)=(sess{1}{2}(:,2)+1)/2;
    beh_dat{1}(:,3)=sess{1}{2}(:,5); %4:RT(converse time~reward) 5:log(RT)
  tot_t=(size(beh_dat{1},1)-3); %3trials Á¦¿Ü 
  %% total alpha, beta
total=[];  %total_sess={};
    for i=1:1000
        [xpar_t, like_t, exitflag_t, output_t] = fitq_Qeach(beh_dat);
        total(i,1)=like_t;
        total(i,2:4)=xpar_t; % beta, alpha, gamma
    end
    find(total(:,1)==min(total(:,1)));
    min_tot=total(find(total(:,1)==min(total(:,1))),:);

    mlik=min(total(:,1)); 

    AIC{ifile,1}=2*mlik + 2*(3);
    BIC{ifile,1}=2*mlik + (3).* log(tot_t);
    BIC{ifile,2}=tot_t;

    min_total(ifile,:)=min_tot(1,:);
    clear sess, clear beh_dat,
end
cd('E:\thesis\Results\alpha_beta_gamma\RL');
save(['Y2_like.mat'], 'min_total', 'AIC', 'BIC');

clc; clear all; close all;
