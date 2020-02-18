clc; clear; close all;

%% Load beh.file
cd('E:\Data\Behavioral analysis\LE17')
load('beh_total.mat');
%% for total session
  tot_t=0; sess={}; beh_dat={};
for ifile=1:size(beh_,1);
%     sess=struct2cell(load(behfiles{ifile}));
    beh_dat{ifile,1}(:,1)=(beh_{ifile,1}(:,1)+1)/2; % CR: -1,1 -->0,1 
    beh_dat{ifile,1}(:,2)=(beh_{ifile,1}(:,2)+1)/2;
%     beh_dat{ifile,1}(:,3)=beh_{ifile,1}(:,5); %4:RT(converse time~reward) 5:log(RT)
  tot_t=tot_t+(size(beh_dat{ifile},1)); %3trials Á¦¿Ü 
end

%% total alpha, beta
total=[];  %total_sess={};
for i=1:100
[xpar_t, like_t, exitflag_t, output_t] = fitq_delta(beh_dat);
total(i,1)=like_t;
total(i,2:4)=xpar_t; % beta, alpha, gamma
end

find(total(:,1)==min(total(:,1)))
min_tot=total(find(total(:,1)==min(total(:,1))),:);

mlik=min(total(:,1)); 

AIC=2*mlik + 2*(3)
BIC=2*mlik + (3).* log(tot_t)

cd('E:\Data\Behavioral analysis\LE17');
save(['LE17_like.mat'], 'min_tot', 'AIC', 'BIC');
