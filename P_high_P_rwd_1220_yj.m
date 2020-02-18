clear; close all;
%%
main_path = 'E:\Data\Behavioral analysis';
sub_path = {'\LE9','\LE16','\LE17'};
cd(main_path)

session_idx = {[2:5,8:17,19:23],[4:16,19:21,24],[6,7,10,12,15]};
%%
block_length = []; P_high = []; P_reward = []; P_high_first = []; P_high_last = []; P_reward_first = []; P_reward_last = [];
for irat = 1:3
    cd([main_path sub_path{irat}])
    load ('total_cond.mat');
    load('beh_total.mat');
    
    block=cumsum(totalblock(:,3:3:12),2);
    % mean block length
    b_t = totalblock(:,3:3:12);
    block_length = [block_length reshape(b_t(session_idx{irat},:),1,[])];
    for isession=session_idx{irat}
        % L:R, Lt: -1, Rt: 1, rewarded: 1, unrewarded: -1
        choice = beh_{isession,1}(:,1);
        reward = beh_{isession,1}(:,2);
        temp_block = block(isession,:);
        % P_high
        first_left_high = totalblock(isession,1)>0.5;
        first_left_high = first_left_high*2-1;
        high_t2=sum(choice(temp_block(1)+1:temp_block(2),:)==first_left_high);
        high_t3=sum(choice(temp_block(2)+1:temp_block(3),:)==-first_left_high);
        high_t4=sum(choice(temp_block(3)+1:length(choice),:)==first_left_high);
        high_t = (high_t2+high_t3+high_t4)/length(temp_block(1)+1:length(choice))*100;
        P_high = [P_high high_t];
        % P_rwd
        reward_t = sum(reward(temp_block(1)+1:length(choice))==1)/length(temp_block(1)+1:length(choice))*100;
        P_reward = [P_reward reward_t];
        
        % 1st 10 trials
        high_tt2=sum(choice(temp_block(1)+1:temp_block(1)+10,:)==first_left_high);
        high_tt3=sum(choice(temp_block(2)+1:temp_block(2)+10,:)==-first_left_high);
        high_tt4=sum(choice(temp_block(3)+1:temp_block(3)+10,:)==first_left_high);
        high_tt = (high_tt2+high_tt3+high_tt4)/30*100;
        P_high_first = [P_high_first high_tt];
        reward_tt = sum(reward([temp_block(1)+1:temp_block(1)+10 temp_block(2)+1:temp_block(2)+10 temp_block(3)+1:temp_block(3)+10])==1)...
            /30*100;
        P_reward_first = [P_reward_first reward_tt];
        
        % last 20 trialss
        high_ttt2=sum(choice(temp_block(2)-19:temp_block(2),:)==first_left_high);
        high_ttt3=sum(choice(temp_block(3)-19:temp_block(3),:)==-first_left_high);
        high_ttt4=sum(choice(length(choice)-19:length(choice),:)==first_left_high);
        high_ttt = (high_ttt2+high_ttt3+high_ttt4)/60*100;
        P_high_last = [P_high_last high_ttt];
        reward_ttt = sum(reward([temp_block(2)-19:temp_block(2) temp_block(3)-19:temp_block(3) length(choice)-19:length(choice)])==1)...
            /60*100;
        P_reward_last = [P_reward_last reward_ttt];
    end
end

%% result
fprintf('Trial number per block: %.3f +- %.3f (SD)\n',mean(block_length),std(block_length))
fprintf('P_high (all) : %.3f +- %.3f (SD) \n',mean(P_high),std(P_high))
fprintf('P_reward (all) : %.3f +- %.3f (SD) \n',mean(P_reward),std(P_reward))
fprintf('P_high (first 10 trials) : %.3f +- %.3f (SD) \n',mean(P_high_first),std(P_high_first))
fprintf('P_reward (first 10 trials): %.3f +- %.3f (SD) \n',mean(P_reward_first),std(P_reward_first))
fprintf('P_high (last 20 trials) : %.3f +- %.3f (SD) \n',mean(P_high_last),std(P_high_last))
fprintf('P_reward (last 20 trials): %.3f +- %.3f (SD) \n',mean(P_reward_last),std(P_reward_last))
