% LE9, LE16, LE17

%% fig 1D
clear; close all;

cd('E:\Data\Behavioral analysis\LE9');
load ('E:\Data\Behavioral analysis\LE9\total_cond.mat'); %txt말고 mat
B=totalblock(:,3:3:12); B=cumsum(B,2);
load('E:\Data\Behavioral analysis\LE9\beh_total');
load(['E:\Data\DH_Results\alpha_beta_gamma_1\LE9_like']);

ibeh = 5;

choices= (beh_{ibeh,1}(:,1)+1).*0.5;
rewards= (beh_{ibeh,1}(:,2)+1).*0.5;
%% PL by RW model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HH = length(choices); value_t = cell(HH,1);

alpha_t=min_tot(1,3);
beta_t=min_tot(1,2);
gamma=min_tot(1,4);
S_=min_tot(1,5); S0_=min_tot(1,6);

ebeta=beta_t.*ones(HH+1,1);
v_right_t = 0; v_left_t = 0;  v_left_sr=0;  v_right_sr=0;%0.5&
value_t{ibeh} = zeros(HH,2);

for j = 1:HH
    pright(j,1)=exp(beta_t.*(v_right_t+gamma))/(exp(beta_t.*(v_right_t+gamma))+exp(beta_t.*v_left_t));
    pleft(j,1)=1-pright(j,1);
    value_t{ibeh}(j+1,1) = v_left_t + v_left_sr;         value_t{ibeh}(j+1,2) = v_right_t + v_right_sr;
    if choices(j,1) > 0
        v_right_t = v_right_t + alpha_t*(rewards(j,1)-v_right_t);
        value_t{ibeh}(j+1,2) = v_right_t;
        v_left_sr = 0;
        v_right_sr = S_*(rewards(j,1)==1)+ S0_*(rewards(j,1)==0) ;
    else
        v_left_t = v_left_t + alpha_t*(rewards(j,1)-v_left_t);
        value_t{ibeh}(j+1,1) = v_left_t;
        v_right_sr = 0;
        v_left_sr = S_*(rewards(j,1)==1)+ S0_*(rewards(j,1)==0);
    end
end

dQ=-diff(value_t{ibeh}, 1, 2); %PL (diff=2-1)
e_pl=[];
e_pl = 1./(1+exp(-ebeta.*dQ(1:HH+1,1)));

%     e_PL=1-e_pl;        %e_PR
%% actual PL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Left_data=[]; E_pl=[];
for i=1:HH
    if i < 5
        SS = find(choices(1:i+4,1)+1 == 1);
        PL = length(SS)/(i+4);
        Left_data(i,1)=[PL];
        
        pl = mean(e_pl(1:i+4)');
        E_pl(i,1)=[pl];
        
    elseif i >= 5 && i <= (HH-5)
        SS = find(choices(i-4:i+5,1)+1 == 1);
        PL = length(SS)/10;
        Left_data(i,1)=[PL];
        
        pl = mean(e_pl(i-4:i+5)');
        E_pl(i,1)=[pl];
        
    elseif i > (HH-5)
        SS = find(choices(i:HH,1)+1 == 1);
        PL = length(SS)/(HH-i);
        Left_data(i,1)=[PL];
        
        pl = mean(e_pl(i:HH)');
        E_pl(i,1)=[pl];
    end
end

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('PaperUnits','Centimeters','PaperPosition',[2 2 19.05 14.29]);
subplot(7,2,[1 3]);

plot(1:HH,Left_data(:,1), 'linewidth', 0.8, 'color',[0 0 0]); %
xlabel('Trial number', 'FontSize', 12);
ylabel('\itP_{\itL}', 'FontSize', 12);
xlim([0 HH]);
ylim([0 1]);
hold on;

for i = 1:HH
    if choices(i,1) == 0; 
        if rewards(i,1) == 0; 
            subplot(7,2,[1 3]); plot([i i], [1.025 1], 'color', [0.2 0.2 0.2], 'linewidth', 0.3);
        else
            subplot(7,2, [1 3]); plot([i i], [1.05 1], 'color', [0.2 0.2 0.2], 'linewidth', 0.3);
        end
    else % Rt
        if rewards(i,1) == 0; 
            subplot(7,2,[1,3]); plot([i i], [-0.025 0], 'color', [0.2 0.2 0.2], 'linewidth', 0.3);
        else
            subplot(7,2,[1,3]); plot([i i], [-0.05 0], 'color', [0.2 0.2 0.2], 'linewidth', 0.3);
        end
    end
    hold on;
end

for k=1:3
    plot([B(ibeh,k) B(ibeh,k)], [0:1], 'linewidth', 0.5, 'color',[0.2 0.2 0.2]);
end
t=totalblock;
text(1.2,1.2,[num2str(t(ibeh,1)),':',num2str(t(ibeh,2))], 'Fontsize',8);

plot(1:HH, E_pl(:,1), 'linewidth', 0.8, 'color',[0.6 0.6 0.6])
axis([0 HH+1 -0.05 1.05]);
set(gca,'XTick' , [0:40:HH+1],'FontSize', 11);
set(gca,'YTick' , [0:0.5:1],'FontSize', 11);
hold on;

%% save
cd('D:\KAIST\Grad\SNL\LDH_yj\Figures');
% print('-dpdf', ['Fig1_D_YJ_P_left_' num2str(ibeh)]);
% print('-dtiff', ['Fig1_D_YJ_P_left_' num2str(ibeh)]);

cd('D:\KAIST\Grad\SNL\LDH_yj\Figures\Final')
print(gcf, '-depsc','-painters',['Fig1_D_YJ_P_left_LE9_' num2str(ibeh) '.ai']);
