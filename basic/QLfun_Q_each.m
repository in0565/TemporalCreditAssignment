% QLfun.m : Q learning model : Daeyeol Lee

function like=QLfun_Q_each(xpar, beh_dat)
beta=xpar(1); 
alpha=xpar(2); 
gamma=xpar(3);
% S_=xpar(4);
% S0_=xpar(5);

% no_day = length(beh_dat);
like=0;
for i = 1:1%no_day
    
    choices = beh_dat{1}(:,1);
    rewards = beh_dat{1}(:,2);
   
    nt = length(choices);
    v_right = 0;  v_left = 0; 
    v_right_sr=0;  v_left_sr=0;
    for k=4:nt,
         pright=exp(beta*(v_right+gamma))/(exp(beta*(v_right+gamma))+exp(beta.*(v_left)));
%        pright=exp(beta.*(v_right+v_right_sr+gamma))/(exp(beta.*(v_right+v_right_sr+gamma))+exp(beta.*(v_left+v_left_sr))); %+gamma
        pleft=1-pright;
        
        if pright==0, pright=10^-10; end;
        if pleft==0, pleft=10^-10; end;            
        
        if choices(k)>0, 
            logp=log(pright);
        else logp=log(pleft);
        end;

%         if k>xtrial, loglike=loglike-logp; end;        
        like=like-logp;
        
        if choices(k)>0,
            v_right=v_right+alpha.*(rewards(k)-v_right); %+q_t
%             v_right_sr = S_*(rewards(k)==1)+ S0_*(rewards(k)==0) ; %+ S0*(rewards(k)==0) ; 
              v_left_sr = 0;
        else
            v_left=v_left+alpha.*(rewards(k)-v_left); %+q_t
%              v_left_sr = S_*(rewards(k)==1)+ S0_*(rewards(k)==0) ; %+ S0*(rewards(k)==0);  
                v_right_sr = 0; 
        end;
    end;
end