% QLfun.m : Q learning model : Daeyeol Lee

function like=QLfun_delta(xpar, beh_dat)
beta=xpar(1); 
alpha=xpar(2); 
gamma=xpar(3);
 
no_day = length(beh_dat);
like=0;
for i = 1:no_day
    

    choices = beh_dat{i,1}(:,1);
    rewards = beh_dat{i,1}(:,2);
   
    nt = length(choices);
    v_right = 0; %0.5; 
    v_left = 0; %0.5; 
    for k=4:nt,
         pright=exp(beta*(v_right+gamma))/(exp(beta*(v_right+gamma))+exp(beta.*(v_left)));
%            pright=exp(beta.*(v_right+v_right_sr+gamma_1(i)))/(exp(beta.*(v_right+v_right_sr+gamma_1(i)))+exp(beta.*(v_left+v_left_sr))); %+gamma
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
        else
            v_left=v_left+alpha.*(rewards(k)-v_left); %+q_t
        end;
    end;
end