function [xpar, like, exitflag, output] = fitq_Qeach(beh_dat)  %(beh_dat)

% fitq.m
% 3-parameter RL vs. Q learning
% April 14, 2007
% Daeyeol Lee

maxit=20000;
maxeval=20000;
maxlike=1; % number of searches

op=optimset('fminsearch');
op.MaxIter=maxit;
op.MaxFunEvals=maxeval;      

mlik=-1;        
for klike=1:maxlike,
    ipar=rand(1,3); %(1,3)
    [xpar like exitflag output]=fminsearch(@QLfun_Q_each, ipar, op, beh_dat);               
%     if mlik<0  || like < mlik
%         mlik=like; 
%     end;
end;