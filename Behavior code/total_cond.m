clc; clear all; close all;

%% Make Beh-file each day with block-shift
cd('E:\Data\Behavioral analysis\LE9 Cond')
totalblock = load('LE9_cond.txt'); 
cd('E:\Data\Behavioral analysis\LE9')
behfiles = FindFiles('*.txt');

beh_={}; beh_B=[]; block_tot=[];
for ifile=1:length(behfiles)
            clear beh_, clear beh_B,
    [path{ifile}, name{ifile}, ext]=fileparts(behfiles{ifile});
    cd(path{ifile});
    beh_{1,1}=name{ifile};  %%index1: beh_name
%% Load txt. files w/o error    
    beh_ori=load(behfiles{ifile}); 
    beh=beh_ori(:,3:4);
    Err=find(beh_ori(:,2)==1);
    beh(Err,:)=[]; 
%% Block num. w/o error
    block=cumsum(totalblock(ifile, 3:3:12),2);
  for ii=1:length(Err)
      block= block - double(Err(ii,1) < block);
  end
  block_=block(1,1);
  block_(1,2:4) = diff(block);
   block_tot=[block_tot; block_];
end

totalblock(:,3:3:12)=block_tot;
cd('E:\Data\Behavioral analysis\LE9')
save(['total_cond.mat'], 'totalblock') 