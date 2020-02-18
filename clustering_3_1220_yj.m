close all; clear all;
%%
load('E:\Data\Neural analysis\invalid_sessions_1220.mat')
cd('E:\Data\CRX\MeanFR\reward01');
region = 'IL';
load([region '_cell.mat'], 'celldata');
H=length(celldata);
type={}; PYR=[]; INT=[]; TAN=[]; UIN=[]; PYR_t=[]; INT_t=[]; TAN_t=[]; UIN_t=[];
%%
for i=1:H;
   %% old standard classification
   [path, name]=fileparts(char(celldata(i,1)));
   ses_name = strsplit(path,'\');
%    if isequal(ses_name{end},'2016-03-30_14-02-24') || isequal(ses_name{end},'2016-09-05_14-20-23')
   if ismember(path,invalid_sessions_1220)
       continue
   end
   if cell2mat(celldata(i,2)) >0.24 && cell2mat(celldata(i,3)) < 8.83
       PYR=celldata(i,:);
       PYR_t=[PYR_t;PYR];

   else
       INT=celldata(i,:);
       INT_t=[INT_t;INT];
       
    end
end
 type = {PYR_t, INT_t, TAN_t, UIN_t};
  
save(['celltype_' region '_YJ_1220.mat'], 'type');

