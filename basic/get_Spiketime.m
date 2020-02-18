function SpikeTime = get_Spiketime(FileName)
%% Get a spiketime file data, which looks like "?_t?".
%% The file is a binary data, consisted of time data from Cheetah Acquisition System.


Record_size = 4; %% One spike has 4 byte capacity

% disp(['Load..',FileName]); 

fid = fopen(FileName,'r','b');
fseek(fid, 0, 'eof');
File_size = ftell(fid);

fseek(fid,0,'bof');
while 1
    header = fgetl(fid);
    
 if (strcmp(header,'%%ENDHEADER'))
        break;
 end
end

% pause;

Header_size = ftell(fid);
N_records = (File_size - Header_size)/Record_size;
SpikeTime = zeros(N_records,1);
 
for i=1:N_records
    SpikeTime(i) = fread(fid, 1,'ulong');
end

%% SpikeTime is recorded per 100 microsecond

fclose(fid);
