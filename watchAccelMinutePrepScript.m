%% minute-wise ACCEL data prep
% 

minute_sampling_rate = 1/60;
accel_sampling_rate = 62.5;

%load(fullfile(destFolder,'prepDays_accel.mat'));
varnames = whos('dep*_accel');

saveAsFile = fullfile(destFolder,'prepDays_accel_minute.mat');
save(saveAsFile, 'depID','minute_sampling_rate','accel_sampling_rate','-v7.3');

for v = 1:length(varnames)
    
    eval(['accel_table' '=' varnames(v).name ';']);

    if height(accel_table)>0
        
        
        % Time Synchronization
        % Aligning timestamps from watch to sampling rate based time slots
        % Finding missing data packets

        hour_range = accel_table.timeIndex(1).Hour:accel_table.timeIndex(end).Hour;
        minute_range = 0:59;
        second_range = 0:1/minute_sampling_rate:60-(1/minute_sampling_rate);

        table_length = length(second_range)*length(minute_range)*length(hour_range);

        new_accel_timeIndex = accel_table.timeIndex(1:accel_sampling_rate*60:end);
        new_accel_mean = -1*ones(table_length,1);
        new_accel_var = -1*ones(table_length,1);
        new_accel_teager = -1*ones(table_length,1);
        
        for k = 0:table_length-1 %1:accel_sampling_rate*60:height(accel_table)
            kaccx = accel_table.accx(1+(k*accel_sampling_rate*60):min((k+1)*accel_sampling_rate*60,height(accel_table)));
            kaccy = accel_table.accy(1+(k*accel_sampling_rate*60):min((k+1)*accel_sampling_rate*60,height(accel_table)));
            kaccz = accel_table.accz(1+(k*accel_sampling_rate*60):min((k+1)*accel_sampling_rate*60,height(accel_table)));
            if (all(kaccx>-30) && all(kaccy>-30) && all(kaccz>-30))
                kmag = sqrt((kaccx.*kaccx)+(kaccy.*kaccy)+(kaccz.*kaccz));
                new_accel_mean(k+1) = mean(kmag);
                new_accel_var(k+1) = var(kmag);
                new_accel_teager(k+1) = mean(rTeagerCompute(kmag,2)); 
            end
        end
        
        
        
        % new_agi_* contain the synced data
        
        clear accel_table khour kminute ksecond hour_range minute_range second_range table_length;
        clear k kaccx kaccy kaccz kmag;
        
        % Filter-out noisy clipped data
        % Median filtering for reducing speckle noise
        % Impute missing data
        % Imputed with the local mean value

        
        ind_accel = find(new_accel_teager>=-1);
        % find the data points with atleast one-minute gaps in between
        indind_accel = find(ind_accel(2:end)-ind_accel(1:end-1)>1 & ...
            ind_accel(2:end)-ind_accel(1:end-1)<minute_sampling_rate*60*3)+1;
        zind_accel_1 = ind_accel(indind_accel);
        zind_accel_2 = ind_accel(indind_accel-1);
        for k = 1:length(zind_accel_1)
            new_accel_teager(zind_accel_2(k)+1:zind_accel_1(k)-1) = ...
                mean([new_accel_teager(zind_accel_2(k)),new_accel_teager(zind_accel_1(k))]);
            new_accel_mean(zind_accel_2(k)+1:zind_accel_1(k)-1) = ...
                mean([new_accel_mean(zind_accel_2(k)),new_accel_mean(zind_accel_1(k))]);
            new_accel_var(zind_accel_2(k)+1:zind_accel_1(k)-1) = ...
                mean([new_accel_var(zind_accel_2(k)),new_accel_var(zind_accel_1(k))]);
        end
                


        clear ind_accel indind_accel;
        clear zind_accel_1 zind_accel_2;

        % save and store pre-processed data
        processedtable = table;
        processedtable.timeIndex = new_accel_timeIndex; 
        processedtable.teager = new_accel_teager;
        processedtable.mean = new_accel_mean;
        processedtable.var = new_accel_var;
    else
        processedtable = table;
        processedtable.timeIndex = '9/17/1970'; 
        processedtable.teager = -1;
        processedtable.mean = -1;
        processedtable.var = -1;
        processedtable(:,:) = [];
    end

    dateNumber = split(varnames(v).name,{'dep','_w','_d','_accel'});
    tname = ['dep' num2str(depID) '_w' dateNumber{3} '_d' dateNumber{4} '_accel_minute'];
    eval([ tname '=' 'processedtable;']);
    clear processedtable dateNumber;
    
    save(saveAsFile, tname,'-append');
    %eval(['clear ' tname ';']);
    eval(['clear ' varnames(v).name ';']);
    clear tname new_accel_* ;
    
    
end

clear v k varnames saveAsFile;