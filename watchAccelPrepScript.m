%% ACCEL data prep
% 

varnames = whos('w*_d*_accel');
accel_sampling_rate = 62.5;

saveAsFile = fullfile(destFolder,'prepDays_accel.mat');
save(saveAsFile, 'depID','accel_sampling_rate','-v7.3');

for i = 1:length(varnames)
    
    eval(['acc_table' '=' varnames(i).name ';']);

    % Magnitude Adjustment  
    % Clipping out-of-range values
    if height(acc_table)>0
        acc_table.x(abs(acc_table.x)>29.99) = sign(acc_table.x(abs(acc_table.x)>29.99))*30;
        acc_table.y(abs(acc_table.y)>29.99) = sign(acc_table.y(abs(acc_table.y)>29.99))*30;
        acc_table.z(abs(acc_table.z)>29.99) = sign(acc_table.z(abs(acc_table.z)>29))*30;

        % Time Synchronization
        % Aligning timestamps from watch to sampling rate based time slots
        % Finding missing data packets

        hour_range = acc_table.timeIndex(1).Hour:acc_table.timeIndex(end).Hour;
        minute_range = 0:59;
        second_range = 0:1/accel_sampling_rate:60-(1/accel_sampling_rate);
        d1year = acc_table.timeIndex(1).Year;
        d1month = acc_table.timeIndex(1).Month;
        d1day = acc_table.timeIndex(1).Day;

        table_length = length(second_range)*length(minute_range)*length(hour_range);
        d1hour = zeros(table_length,1);
        d1minute = zeros(table_length,1);
        d1second = zeros(table_length,1);
        hlen = (length(second_range)*length(minute_range));
        mlen = length(second_range);
        for k = 1:length(hour_range)
            d1hour((hlen*(k-1))+1:(hlen*(k-1))+hlen) = hour_range(k);
            for j = 1:length(minute_range)
                d1minute((hlen*(k-1))+(mlen*(j-1))+1:(hlen*(k-1))+(mlen*(j-1))+mlen) ...
                    = minute_range(j);
                d1second((hlen*(k-1))+(mlen*(j-1))+1:(hlen*(k-1))+(mlen*(j-1))+mlen) ...
                    = second_range;
            end
        end

        % value assigned to missing data points = -30
        new_acc_timeIndex = datetime([d1year*ones(table_length,1),d1month*ones(table_length,1),...
            d1day*ones(table_length,1), d1hour, d1minute, d1second]);
        new_acc_x = -30*ones(table_length,1);
        new_acc_y = -30*ones(table_length,1);
        new_acc_z = -30*ones(table_length,1);


        khour = acc_table.timeIndex.Hour - new_acc_timeIndex(1).Hour;
        kminute = acc_table.timeIndex.Minute;
        ksecond = round(max(acc_table.timeIndex.Second-...
            (1/(2*accel_sampling_rate)),0)*accel_sampling_rate);

        kindex = khour*hlen + kminute*mlen + ksecond +1;

        new_acc_x(kindex) = acc_table.x;
        new_acc_y(kindex) = acc_table.y;
        new_acc_z(kindex) = acc_table.z;

        % new_agi_* contain the synced data

        clear acc_table khour kminute ksecond hour_range minute_range second_range table_length;
        clear d1year d1month d1day d1hour d1minute d1second hlen mlen kindex;

        % Filter-out noisy clipped data
        % Median filtering for reducing speckle noise

        zind_accx = find(abs(new_acc_x)==30);
        for k = 1:length(zind_accx)
            if (zind_accx(k)>2 && zind_accx(k)<length(new_acc_x)-2)
                new_acc_x(zind_accx(k)) = median(new_acc_x(zind_accx(k)-2:zind_accx(k)+2));
            end
        end

        zind_accy = find(abs(new_acc_y)==30);
        for k = 1:length(zind_accy)
            if (zind_accy(k)>2 && zind_accy(k)<length(new_acc_y)-2)
                new_acc_y(zind_accy(k)) = median(new_acc_y(zind_accy(k)-2:zind_accy(k)+2));
            end
        end

        zind_accz = find(abs(new_acc_z)==30);
        for k = 1:length(zind_accz)
            if (zind_accz(k)>2 && zind_accz(k)<length(new_acc_z)-2)
                new_acc_z(zind_accz(k)) = median(new_acc_z(zind_accz(k)-2:zind_accz(k)+2));
            end
        end

        % Impute missing data
        % Imputed with the local mean value

        ind_accx = find(new_acc_x~=-30);
        % find the data points with atleast one-minute gaps in between
        indind_accx = find(ind_accx(2:end)-ind_accx(1:end-1)>1 & ...
            ind_accx(2:end)-ind_accx(1:end-1)<accel_sampling_rate*60*1)+1;
        zind_accx_1 = ind_accx(indind_accx);
        zind_accx_2 = ind_accx(indind_accx-1);
        for k = 1:length(zind_accx_1)
            new_acc_x(zind_accx_2(k)+1:zind_accx_1(k)-1) = ...
                mean([new_acc_x(zind_accx_2(k)),new_acc_x(zind_accx_1(k))]);
        end

        ind_accy = find(new_acc_y~=-30);
        indind_accy = find(ind_accy(2:end)-ind_accy(1:end-1)>1 & ...
            ind_accy(2:end)-ind_accy(1:end-1)<accel_sampling_rate*60*1)+1;
        zind_accy_1 = ind_accy(indind_accy);
        zind_accy_2 = ind_accy(indind_accy-1);
        for k = 1:length(zind_accy_1)
            new_acc_y(zind_accy_2(k)+1:zind_accy_1(k)-1) = ...
                mean([new_acc_y(zind_accy_2(k)),new_acc_y(zind_accy_1(k))]);
        end

        ind_accz = find(new_acc_z~=-30);
        indind_accz = find(ind_accz(2:end)-ind_accz(1:end-1)>1 & ...
            ind_accz(2:end)-ind_accz(1:end-1)<accel_sampling_rate*60*1)+1;
        zind_accz_1 = ind_accz(indind_accz);
        zind_accz_2 = ind_accz(indind_accz-1);
        for k = 1:length(zind_accz_1)
            new_acc_z(zind_accz_2(k)+1:zind_accz_1(k)-1) = ...
                mean([new_acc_z(zind_accz_2(k)),new_acc_z(zind_accz_1(k))]);
        end

        clear ind_accx ind_accy ind_accz indind_accx indind_accy indind_accz;
        clear zind_accx zind_accx_1 zind_accx_2 zind_accy zind_accy_1 zind_accy_2 zind_accz zind_accz_1 zind_accz_2;

        % save and store pre-processed data
        % as p2dS_dDD
        processedtable = table;
        processedtable.timeIndex = new_acc_timeIndex; 
        processedtable.accx = new_acc_x;
        processedtable.accy = new_acc_y;
        processedtable.accz = new_acc_z;
    else
        processedtable = table;
        processedtable.timeIndex = '9/17/1970'; 
        processedtable.accx = 0;
        processedtable.accy = 0;
        processedtable.accz = 0;
        processedtable(:,:) = [];
    end

    dateNumber = split(varnames(i).name,{'w','_d','_accel'});
    tname = [num2str(depID) '_w' dateNumber{2} '_d' dateNumber{3} '_accel'];
    eval([ tname '=' 'processedtable;']);
    clear processedtable dateNumber;
    
    save(saveAsFile, tname,'-append');
    eval(['clear ' tname ';']);
    clear tname new_acc_timeIndex new_acc_x new_acc_y new_acc_z;
    
end
clear i k varnames saveAsFile;