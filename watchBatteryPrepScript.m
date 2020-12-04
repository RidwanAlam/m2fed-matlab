%% BATTERY Data Prep
%

varnames = whos('w*_d*_battery');
battery_sampling_rate = 1/60; % every 120 seconds

saveAsFile = fullfile(destFolder,'prepDays_battery.mat');
save(saveAsFile, 'depID','battery_sampling_rate','-v7.3');

for i = 1:length(varnames)
    
    eval(['battery_table' '=' varnames(i).name ';']);

    % Magnitude Adjustment  
    % Clipping out-of-range values
    if height(battery_table)>0
        battery_table.power = max(min(battery_table.power,1),0);
        battery_table.charging = double(battery_table.charging~=502);
        
        
        % Time Synchronization
        % Aligning timestamps from watch to sampling rate based time slots
        % Finding missing data packets

        hour_range = battery_table.timeIndex(1).Hour:battery_table.timeIndex(end).Hour;
        minute_range = 0:59;
        second_range = 0:1/battery_sampling_rate:60-(1/battery_sampling_rate);
        d1year = battery_table.timeIndex(1).Year;
        d1month = battery_table.timeIndex(1).Month;
        d1day = battery_table.timeIndex(1).Day;

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
        
        % value assigned to missing data points = -110
        new_battery_timeIndex = datetime([d1year*ones(table_length,1),d1month*ones(table_length,1),...
            d1day*ones(table_length,1), d1hour, d1minute, d1second]);
        new_battery_charging = -1*ones(table_length,1);
        new_battery_power = -1*ones(table_length,1);
        %new_beacon_rssi = -110*ones(table_length,1);
        
        
        khour = battery_table.timeIndex.Hour - new_battery_timeIndex(1).Hour;
        kminute = battery_table.timeIndex.Minute;
        ksecond = round(max(battery_table.timeIndex.Second-...
            (1/(2*battery_sampling_rate)),0)*battery_sampling_rate);
        
        kindex = khour*hlen + kminute*mlen + ksecond +1;
        [kunique,kindexi,~]=unique(kindex);
        
        for k = 1:length(kunique)
            kindexk = kunique(k);
            if k<length(kunique)
                kpower = battery_table.power(kindexi(k):kindexi(k+1)-1);
                kcharging = battery_table.charging(kindexi(k):kindexi(k+1)-1);
            else 
                kpower = battery_table.power(kindexi(k):length(kindex));
                kcharging = battery_table.charging(kindexi(k):length(kindex));
            end
            new_battery_charging(kindexk) = double(all(kcharging));
            new_battery_power(kindexk) = mean(kpower);            
        end
        
        
        % new_agi_* contain the synced data
        
        clear battery_table khour kminute ksecond hour_range minute_range second_range;
        clear d1year d1month d1day d1hour d1minute d1second hlen mlen;
        clear kindex kunique kindexi kindexk k kpower kcharging;
        
        % Filter-out noisy clipped data
        % Median filtering for reducing speckle noise
        


        % Impute missing data
        % Imputed with the local mean value

        
        ind_battery = find(new_battery_power>=0);
        % find the data points with atleast one-minute gaps in between
        indind_battery = find(ind_battery(2:end)-ind_battery(1:end-1)>1 & ...
            ind_battery(2:end)-ind_battery(1:end-1)<battery_sampling_rate*60*5)+1;
        zind_battery_1 = ind_battery(indind_battery);
        zind_battery_2 = ind_battery(indind_battery-1);
        for k = 1:length(zind_battery_1)
            new_battery_power(zind_battery_2(k)+1:zind_battery_1(k)-1) = ...
                mean([new_battery_power(zind_battery_2(k)),new_battery_power(zind_battery_1(k))]);
            new_battery_charging(zind_battery_2(k)+1:zind_battery_1(k)-1) = ...
                mean([new_battery_charging(zind_battery_2(k)),new_battery_charging(zind_battery_1(k))]);
        end
                


        clear ind_battery indind_battery;
        clear zind_battery zind_battery_1 zind_battery_2;

        % save and store pre-processed data
        % as p2dS_dDD
        processedtable = table;
        processedtable.timeIndex = new_battery_timeIndex; 
        processedtable.power = new_battery_power;
        processedtable.charging = new_battery_charging;
    else
        processedtable = table;
        processedtable.timeIndex = '9/17/1970'; 
        processedtable.power = -1;
        processedtable.charging = -1;
        processedtable(:,:) = [];
    end

    dateNumber = split(varnames(i).name,{'w','_d','_battery'});
    tname = [num2str(depID) '_w' dateNumber{2} '_d' dateNumber{3} '_battery'];
    eval([ tname '=' 'processedtable;']);
    clear processedtable datenumber;
    
    save(saveAsFile, tname,'-append');
    eval(['clear ' tname ';']);
    clear tname new_battery_* ;
    
end

clear i saveAsFile varnames;
%%