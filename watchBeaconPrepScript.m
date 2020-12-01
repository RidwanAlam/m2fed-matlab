%% BEACON Data Prep
%

varnames = whos('w*_d*_beacon');
beacon_sampling_rate = 1/60; % every 120 seconds

saveAsFile = fullfile(destFolder,'prepDays_beacon.mat');
save(saveAsFile, 'depID','beacon_sampling_rate','-v7.3');

for i = 1:length(varnames)
    
    eval(['beacon_table' '=' varnames(i).name ';']);

    % Magnitude Adjustment  
    % Clipping out-of-range values
    if height(beacon_table)>0
        beacon_table.rssi = max(min(beacon_table.rssi,-10),-110);
        [~,beacon_table.room] = ismember(beacon_table.beaconID, beacons);
        
        
        % Time Synchronization
        % Aligning timestamps from watch to sampling rate based time slots
        % Finding missing data packets

        hour_range = beacon_table.timeIndex(1).Hour:beacon_table.timeIndex(end).Hour;
        minute_range = 0:59;
        second_range = 0:1/beacon_sampling_rate:60-(1/beacon_sampling_rate);
        d1year = beacon_table.timeIndex(1).Year;
        d1month = beacon_table.timeIndex(1).Month;
        d1day = beacon_table.timeIndex(1).Day;

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
        new_beacon_timeIndex = datetime([d1year*ones(table_length,1),d1month*ones(table_length,1),...
            d1day*ones(table_length,1), d1hour, d1minute, d1second]);
        for r = 1:length(beacons)
            eval(['new_beacon_room_' num2str(r) '=' '-120*ones(table_length,1);']);
        end
        %new_beacon_room = -120*ones(table_length,1);
        
        
        khour = beacon_table.timeIndex.Hour - new_beacon_timeIndex(1).Hour;
        kminute = beacon_table.timeIndex.Minute;
        ksecond = round(max(beacon_table.timeIndex.Second-...
            (1/(2*beacon_sampling_rate)),0)*beacon_sampling_rate);
        
        kindex = khour*hlen + kminute*mlen + ksecond +1;
        [kunique,kindexi,~]=unique(kindex);
        
        for k = 1:length(kunique)
            kindexk = kunique(k);
            if k<length(kunique)
                krooms = beacon_table.room(kindexi(k):kindexi(k+1)-1);
                krssi = beacon_table.rssi(kindexi(k):kindexi(k+1)-1);
            else 
                krooms = beacon_table.room(kindexi(k):length(kindex));
                krssi = beacon_table.rssi(kindexi(k):length(kindex));
            end
            for r = 1:length(beacons)
                rrssi = krssi(krooms==r);
                if ~isempty(rrssi)
                    rssi = median(rrssi);
                else
                    rssi = -120;
                end
                eval(['new_beacon_room_' num2str(r) '(' num2str(kindexk) ')' ...
                    '=' num2str(rssi) ';']);
            end
        end
        
        
        clear beacon_table khour kminute ksecond kindex kunique kindexi hour_range minute_range second_range;
        clear d1year d1month d1day d1hour d1minute d1second hlen mlen;
        clear rssi rrssi r k kindexk krooms krssi;
        
        % Filter-out noisy clipped data
        % Median filtering for reducing speckle noise
        


        % Impute missing data
        % Imputed with the local mean value

        for r =1:length(beacons)
            new_beacon_room = eval(['new_beacon_room_' num2str(r)]);
            ind_beacon = find(new_beacon_room>=-110);
            % find the data points with atleast one-minute gaps in between
            indind_beacon = find(ind_beacon(2:end)-ind_beacon(1:end-1)>1 & ...
                ind_beacon(2:end)-ind_beacon(1:end-1)<beacon_sampling_rate*60*2.5)+1;
            zind_beacon_1 = ind_beacon(indind_beacon);
            zind_beacon_2 = ind_beacon(indind_beacon-1);
            for k = 1:length(zind_beacon_1)
                new_beacon_room(zind_beacon_2(k)+1:zind_beacon_1(k)-1) = ...
                    new_beacon_room(zind_beacon_2(k));
                    %mean([new_beacon_room(zind_beacon_2(k)),new_beacon_room(zind_beacon_1(k))]);
            end
            eval(['new_beacon_room_' num2str(r) '=' 'new_beacon_room;']);
        end


        clear ind_beacon indind_beacon;
        clear zind_beacon zind_beacon_1 zind_beacon_2;

        % save and store pre-processed data
        % as p2dS_dDD
        processedtable = table;
        processedtable.timeIndex = new_beacon_timeIndex; 
        for r = 1:length(beacons)
            eval(['processedtable.room' num2str(r) '=' 'new_beacon_room_' num2str(r) ';']);
        end
    else
        processedtable = table;
        processedtable.timeIndex = '9/17/1970'; 
        for r = 1:length(beacons)
            eval(['processedtable.room' num2str(r) '=' '-120;']);
        end
        processedtable(:,:) = [];
    end

    dateNumber = split(varnames(i).name,{'w','_d','_beacon'});
    tname = ['dep' num2str(depID) '_w' dateNumber{2} '_d' dateNumber{3} '_beacon'];
    eval([ tname '=' 'processedtable;']);
    clear processedtable dateNumber;
    
    save(saveAsFile, tname,'-append');
    eval(['clear ' tname ';']);
    clear tname new_beacon_* ;
    
end

clear i saveAsFile varnames;