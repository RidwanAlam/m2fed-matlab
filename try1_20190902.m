%%

depID = 30;
d = find(cell2mat(depIDs)==depID);
watchIDs = watchIDss{d};
beacons = beaconIDs{d};

depStartDate = [depStartDates{d} ' 00:00:00']; % deployment start-date 
depEndDate = [depEndDates{d} ' 23:59:59'];


startUnixTime = posixtime(datetime(depStartDate,'TimeZone','America/Los_Angeles'))*1000;
endUnixTime = posixtime(datetime(depEndDate,'TimeZone','America/Los_Angeles'))*1000;

startDate = datetime((startUnixTime)./1000,'ConvertFrom','posixtime',...
    'TimeZone','America/Los_Angeles');
endDate = datetime((endUnixTime)./1000,'ConvertFrom','posixtime',...
    'TimeZone','America/Los_Angeles');

startDateNum = floor(datenum(startDate));

startdayindex = 0;
enddayindex = floor(datenum(endDate))-startDateNum;


dataHome = 'C:\Users\Ridwan\Documents\MATLAB\M2FED Data';
sourceFolder = fullfile(dataHome,['Processed Data\deployment_',num2str(depID),'\data_mat']);
destFolder = fullfile(dataHome,['dep_',num2str(depID)]);

saveAsFile = fullfile(destFolder,'info.mat');
save(saveAsFile, 'depID','watchIDs','beacons','depStartDate','depEndDate','-v7.3');

%%
% data available as "wWW_dDD_mode" in "prepDays_modes.mat"

load(fullfile(destFolder,'prepDays_accel.mat'));
run watchAccelMinutePrepScript.m;
%load(fullfile(destFolder,'prepDays_accel_minute.mat'));
load(fullfile(destFolder,'prepDays_beacon.mat'));
load(fullfile(destFolder,'prepDays_battery.mat'));


%%
% % minute-wise motion status
% % Minute status
% 
% minute_sampling_rate = 1/60;
% 
% saveAsFile = fullfile(destFolder,'prepDays_accel_minute.mat');
% save(saveAsFile, 'depID','minute_sampling_rate','-v7.3');
% 
% varnames = whos('dep*_accel');
% 
% for v = 1:length(varnames)
%     
%     eval(['accel_table' '=' varnames(v).name ';']);
% 
%     if height(accel_table)>0
%         
%         
%         % Time Synchronization
%         % Aligning timestamps from watch to sampling rate based time slots
%         % Finding missing data packets
% 
%         hour_range = accel_table.timeIndex(1).Hour:accel_table.timeIndex(end).Hour;
%         minute_range = 0:59;
%         second_range = 0:1/minute_sampling_rate:60-(1/minute_sampling_rate);
% 
%         table_length = length(second_range)*length(minute_range)*length(hour_range);
% 
%         new_accel_timeIndex = accel_table.timeIndex(1:accel_sampling_rate*60:end);
%         new_accel_mean = -1*ones(table_length,1);
%         new_accel_var = -1*ones(table_length,1);
%         new_accel_teager = -1*ones(table_length,1);
%         
%         for k = 0:table_length-1 %1:accel_sampling_rate*60:height(accel_table)
%             kaccx = accel_table.accx(1+(k*accel_sampling_rate*60):min((k+1)*accel_sampling_rate*60,height(accel_table)));
%             kaccy = accel_table.accy(1+(k*accel_sampling_rate*60):min((k+1)*accel_sampling_rate*60,height(accel_table)));
%             kaccz = accel_table.accz(1+(k*accel_sampling_rate*60):min((k+1)*accel_sampling_rate*60,height(accel_table)));
%             if (all(kaccx>-30) && all(kaccy>-30) && all(kaccz>-30))
%                 kmag = sqrt((kaccx.*kaccx)+(kaccy.*kaccy)+(kaccz.*kaccz));
%                 new_accel_mean(k+1) = mean(kmag);
%                 new_accel_var(k+1) = var(kmag);
%                 new_accel_teager(k+1) = mean(rTeagerCompute(kmag,2)); 
%             end
%         end
%         
%         
%         
%         % new_agi_* contain the synced data
%         
%         clear accel_table khour kminute ksecond hour_range minute_range second_range table_length;
%         clear k ksccx kaccy kaccz kmag;
%         
%         % Filter-out noisy clipped data
%         % Median filtering for reducing speckle noise
%         % Impute missing data
%         % Imputed with the local mean value
% 
%         
%         ind_accel = find(new_accel_teager>=-1);
%         % find the data points with atleast one-minute gaps in between
%         indind_accel = find(ind_accel(2:end)-ind_accel(1:end-1)>1 & ...
%             ind_accel(2:end)-ind_accel(1:end-1)<minute_sampling_rate*60*3)+1;
%         zind_accel_1 = ind_accel(indind_accel);
%         zind_accel_2 = ind_accel(indind_accel-1);
%         for k = 1:length(zind_accel_1)
%             new_accel_teager(zind_accel_2(k)+1:zind_accel_1(k)-1) = ...
%                 mean([new_accel_teager(zind_accel_2(k)),new_accel_teager(zind_accel_1(k))]);
%             new_accel_mean(zind_accel_2(k)+1:zind_accel_1(k)-1) = ...
%                 mean([new_accel_mean(zind_accel_2(k)),new_accel_mean(zind_accel_1(k))]);
%             new_accel_var(zind_accel_2(k)+1:zind_accel_1(k)-1) = ...
%                 mean([new_accel_var(zind_accel_2(k)),new_accel_var(zind_accel_1(k))]);
%         end
%                 
% 
% 
%         clear ind_accel indind_accel;
%         clear zind_accel_1 zind_accel_2;
% 
%         % save and store pre-processed data
%         processedtable = table;
%         processedtable.timeIndex = new_accel_timeIndex; 
%         processedtable.teager = new_accel_teager;
%         processedtable.mean = new_accel_mean;
%         processedtable.var = new_accel_var;
%     else
%         processedtable = table;
%         processedtable.timeIndex = '9/17/1970'; 
%         processedtable.teager = -1;
%         processedtable.mean = -1;
%         processedtable.var = -1;
%         processedtable(:,:) = [];
%     end
% 
%     dateNumber = split(varnames(v).name,{'dep','_w','_d','_accel'});
%     tname = ['dep' num2str(depID) '_w' dateNumber{3} '_d' dateNumber{4} '_accel_minute'];
%     eval([ tname '=' 'processedtable;']);
%     clear processedtable;
%     
%     save(saveAsFile, tname,'-append');
%     eval(['clear ' tname ';']);
%     clear tname new_accel_* ;
%         
%     
% end
% 
% 
%%

% Save figures for each day and each watch
% Helpful for manual inspection

mmodes = {'accel_minute','battery','beacon'};

for dayindex = startdayindex:1:enddayindex
    
    f = figure(1);
    title(['Day ' num2str(dayindex)]);
    for w = 1:length(watchIDs)
        watchID = watchIDs(w);
        for m = 1:length(mmodes)
            mode = mmodes{m};
            modeDataName = ['dep' num2str(depID) '_w' num2str(watchID,'%02d') '_d' num2str(dayindex,'%02d') '_' mode];
            if exist(modeDataName,'var')==1
                modeData = eval(modeDataName);
                times = modeData.timeIndex;
                if mode=="accel_minute"
                    values = modeData.teager;
                    if exist('ax1','var')==0
                        ax1 = subplot(3,1,1);
                    end
                    plot(ax1,times,values,'Color',colors{w});hold on;
                elseif mode=="beacon"
                    [rrssi,curroom] = max(table2array(modeData(:,2:end)),[],2);
                    curroom(rrssi==-120) = -1;
                    
                    rrssi = max(rrssi+120,0.1);
                    dist = 2*abs(log(100)-log(rrssi));
                    prob = min(1,2*exp(-dist));
                    
                    values = double(prob>0.1).*curroom;                   
                    if exist('ax2','var')==0
                        ax2 = subplot(3,1,2);
                    end
                    plot(ax2,times,values,'Color',colors{w},...
                        'LineStyle','none','Marker','o','MarkerSize',2);hold on;

                elseif mode=="battery"
                    values = max(-1,(modeData.power + modeData.charging));
                    if exist('ax3','var')==0
                        ax3 = subplot(3,1,3);
                    end
                    plot(ax3,times,values,'Color',colors{w});hold on;
                end
                
            end
        end
    end
    
    filename = ['..\M2FED Data\dep_' num2str(depID,'%02d') '\figs\f' num2str(dayindex,'%02d') ];
    print([filename '.png'],'-dpng');savefig(f,[filename '.fig']);
    close(f);
    clear ax1 ax2 ax3 f;
end

%%

for dayindex = startdayindex:1:enddayindex
    
    daytime_now = startDate+dayindex;

    new_timeIndex = daytime_now + minutes(0:24*60-1)';
    daytable = table;
    daytable.timeIndex = new_timeIndex;
    
    for w = 1:length(watchIDs)
        watchID = watchIDs(w);
         
        eval(['valid_accel' '=' '-1*ones(size(new_timeIndex));']);    
        eval(['valid_beacon' '=' '-1*ones(size(new_timeIndex));']);    
        eval(['valid_battery' '=' '-1*ones(size(new_timeIndex));']);  
        eval(['valid_charging' '=' '-1*ones(size(new_timeIndex));']);
        
        mmodes = {'accel_minute','battery','beacon'};

        for m = 1:length(mmodes)
            mode = mmodes{m};
            modeDataName = ['dep' num2str(depID) '_w' num2str(watchID,'%02d') '_d' num2str(dayindex,'%02d') '_' mode];
            if exist(modeDataName,'var')==1
                
                modeData = eval(modeDataName);
                times = modeData.timeIndex;
                kindex = (times.Hour*60)+(times.Minute)+1;
                if mode=="accel_minute"
                    %valid = eval(['valid_accel_w' num2str(w,'%02d')]);
                    values = modeData.teager;
                    valid_accel(kindex(values>=0 & values<0.01)) = 0;
                    valid_accel(kindex(values>=0.01)) = 1;
                    %valid(kindex) = double(values>=0);
                    %eval(['valid_accel_w' num2str(w,'%02d') '=' 'valid;']);
                    
%                     if exist('ax1','var')==0
%                         ax1 = subplot(3,1,1);
%                     end
%                     plot(ax1,times,values,'Color',colors{w});hold on;
                elseif mode=="beacon"
                    %valid = eval(['valid_beacon_w' num2str(w,'%02d')]);
                    
                    [rrssi,curroom] = max(table2array(modeData(:,2:end)),[],2);
                    curroom(rrssi==-120) = -1;
                    
                    rrssi = max(rrssi+120,0.1);
                    dist = 2*abs(log(100)-log(rrssi));
                    prob = min(1,2*exp(-dist));
                    
                    %values = double(prob>0.1).*curroom;
                    %valid_beacon(kindex(curroom>0)) = 0;
                    valid_beacon(kindex(curroom>0)) = 1;
                    
                    
                    %valid(kindex) = double(curroom>0); 
                    %eval(['valid_beacon_w' num2str(w,'%02d') '=' 'valid;']);
%                     if exist('ax2','var')==0
%                         ax2 = subplot(3,1,2);
%                     end
%                     plot(ax2,times,values,'Color',colors{w},...
%                         'LineStyle','none','Marker','o','MarkerSize',2);hold on;

                elseif mode=="battery"
                    %valid = eval(['valid_battery_w' num2str(w,'%02d')]);
                    %values = max(-1,(modeData.power + modeData.charging));
                    valid_charging(kindex(modeData.charging==1)) = 1;
                    valid_charging(kindex(modeData.charging==0)) = 0;
                    
                    %valid_battery(kindex(modeData.charging==1)) = 0;
                    valid_battery(kindex(modeData.charging>=0)) = modeData.power(modeData.charging>=0);
                    %double(values>=0 & values<=1);                    
                    %eval(['valid_battery_w' num2str(w,'%02d') '=' 'valid;']);
                    
%                     if exist('ax3','var')==0
%                         ax3 = subplot(3,1,3);
%                     end
%                     plot(ax3,times,values,'Color',colors{w});hold on;
                end
%             else
                
            end
        end
        
        % charging: battery = accel = beacon = 0 / dead?
        valid_accel(valid_charging==1) = 0;
        valid_beacon(valid_charging==1) = 0;
        
        
        eval(['daytable.accel_w' num2str(w,'%02d') '=' 'valid_accel;']); 
        eval(['daytable.beacon_w' num2str(w,'%02d') '=' 'valid_beacon;']); 
        eval(['daytable.battery_w' num2str(w,'%02d') '=' 'valid_battery;']); 
        eval(['daytable.charging_w' num2str(w,'%02d') '=' 'valid_charging;']); 
        
        %      
        %activeStat = -1*ones(size(new_timeIndex));%eval(['active_w' num2str(w,'%02d')]);
        %unknown = 
        
        
    end
    
    
    
    eval(['d' num2str(dayindex,'%02d') '_valid_minutes' '=' 'daytable;']);
    clear daytable new_timeIndex valid_* valid values times
    
%     filename = ['..\M2FED Data\dep_' num2str(depID,'%02d') '\figs\f' num2str(dayindex,'%02d') ];
%     print([filename '.png'],'-dpng');savefig(f,[filename '.fig']);
%     close(f);
%     clear ax1 ax2 ax3 f;
end


%%
varnames = whos('d*_valid_minutes');
for v = 1:length(varnames)
    daytable = eval([varnames(v).name]);%d03_valid_minutes;
    f = figure;
    subplot(4,1,1);
    plot(daytable.timeIndex,daytable.ema_w01,'g');hold on;
    plot(daytable.timeIndex,0.7*daytable.accel_w01,'m');hold on;
    plot(daytable.timeIndex,0.4*daytable.beacon_w01,'b');
    ylabel(['Subj-1']);
    legend('Acc','Loc','EMA');
    subplot(4,1,2);
    plot(daytable.timeIndex,daytable.ema_w02,'g');hold on;
    plot(daytable.timeIndex,0.7*daytable.accel_w02,'m');hold on;
    plot(daytable.timeIndex,0.4*daytable.beacon_w02,'b');    
    ylabel(['Subj-2']);
    subplot(4,1,3);
    plot(daytable.timeIndex,daytable.ema_w03,'g');hold on;
    plot(daytable.timeIndex,0.7*daytable.accel_w03,'m');hold on;
    plot(daytable.timeIndex,0.4*daytable.beacon_w03,'b');
    ylabel(['Subj-3']);
    subplot(4,1,4);
    plot(daytable.timeIndex,daytable.ema_w04,'g');hold on;
    plot(daytable.timeIndex,0.7*daytable.accel_w04,'m');hold on;
    plot(daytable.timeIndex,0.4*daytable.beacon_w04,'b');    
    ylabel(['Subj-4']);
    filename = ['..\M2FED Data\dep_' num2str(depID,'%02d') '\figs\daily_valid_' num2str(v,'%02d') ];
    print([filename '.png'],'-dpng');savefig(f,[filename '.fig']);
    close(f);
end



%%
% interpolate device power
varnames = whos('d*_valid_minutes');
for v = 1:length(varnames)
    minutetable = eval([varnames(v).name]);

    for w = 1:length(watchIDs)

        valid_battery = eval(['minutetable.battery_w' num2str(w,'%02d')]);
        valid_beacon = eval(['minutetable.beacon_w' num2str(w,'%02d')]);
        valid_accel = eval(['minutetable.accel_w' num2str(w,'%02d')]);
        valid_charging = eval(['minutetable.charging_w' num2str(w,'%02d')]);

        clear deviceStatus;

        deviceStatus.power = zeros(size(valid_battery));
        deviceStatus.power(valid_charging==0) = valid_battery(valid_charging==0);
        deviceStatus.power(valid_charging==1) = valid_battery(valid_charging==1);
        if valid_charging(end)==-1
            endInd = find(valid_charging>=0,1,'last');
            if ~isempty(endInd)
                deviceStatus.power(end) = max(0,deviceStatus.power(endInd)-...
                    ((length(valid_battery)-endInd)/500));
            end
        end

        deviceStatus.onCharger = zeros(size(valid_battery));
        deviceStatus.onCharger(valid_charging==1) = 1;

        unknownInd = find(valid_charging==-1);
        if ~isempty(unknownInd)
            uindex = find((unknownInd(2:end)-unknownInd(1:end-1))>1);
            if isempty(uindex)
                uindex1 = unknownInd(1);
                uindex2 = unknownInd(end);
            else
                uindex1 = [unknownInd(1);unknownInd(uindex+1)];
                uindex2 = [unknownInd(uindex);unknownInd(end)];
            end

            for u=1:length(uindex1)
                initialCh = deviceStatus.power(max(1,uindex1(u)-1));
                endCh = deviceStatus.power(min(length(valid_battery),uindex2(u)+1));
                if uindex2(u)==uindex1(u)
                    deviceStatus.power(uindex1(u)) = mean([initialCh,endCh]);
                else 
                    if (endCh<=initialCh)
                        deviceStatus.power(uindex1(u):uindex2(u)) = ...
                            linspace(initialCh,endCh,uindex2(u)-uindex1(u)+1);
                    else
                        endCh = initialCh - (uindex2(u)-uindex1(u))/500;
                        deviceStatus.power(uindex1(u):uindex2(u)) = ... %initialCh/2;
                            linspace(initialCh,endCh,uindex2(u)-uindex1(u)+1);
                        deviceStatus.power(uindex1(u):uindex2(u)) = ...
                            max(0,deviceStatus.power(uindex1(u):uindex2(u)));
                    end
                end            
            end
        end

        deviceStatus.outOfCharge = zeros(size(valid_battery));    
        deviceStatus.outOfCharge(deviceStatus.power<0.01) = 1;

        deviceStatus.inSleepMode = zeros(size(valid_battery)); 
        deviceStatus.inSleepMode((deviceStatus.outOfCharge==0) && ...
            (valid_accel==-1))



        % possible dead index
    %     uindex3 = uindex1((uindex2-uindex1)>100);
    %     uindex4 = uindex2((uindex2-uindex1)>100);
    %     
    %     % other indices
    %     uindex5 = uindex1((uindex2-uindex1)<=100);
    %     uindex6 = uindex2((uindex2-uindex1)<=100);
    %     

    %     deadIndex = find((valid_battery==-1) & (valid_beacon==-1) & (valid_accel==-1));
    %     
    %     cindex = uindex1(valid_charging(uindex1-1)==1);
    %     
    %     lastKnownCharge = zeros(size(valid_battery));
    %     lastKnownCharge(valid_charging==0) = valid_battery(valid_charging==0);
    %     powerInd = find(valid_charging)
    %     for t = 1:length(chargeInd)
    %         if valid_charging(chargeInd(t)-1)
    %     lastKnownCharge(valid_charging==1) = valid_battery(valid_charging==0);
    %     lastKnownCharge
    %     lastKnownCharge(valid_battery<0) = valid
    %     
    %     deadIndex = find((valid_battery==-1) & (valid_beacon==-1) & (valid_accel==-1));
    %     
    %     
    %     
    %     lastKnownCharge = zeros(size(hour_timeIndex));
    %     lastKnownLoc = zeros(size(hour_timeIndex));
    %     
    %     for h = 1:length(hour_timeIndex)
    %         chargeArray = valid_battery(minutetable.timeIndex.Hour==hour_timeIndex(h).Hour);
    %         lastKnownCharge(h) = min(chargeArray(chargeArray>0));
    %         beaconArray = valid_beacon(minutetable.timeIndex.Hour==hour_timeIndex(h).Hour);
    %         lastKnownLoc(h) = any(beaconArray(beaconArray>0));       
    %     end
    %     
    %     deviceStats.onCharge = zeros(size(hour_timeIndex));
    %     
    %     deviceStats(valid_battery==0) = 0;
    %     
        
    
        dateNumber = split(varnames(v).name,{'d','_valid_minutes'});
        eval(['deviceStats_w' num2str(w,'%02d') '_d' dateNumber{2} '=' 'deviceStats;']);


    end

end

%%
% load EMA endtimes for each subject from excel
% as depXX_ema_sSS

% ema valid minutes

varnames = whos('d*_valid_minutes');

for v = 1:length(varnames)
    
    daytable = eval([varnames(v).name]);
    new_timeIndex = daytable.timeIndex;
    new_timeIndex.TimeZone='';
    datenum = split(varnames(v).name,{'d','_valid_'});
    
    for w = 1:length(watchIDs)
        
        ema_endtimes = eval(['dep' num2str(depID) '_ema_s' num2str(w,'%02d')]);
        offset = minutes(15);
        valid_ema = zeros(size(new_timeIndex));
        
        ema_endtimes_day = ema_endtimes((ema_endtimes>=new_timeIndex(1)-offset) & ...
            (ema_endtimes<=new_timeIndex(end)+offset));
        
        for d = 1:length(ema_endtimes_day)
            startIndex = find(new_timeIndex<=...
                max(new_timeIndex(1),ema_endtimes_day(d)-offset),1,'last');
            endIndex = find(new_timeIndex>=...
                min(new_timeIndex(end),ema_endtimes_day(d)+offset),1,'first');
            valid_ema(startIndex:endIndex) = 1;
        end
        eval(['daytable.ema_w' num2str(w,'%02d') '=' 'valid_ema;']);   
        eval(['valid_ema_d' datenum{2} '_w' num2str(w,'%02d') '=' 'valid_ema;']);
        
    end
    
    eval([varnames(v).name '=' 'daytable;']);
end

%%

% times when in-home for sure:
% subj S, source (Sensor/EMA), start time, end time

inhomeTable_Subject = [];
inhomeTable_Source = [];
inhomeTable_StartTime = [];
inhomeTable_EndTime = [];

varnames = whos('d*_valid_minutes');
for v = 1:length(varnames)
    daytable = eval([varnames(v).name]);
    datenum = split(varnames(v).name,{'d','_valid_'});
    
    for w = 1:length(watchIDs)
        % source: EMA 
        source_ema = eval(['daytable.ema_w' num2str(w,'%02d')]);
        if any(source_ema)
            ema_start_index = min(find(source_ema(2:end)-source_ema(1:end-1)>0)+1,length(source_ema));
            %ema_end_index = min(find(source_ema(2:end)-source_ema(1:end-1)<0)+1,length(source_ema));
            inhomeTable_Subject = [inhomeTable_Subject; repelem(w,length(ema_start_index))'];
            inhomeTable_Source = [inhomeTable_Source; repelem({'EMA'},length(ema_start_index))'];
            inhomeTable_StartTime = [inhomeTable_StartTime;daytable.timeIndex(ema_start_index)];
            inhomeTable_EndTime = [inhomeTable_EndTime;...
                daytable.timeIndex(min(ema_start_index+30,length(source_ema)))];
        end
        % source: Sensor
        valid_accel = eval(['daytable.accel_w' num2str(w,'%02d')]);
        valid_beacon = eval(['daytable.beacon_w' num2str(w,'%02d')]);
        source_sensor = double((valid_accel>0) & (valid_beacon>0));
        if any(source_sensor)
            sensor_start_index = min(find(source_sensor(2:end)-source_sensor(1:end-1)>0)+1,length(source_sensor));
            sensor_end_index = [];
            for i = 1:length(sensor_start_index)
                if sensor_start_index(i)==length(source_sensor)
                    sensor_end_index = [sensor_end_index; sensor_start_index(i)];
                    break;
                else
                    temp_index = find(source_sensor(sensor_start_index(i)+1:end)-...
                        source_sensor(sensor_start_index(i):end-1)<0,1,'first')+...
                        sensor_start_index(i);
                    if isempty(temp_index)
                        temp_index = length(source_sensor);
                        sensor_end_index = [sensor_end_index;temp_index];
                        break;
                    else
                        temp_index = min(temp_index,length(source_sensor));
                        sensor_end_index = [sensor_end_index;temp_index];
                    end
                end
            end
            
            inhomeTable_Subject = [inhomeTable_Subject; repelem(w,length(sensor_start_index))'];
            inhomeTable_Source = [inhomeTable_Source; repelem({'Sensor'},length(sensor_start_index))'];
            inhomeTable_StartTime = [inhomeTable_StartTime;daytable.timeIndex(sensor_start_index)];
            inhomeTable_EndTime = [inhomeTable_EndTime;daytable.timeIndex(sensor_end_index)];
        end        
    end
end

inhomeTable = table; %('VariableNames',{'Subject','Source','StartTime','EndTime'});
inhomeTable.Subject = inhomeTable_Subject;
inhomeTable.Source = inhomeTable_Source;
inhomeTable.StartTime = inhomeTable_StartTime;
inhomeTable.EndTime = inhomeTable_EndTime;
clear inhomeTable_* v varnames daytable datenum w source_* valid_accel valid_beacon *_start_index *_end_index temp_index;

% write to file
inhomeFileName = fullfile(destFolder, ['dep' num2str(depID) '_inhome_all.csv']);
writetable(inhomeTable,inhomeFileName);      

%%    



















