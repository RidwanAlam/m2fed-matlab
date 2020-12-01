%%

dataHome = 'C:\Users\Ridwan\Documents\MATLAB\M2FED Data';

for d = 3:length(depIDs)

depID = depIDs{d};
watchIDs = watchIDss{d};
%modes = {'accel','battery','beacon'};
beacons = beaconIDs{d};

depStartDate = [depStartDates{d} ' 00:00:00']; % deployment start-date 
depEndDate = [depEndDates{d} ' 23:59:59'];


% depID = 30;
% d = find(cell2mat(depIDs)==depID);
% watchIDs = watchIDss{d};
% beacons = beaconIDs{d};
% depStartDate = [depStartDates{d} ' 00:00:00']; % deployment start-date 
% depEndDate = [depEndDates{d} ' 23:59:59'];


startUnixTime = posixtime(datetime(depStartDate,'TimeZone','America/Los_Angeles'))*1000;
endUnixTime = posixtime(datetime(depEndDate,'TimeZone','America/Los_Angeles'))*1000;

startDate = datetime((startUnixTime)./1000,'ConvertFrom','posixtime',...
    'TimeZone','America/Los_Angeles');
endDate = datetime((endUnixTime)./1000,'ConvertFrom','posixtime',...
    'TimeZone','America/Los_Angeles');

startDateNum = floor(datenum(startDate));

startdayindex = 0;
enddayindex = floor(datenum(endDate))-startDateNum;

sourceFolder = fullfile(dataHome,['Processed Data\deployment_',num2str(depID),'\data_mat']);
destFolder = fullfile(dataHome,['dep_',num2str(depID)]);

% saveAsFile = fullfile(destFolder,'info.mat');
% save(saveAsFile, 'depID','watchIDs','beacons','depStartDate','depEndDate','-v7.3');

%
% data available as "wWW_dDD_mode" in "prepDays_modes.mat"

load(fullfile(destFolder,'prepDays_accel.mat'));
run watchAccelMinutePrepScript.m;
%load(fullfile(destFolder,'prepDays_accel_minute.mat'));
load(fullfile(destFolder,'prepDays_beacon.mat'));
load(fullfile(destFolder,'prepDays_battery.mat'));
load(fullfile(destFolder,'prepDays_ema.mat'));

%


%

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
            modeDataName = ['dep' num2str(depID) '_w' num2str(watchID,'%02d')...
                '_d' num2str(dayindex,'%02d') '_' mode];
            if exist(modeDataName,'var')==1
                
                modeData = eval(modeDataName);
                times = modeData.timeIndex;
                kindex = (times.Hour*60)+(times.Minute)+1;
                if mode=="accel_minute"
                    values = modeData.teager;
                    valid_accel(kindex(values>=0 & values<0.01)) = 0;
                    valid_accel(kindex(values>=0.01)) = 1;
                    
                elseif mode=="beacon"
                    [rrssi,curroom] = max(table2array(modeData(:,2:end)),[],2);
                    curroom(rrssi==-120) = -1;                    
                    rrssi = max(rrssi+120,0.1);
                    dist = 2*abs(log(100)-log(rrssi));
                    prob = min(1,2*exp(-dist));
                    
                    %values = double(prob>0.1).*curroom;
                    %valid_beacon(kindex(curroom>0)) = 0;
                    valid_beacon(kindex(curroom>0)) = 1;
                    clear rrssi curroom dist prob;
                    
                elseif mode=="battery"
                    valid_charging(kindex(modeData.charging==1)) = 1;
                    valid_charging(kindex(modeData.charging==0)) = 0;
                    
                    valid_battery(kindex(modeData.charging>=0)) = ...
                        modeData.power(modeData.charging>=0);
                    
                end
                
            end
        end
        clear m modeData mode modeDataName values times;
        % charging: battery = accel = beacon = 0 / dead?
        valid_accel(valid_charging==1) = 0;
        valid_beacon(valid_charging==1) = 0;
                
        eval(['daytable.accel_w' num2str(w,'%02d') '=' 'valid_accel;']); 
        eval(['daytable.beacon_w' num2str(w,'%02d') '=' 'valid_beacon;']); 
        eval(['daytable.battery_w' num2str(w,'%02d') '=' 'valid_battery;']); 
        eval(['daytable.charging_w' num2str(w,'%02d') '=' 'valid_charging;']); 
            
    end
    
    eval(['d' num2str(dayindex,'%02d') '_valid_minutes' '=' 'daytable;']);
    clear daytable new_timeIndex valid_* valid w new_timeIndex daytime_now kindex;
    eval(['clear dep*_d' num2str(dayindex,'%02d') '_*']);

end
clear dayindex;


% load EMA endtimes for each subject from excel
% as depXX_ema_sSS

% ema valid minutes

varnames = whos('d*_valid_minutes');

for v = 1:length(varnames)
    
    daytable = eval([varnames(v).name]);
    new_timeIndex = daytable.timeIndex;
    new_timeIndex.TimeZone='';
    %datenum = split(varnames(v).name,{'d','_valid_'});
    
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
        %eval(['valid_ema_d' datenum{2} '_w' num2str(w,'%02d') '=' 'valid_ema;']);
        
    end
    
    eval([varnames(v).name '=' 'daytable;']);
    clear w daytable new_timeIndex valid_ema startIndex endIndex offset ema_endtimes* d
end

clear varnames v dep*_ema_s* 


% save VALID 
saveAsFile = fullfile(destFolder,'validDays_all.mat');
save(saveAsFile, 'depID','watchIDs','beacons','d*_valid_minutes','-v7.3');
clear saveAsFile

mkdir(fullfile(destFolder,'figs'));
varnames = whos('d*_valid_minutes');
for v = 1:length(varnames)
    daytable = eval([varnames(v).name]);%d03_valid_minutes;
    f = figure;
    for w = 1:length(watchIDs)
        subplot(length(watchIDs),1,w);
        temp_ema_data = eval(['daytable.ema_w' num2str(w,'%02d')]);
        temp_accel_data = eval(['daytable.accel_w' num2str(w,'%02d')]);
        temp_beacon_data = eval(['daytable.beacon_w' num2str(w,'%02d')]);
        plot(daytable.timeIndex,temp_ema_data,'g');hold on;
        plot(daytable.timeIndex,0.7*temp_accel_data,'m');hold on;
        plot(daytable.timeIndex,0.4*temp_beacon_data,'b');
        ylabel(['Subj-' num2str(w)]);
        legend('Acc','Loc','EMA');
        clear temp_*_data
    end
    
    %     subplot(4,1,2);
    %     plot(daytable.timeIndex,daytable.ema_w02,'g');hold on;
    %     plot(daytable.timeIndex,0.7*daytable.accel_w02,'m');hold on;
    %     plot(daytable.timeIndex,0.4*daytable.beacon_w02,'b');    
    %     ylabel(['Subj-2']);
    %     subplot(4,1,3);
    %     plot(daytable.timeIndex,daytable.ema_w03,'g');hold on;
    %     plot(daytable.timeIndex,0.7*daytable.accel_w03,'m');hold on;
    %     plot(daytable.timeIndex,0.4*daytable.beacon_w03,'b');
    %     ylabel(['Subj-3']);
    %     subplot(4,1,4);
    %     plot(daytable.timeIndex,daytable.ema_w04,'g');hold on;
    %     plot(daytable.timeIndex,0.7*daytable.accel_w04,'m');hold on;
    %     plot(daytable.timeIndex,0.4*daytable.beacon_w04,'b');    
    %     ylabel(['Subj-4']);
    
    filename = fullfile(destFolder, 'figs', ['daily_valid_' num2str(v,'%02d')]);
    print([filename '.png'],'-dpng');savefig(f,[filename '.fig']);
    close(f);
    
end
clear varnames daytable f filename;





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
clear d*_valid_minutes i

% write to file
inhomeFileName = fullfile(destFolder, ['dep' num2str(depID) '_inhome_all.csv']);
writetable(inhomeTable,inhomeFileName);      
clear inhomeTable inhomeFileName;


clear depID watchID watchIDs beacons depStartDate depEndDate
clear startUnixTime endUnixTime startDate endDate startDateNum startdayindex enddayindex
clear sourceFolder destFolder

end
%%    
