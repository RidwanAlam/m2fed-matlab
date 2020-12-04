%%

ema_offset = minutes(15);


dataDir = 'C:\Users\Ridwan\Documents\MATLAB\M2FED Data';
clear datenum;
for d = 1:height(Deployments)

    depID = Deployments.DeploymentID(d);
    watchIDs = Deployments.WatchIDs{d};
    %modes = {'accel','battery','beacon'};
    beacons = Deployments.BeaconIDs{d};

    startDate = Deployments.StartDate(d); % deployment start-date 
    endDate = Deployments.EndDate(d);
%     startUnixTime = posixtime(datetime(startDate,'TimeZone','America/Los_Angeles'))*1000;
%     endUnixTime = posixtime(datetime(endDate,'TimeZone','America/Los_Angeles'))*1000;

    startDateNum = floor(datenum(startDate));

    startdayindex = 0;
    enddayindex = floor(datenum(endDate))-startDateNum;


    % dataDir = 'C:\Users\Ridwan\Documents\MATLAB\M2FED Data';
    % sourceFolder = fullfile(rawDataDir, depID ,"data_mat");
    destFolder = fullfile(dataDir,depID);
%
    % data available as "wWW_dDD_mode" in "prepDays_modes.mat"

    %load(fullfile(destFolder,'prepDays_accel.mat'));
    %run watchAccelMinutePrepScript.m;
    load(fullfile(destFolder,'prepDays_accel_minute.mat'));
    load(fullfile(destFolder,'prepDays_beacon.mat'));
    load(fullfile(destFolder,'prepDays_battery.mat'));
    
    
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
                modeDataName = [num2str(depID) '_w' num2str(watchID,'%02d')...
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

            eval(['daytable.accel_w' num2str(watchID,'%02d') '=' 'valid_accel;']); 
            eval(['daytable.beacon_w' num2str(watchID,'%02d') '=' 'valid_beacon;']); 
            eval(['daytable.battery_w' num2str(watchID,'%02d') '=' 'valid_battery;']); 
            eval(['daytable.charging_w' num2str(watchID,'%02d') '=' 'valid_charging;']); 

        end

        eval(['d' num2str(dayindex,'%02d') '_valid_minutes' '=' 'daytable;']);
        clear daytable new_timeIndex valid_* valid w new_timeIndex daytime_now kindex;
        eval(['clear D*_d' num2str(dayindex,'%02d') '_*']);

    end
    clear dayindex;


    % load EMA endtimes for each subject from excel
    % as DXX_ema_sSS
    % [char(depID) '_ema_s' tbl.Properties.VariableNames{k}]
    load(fullfile(destFolder,'prepDays_ema.mat'));
    ematables = whos([char(depID) '_ema_s*']);
    subjectIDs = [];
    % ema valid minutes
    

    varnames = whos('d*_valid_minutes');

    for v = 1:length(varnames)

        daytable = eval([varnames(v).name]);
        new_timeIndex = daytable.timeIndex;
        new_timeIndex.TimeZone='';
        %datenum = split(varnames(v).name,{'d','_valid_'});

        for s = 1:length(ematables) %length(watchIDs)

            ema_endtimes = eval(ematables(s).name);
            tempvar = split(ematables(s).name,'_ema_s');
            subjectIDs(s) = str2num(tempvar{2});
            clear tempvar;
            
            valid_ema = zeros(size(new_timeIndex));

            ema_endtimes_day = ema_endtimes((ema_endtimes>=new_timeIndex(1)-ema_offset) & ...
                (ema_endtimes<=new_timeIndex(end)+ema_offset));

            for d1 = 1:length(ema_endtimes_day)
                startIndex = find(new_timeIndex<=...
                    max(new_timeIndex(1),ema_endtimes_day(d1)-ema_offset),1,'last');
                endIndex = find(new_timeIndex>=...
                    min(new_timeIndex(end),ema_endtimes_day(d1)+ema_offset),1,'first');
                valid_ema(startIndex:endIndex) = 1;
            end
            eval(['daytable.ema_s' num2str(subjectIDs(s)) '=' 'valid_ema;']);   
            %eval(['valid_ema_d' datenum{2} '_w' num2str(w,'%02d') '=' 'valid_ema;']);

        end

        eval([varnames(v).name '=' 'daytable;']);
        clear s daytable new_timeIndex valid_ema ...
            startIndex endIndex ema_endtimes* d1
        
    end

    clear varnames v D*_ema_s* ematables 


    % save VALID 
    saveAsFile = fullfile(destFolder,'validDays_all.mat');
    save(saveAsFile, 'depID','watchIDs','subjectIDs','beacons','d*_valid_minutes','-v7.3');
    clear saveAsFile

    mkdir(fullfile(destFolder,'figs'));
    varnames = whos('d*_valid_minutes');
    for v = 1:length(varnames)
        daytable = eval([varnames(v).name]);%d03_valid_minutes;
        f = figure;
        for w = 1:length(watchIDs)
            subplot(length(watchIDs),1,w);
            temp_ema_data = eval(['daytable.ema_s' num2str(subjectIDs(w))]);
            temp_accel_data = eval(['daytable.accel_w' num2str(watchIDs(w),'%02d')]);
            temp_beacon_data = eval(['daytable.beacon_w' num2str(watchIDs(w),'%02d')]);
            plot(daytable.timeIndex,temp_ema_data,'g');hold on;
            plot(daytable.timeIndex,0.7*temp_accel_data,'m');hold on;
            plot(daytable.timeIndex,0.4*temp_beacon_data,'b');
            ylabel(['Subj-' num2str(w)]);
            ylim([-1.2 1.2]);
            legend('EMA','Acc','Loc');
            clear temp_*_data
        end

        filename = fullfile(destFolder, 'figs', ['daily_valid_' num2str(v,'%02d')]);
        print([char(filename) '.png'],'-dpng');savefig(f,[char(filename) '.fig']);
        close(f);

    end
    clear varnames daytable f filename;


    % times when in-home for sure:
    % subj S, watch W, source (Sensor/EMA), start time, end time

    inhomeTable_Subject = [];
    inhomeTable_Watch = [];
    inhomeTable_Source = [];
    inhomeTable_StartTime = [];
    inhomeTable_EndTime = [];

    varnames = whos('d*_valid_minutes');
    for v = 1:length(varnames)
        daytable = eval([varnames(v).name]);
        %datenums = split(varnames(v).name,{'d','_valid_'});

        for w = 1:length(watchIDs)
            % source: EMA 
            source_ema = eval(['daytable.ema_s' num2str(subjectIDs(w),'%02d')]);
            if any(source_ema)
                ema_start_index = min(find(source_ema(2:end)-source_ema(1:end-1)>0)+1,length(source_ema));
                %ema_end_index = min(find(source_ema(2:end)-source_ema(1:end-1)<0)+1,length(source_ema));
                inhomeTable_Subject = [inhomeTable_Subject; repelem(subjectIDs(w),length(ema_start_index))'];
                inhomeTable_Watch = [inhomeTable_Watch; repelem(watchIDs(w),length(ema_start_index))'];
                inhomeTable_Source = [inhomeTable_Source; repelem({'EMA'},length(ema_start_index))'];
                inhomeTable_StartTime = [inhomeTable_StartTime;daytable.timeIndex(ema_start_index)];
                inhomeTable_EndTime = [inhomeTable_EndTime;...
                    daytable.timeIndex(min(ema_start_index+(2*minutes(ema_offset)),length(source_ema)))];
            end
            % source: Sensor
            valid_accel = eval(['daytable.accel_w' num2str(watchIDs(w),'%02d')]);
            valid_beacon = eval(['daytable.beacon_w' num2str(watchIDs(w),'%02d')]);
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
                inhomeTable_Subject = [inhomeTable_Subject; repelem(subjectIDs(w),length(sensor_start_index))'];
                inhomeTable_Watch = [inhomeTable_Watch; repelem(watchIDs(w),length(sensor_start_index))'];
                inhomeTable_Source = [inhomeTable_Source; repelem({'Sensor'},length(sensor_start_index))'];
                inhomeTable_StartTime = [inhomeTable_StartTime;daytable.timeIndex(sensor_start_index)];
                inhomeTable_EndTime = [inhomeTable_EndTime;daytable.timeIndex(sensor_end_index)];
            end        
        end
    end


    inhomeTable = table; %('VariableNames',{'Subject','Source','StartTime','EndTime'});
    inhomeTable.SubjectID = inhomeTable_Subject;
    inhomeTable.WatchID = inhomeTable_Watch;
    inhomeTable.Source = inhomeTable_Source;
    inhomeTable.StartTime = inhomeTable_StartTime;
    inhomeTable.EndTime = inhomeTable_EndTime;
    clear inhomeTable_* v varnames daytable datenum w source_* valid_accel valid_beacon *_start_index *_end_index temp_index;
    clear d*_valid_minutes i

    % write to file
    inhomeFileName = fullfile(destFolder, [num2str(depID) '_inhome_all_20201203.csv']);
    writetable(inhomeTable,inhomeFileName);      
    clear inhomeTable inhomeFileName;


    clear depID watchID watchIDs subjectIDs beacons depStartDate depEndDate
    clear startUnixTime endUnixTime startDate endDate startDateNum startdayindex enddayindex
    clear sourceFolder destFolder


end
%%    
