%%
% initialize day-wise indices for the specific deployment/study


%% Data Loader - Deployment Summary
% Output: Deployments

dataDir = "C:\Users\Ridwan\Documents\MATLAB\M2FED Data";
rawDataDir = fullfile(dataDir,"Raw Data");
depFileName = fullfile(dataDir, "Deployment_summary_9.29.19.xlsx");
Deployments = importDeployments(depFileName,"Summary", [8, 27]);

depStartDates = Deployments.StartDate; % deployment start-date 
depEndDates = Deployments.EndDate + hours(24) - seconds(1);

startUnixTimes = posixtime(datetime(depStartDates,'TimeZone','America/Los_Angeles'))*1000;
endUnixTimes = posixtime(datetime(depEndDates,'TimeZone','America/Los_Angeles'))*1000;

startDates = datetime((startUnixTimes)./1000,'ConvertFrom','posixtime',...
    'TimeZone','America/Los_Angeles');
endDates = datetime((endUnixTimes)./1000,'ConvertFrom','posixtime',...
    'TimeZone','America/Los_Angeles');

Deployments.StartDate = startDates;
Deployments.EndDate = endDates;

clear startDates endDates startUnixTimes endUnixTimes;

watchIDs = {};
beaconIDs = {};

for i = 1:height(Deployments)
    watchIDs{i,1} = str2num(Deployments.WatchIDs(i));
    beaconIDs{i,1} = str2num(Deployments.BeaconIDs(i));
end
Deployments.WatchIDs = watchIDs;
Deployments.BeaconIDs = beaconIDs;    
% Deployments.WatchIDs{4} = [3,8,9];
% Deployments.BeaconIDs{4} = [3,4,5,6,8,10,12];
clear watchIDs beaconIDs i dep*;

%%
%%

modes = {'accel','battery','beacon'};

for d = 1:height(Deployments)

    depID = Deployments.DeploymentID(d);
    watchIDs = Deployments.WatchIDs{d};
    %modes = {'accel','battery','beacon'};
    beacons = Deployments.BeaconIDs{d};

    startDate = Deployments.StartDate(d); % deployment start-date 
    endDate = Deployments.EndDate(d);
    startUnixTime = posixtime(datetime(startDate,'TimeZone','America/Los_Angeles'))*1000;
    endUnixTime = posixtime(datetime(endDate,'TimeZone','America/Los_Angeles'))*1000;

    startDateNum = floor(datenum(startDate));

    startdayindex = 0;
    enddayindex = floor(datenum(endDate))-startDateNum;


    % dataDir = 'C:\Users\Ridwan\Documents\MATLAB\M2FED Data';
    sourceFolder = fullfile(rawDataDir, depID ,"data_mat");

    % data available as "modes_wXX" in "modes.mat"
    load(fullfile(sourceFolder,'accel.mat'));
    load(fullfile(sourceFolder,'beacon.mat'));
    load(fullfile(sourceFolder,'battery.mat'));

    destFolder = fullfile(dataDir,depID);
    mkdir(destFolder);

    saveAsFile = fullfile(destFolder,'info.mat');
    save(saveAsFile, 'depID','watchIDs','beacons','startDate','endDate','modes','-v7.3');
    clear saveAsFile;

    % sort data in day-wise tables
    
    for m = 1:length(modes)    
        mode = modes{m};
        for i = 1:length(watchIDs)
            watchID = watchIDs(i);
            modeData = table; 
            if mode=="accel"
                if exist([mode '_w' num2str(watchID)],'var')
                    modeData.timestamp = eval([mode '_w' num2str(watchID) '(:,1)']);
                    modeData.x = eval([mode '_w' num2str(watchID) '(:,2)']);
                    modeData.y = eval([mode '_w' num2str(watchID) '(:,3)']);
                    modeData.z = eval([mode '_w' num2str(watchID) '(:,4)']);
                end
            elseif mode=="battery"
                if exist([mode '_w' num2str(watchID)],'var')
                    modeData.timestamp = eval([mode '_w' num2str(watchID) '(:,2)']);
                    modeData.power = eval([mode '_w' num2str(watchID) '(:,5)']);
                    modeData.charging = eval([mode '_w' num2str(watchID) '(:,4)']); % 501 = charging, 502 = not;
                end
            elseif mode=="beacon"
                if exist([mode '_w' num2str(watchID)],'var')
                    modeData.timestamp = eval([mode '_w' num2str(watchID) '(:,3)']);
                    modeData.beaconID = eval([mode '_w' num2str(watchID) '(:,1)']);
                    modeData.rssi = eval([mode '_w' num2str(watchID) '(:,5)']);
                    modeData.rssi = double(modeData.rssi);
                    modeData.beaconID = double(modeData.beaconID);
                end
            end

            if ~isempty(modeData)
                times=modeData.timestamp;
                indx = find(times<endUnixTime & times>startUnixTime);
                modeData = modeData(indx,:);
                modeData.timeIndex = datetime((modeData.timestamp)./1000,'ConvertFrom',...
                    'posixtime','TimeZone','America/Los_Angeles');
                modeData.timestamp = [];
            end

            % sort data to day-wise tables
            while(1)
                if height(modeData)==0
                    break;
                end
                tempDateNum = floor(datenum(modeData.timeIndex(1)));
                tempInd = find(floor(datenum(modeData.timeIndex))==tempDateNum);


                newTableName = ['w' num2str(watchID,'%02d') '_d' num2str(tempDateNum-startDateNum,'%02d') '_' mode];
                if (~isempty(tempInd))
                    eval([newTableName '=' 'modeData(tempInd,:);']);
                    eval([newTableName '=' 'sortrows(' newTableName ',{''timeIndex''});']);
                    %eval(['save(''' saveAsFileName ''',''' newTableName ''',''-append'');']);
                    modeData(tempInd,:) = [];
                else
                    eval([newTableName '=' 'table;']);
                end
                clear newTableName;
            end

            clear tempDateNum tempInd indx

        end

    end


    run watchBeaconPrepScript.m;
    run watchBatteryPrepScript.m;
    run watchAccelPrepScript.m;
    run watchAccelMinutePrepScript.m;


    clear i j k r m mode modeData times
    clear accel* *accel battery* *battery beacon_* *beacon 
    clear start* end* sourceFolder destFolder 
    clear depStartDate depEndDate depID beacons watchIDs watchID

end


run importEmaExcel.m;


%%
% 
% % Save figures for each day and each Pebble 
% % Helpful for manual inspection
% 
% colors = {'r','g','b','k'};
% 
% for dayindex = startdayindex:1:enddayindex
%     
%     f = figure(1);
%     for w = 1:length(watchIDs)
%         watchID = watchIDs(w);
%         for m = 1:length(modes)
%             mode = modes{m};
%             modeDataName = ['w' num2str(watchID,'%02d') '_d' num2str(dayindex,'%02d') '_' mode];
%             if exist(modeDataName,'var')==1
%                 modeData = eval(modeDataName);
%                 times = modeData.timeIndex;
%                 if mode=="accel"
%                     values = sqrt(((modeData.x).*(modeData.x))+...
%                         ((modeData.y).*(modeData.y))+((modeData.z).*(modeData.z)));
%                     if exist('ax1','var')==0
%                         ax1 = subplot(3,1,1);
%                     end
%                     plot(ax1,times,values,'Color',colors{w});hold on;
%                 elseif mode=="beacon"
%                     dist = log(100)-log(modeData.rssi + 110);
%                     prob = min(1,1.8*exp(-dist));
%                     beaconID = modeData.beaconID;
%                     values = (prob>0.4).*beaconID;                   
%                     if exist('ax2','var')==0
%                         ax2 = subplot(3,1,2);
%                     end
%                     plot(ax2,times,values,'Color',colors{w});hold on;
% 
%                 elseif mode=="battery"
%                     values = modeData.power;
%                     if exist('ax3','var')==0
%                         ax3 = subplot(3,1,3);
%                     end
%                     plot(ax3,times,values,'Color',colors{w});hold on;
%                 end
%                 
%             end
%         end
%     end
%     
%     filename = ['..\M2FED Data\dep_' num2str(depID,'%02d') '\figs\f' num2str(dayindex,'%02d') ];
%     print([filename '.png'],'-dpng');savefig(f,[filename '.fig']);
%     close(f);
%     clear ax1 ax2 ax3 f;
% end



%%

