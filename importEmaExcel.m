%% Import data from spreadsheet

emaHomeFolder = fullfile(dataDir,'EMA');

for d = 1:height(Deployments)

    depID = Deployments.DeploymentID(d);
    %depID = D310;
    emaFile = fullfile(emaHomeFolder,[char(depID) '_ema.xlsx']);

    % Setup the Import Options
    opts = detectImportOptions(emaFile, 'Sheet','endtimes',...
        'ReadVariableNames',true,'PreserveVariableNames',true,...
        'VariableNamesRange','A1','DataRange','A2');

    % Import the data
    tbl = readtable(emaFile, opts);

    % Convert to output type
    for k = 1:width(tbl)
        tempVar = table2array(tbl(:,k));
        tempVar(isnat(tempVar))=[];
        eval([char(depID) '_ema_s' tbl.Properties.VariableNames{k} '=' 'tempVar;']);
        clear tempVar;
    end

    clear k tbl opts emaFile;

    % save as mat file
    destFolder = fullfile(dataDir,depID);
    saveAsFile = fullfile(destFolder,'prepDays_ema.mat');
    save(saveAsFile, 'depID','*_ema_s*','-v7.3');

    clear *_ema_s* depID saveAsFile;
end