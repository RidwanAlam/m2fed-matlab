%% Import data from spreadsheet

emaHomeFolder = 'C:\Users\Ridwan\Project M2FED\EMA-Data-from-Meiyi-20190810';

depID = 310;
emaFile = fullfile(emaHomeFolder,['dep_' num2str(depID) '_ema.xlsx']);

% Setup the Import Options
opts = detectImportOptions(emaFile, 'Sheet','endtimes', 'Range','A2');

% Import the data
tbl = readtable(emaFile, opts, "UseExcel", false);

% Convert to output type
for k = 1:width(tbl)
    tempVar = table2array(tbl(:,k));
    tempVar(isnat(tempVar))=[];
    eval(['dep' num2str(depID) '_ema_s' num2str(k,'%02d') '=' 'tempVar;']);
    clear tempVar;
end

clear k tbl opts emaFile;

% save as mat file
destFolder = fullfile(dataHome,['dep_',num2str(depID)]);
saveAsFile = fullfile(destFolder,'prepDays_ema.mat');
save(saveAsFile, 'depID','dep*_ema_s*','-v7.3');

clear dep*_ema_s* depID saveAsFile;