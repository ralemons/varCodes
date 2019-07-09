%% Plot the data from the Enviromental Monitor in the lab

%% Read In the Data
clear;

% User Input
data = dataImport_EnvMon();


%% Plot Data in Menory

clear datePlot
numXticks = 30;
numAvg = 1;
% datePlot = datetime(2019,07,08,17,30,00); % Uncomment to plot a line at a date/time
dates = dateRange(data);


% Calculations
data{:,2} = movmean(data{:,2},numAvg);
data{:,3} = movmean(data{:,3},numAvg);
data{:,4} = movmean(data{:,4},numAvg);

Y(1,:) = [(min(data{dates,2})*0.99),(max(data{dates,2})*1.01)];
Y(2,:) = [(min(data{dates,3})*0.99),(max(data{dates,3})*1.01)];
Y(3,:) = [(min(data{dates,4})*0.999),(max(data{dates,4})*1.001)];

if exist('datePlot','var')
    closestData = find( ge(data{:,1},datePlot).*dates == 1 , 1 );
end

titles = {'Temperature Plot','Humidity Plot','Pressure Plot'};
units = {' C',' %RH',' ATM'};

%% Individual Plots
% Temp Plot
figure(1)
plot(data{dates,1},data{dates,2})
xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)
xlim([data{1,1},data{end,1}]);
title(titles{1});

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],Y(1,:),'Color',[1 0 0 .5],'LineStyle','--');
    text(mean(data{dates,1}),max(data{dates,2})+1,[datestr(datePlot), ': ', num2str(data{closestData,2}), units{1}],'FontSize',14,'Color',[1 0 0]);
end

% Humidity Plot
figure(2)
plot(data{dates,1},data{dates,3})
xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)
xlim([data{1,1},data{end,1}]);
title(titles{2});

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],Y(2,:),'Color',[1 0 0 .5],'LineStyle','--');
    text(mean(data{dates,1}),max(data{dates,3})+3,[datestr(datePlot), ': ', num2str(data{closestData,3}), units{2}],'FontSize',14,'Color',[1 0 0]);
end

% Pressure Plot
figure(3)
plot(data{dates,1},data{dates,4})
xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)
xlim([data{1,1},data{end,1}]);
title(titles{3});

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],Y(3,:),'Color',[1 0 0 .5],'LineStyle','--');
    text(mean(data{dates,1}),max(data{dates,4})*1.005,[datestr(datePlot), ': ', num2str(data{closestData,4}), units{3}],'FontSize',14,'Color',[1 0 0]);
end

%% Single Plot

yText = [max(data{dates,2}),max(data{dates,3}),max(data{dates,4})];

figure(99)
for ii = 1:3
subplot(3,1,ii)
plot(data{dates,1},data{dates,ii+1},'LineWidth',3)

title(titles{ii})
xlim([data{1,1},data{end,1}]);
ylim(Y(ii,:))



if ii < 3
    set(gca,'XTick',[]);
end

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],Y(ii,:),'Color',[1 0 0 1],'LineStyle','--');
    text(datePlot+hours(1),yText(ii),[datestr(datePlot), ': ', num2str(data{closestData,4}), units{ii}],'FontSize',14,'Color',[1 0 0]);
end

set(gca,'FontSize',12);

% fillFig(0,0)
fillFig(0.00,-0.05)
end

xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)


%% Helper Functions
function dates = dateRange(data)

t0 = data{1,1};
tF = data{end,1};

disp(['>> First possible date is: ',datestr(t0,'yyyy-mm-dd')]);
low = input('Enter First Date: yyyy-mm-dd ','s');
if ~strcmpi(low,'')
    low = datetime(low,'InputFormat','yyyy-MM-dd');
else
    low = data{1,1};
end

disp(['>> Last possible date is: ',datestr(tF,'yyyy-mm-dd')]);
high = input('Enter Last Date: yyyy-mm-dd ','s');
if ~strcmpi(high,'')
    high = datetime(high,'InputFormat','yyyy-MM-dd');
else
    high = data{end,1};
end

dates = ge(data{:,1},low) & le(data{:,1},high);


end