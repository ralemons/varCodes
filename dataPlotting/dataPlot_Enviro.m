%% Plot the data from the Enviromental Monitor in the lab

%% Read In the Data
clear;

% User Input
data = dataImport_EnvMon();


%% Plot Data in Menory

numXticks = 30;
numAvg = 6;
datePlot = datetime(2018,10,06,16,00,00); % Uncomment to plot a line at a date/time
dates = dateRange(data);


% Calculations
data{:,2} = movmean(data{:,2},numAvg);
data{:,3} = movmean(data{:,3},numAvg);
data{:,4} = movmean(data{:,4},numAvg);

tempY = [(min(data{dates,2})*0.9),(max(data{dates,2})*1.1)];
humidY = [(min(data{dates,3})*0.9),(max(data{dates,3})*1.1)];
pressY = [(min(data{dates,4})*0.99),(max(data{dates,4})*1.01)];

if exist('datePlot','var')
    closestData = find( ge(data{:,1},datePlot).*dates == 1 , 1 );
end

% Temp Plot
figure(1)
plot(data{dates,1},data{dates,2})
xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)
% xlim([data{1,1}-2,data{end,1}+2]);
ylim(tempY);
title('Temperature Plot');

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],tempY,'Color',[1 0 0 .5],'LineStyle','--');
    text(mean(data{dates,1}),max(data{dates,2})+1,[datestr(datePlot), ': ', num2str(data{closestData,2}), ' C'],'FontSize',14,'Color',[1 0 0]);
end

% Humidity Plot
figure(2)
plot(data{dates,1},data{dates,3})
xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)
% xlim([data{1,1}-2,data{end,1}+2]);
ylim(humidY);
title('Humidity Plot');

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],humidY,'Color',[1 0 0 .5],'LineStyle','--');
    text(mean(data{dates,1}),max(data{dates,3})+3,[datestr(datePlot), ': ', num2str(data{closestData,3}), ' %RH'],'FontSize',14,'Color',[1 0 0]);
end

% Pressure Plot
figure(3)
plot(data{dates,1},data{dates,4})
xticks(data{1:floor(length(data{dates,1})/numXticks):end,1})
xtickangle(60)
% xlim([data{1,1}-2,data{end,1}+2]);
ylim(pressY);
title('Pressure Plot');

% Date Line
if exist('datePlot','var')
    line([datePlot,datePlot],pressY,'Color',[1 0 0 .5],'LineStyle','--');
    text(mean(data{dates,1}),max(data{dates,4})*1.005,[datestr(datePlot), ': ', num2str(data{closestData,4}), ' ATM'],'FontSize',14,'Color',[1 0 0]);
end


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