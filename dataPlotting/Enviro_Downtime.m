%% Find number of points in each day, Ie. when logger is down

data = readtable('/home/randy/Documents/research/SLAC/Summer18/CEP/LoggerData_20180803-20190116.csv');
data = data{:,1};

tlower = datetime(2018,08,03,00,00,00);
tupper = datetime(2018,08,03,23,59,59);

dur = round(days( data(end)-data(1) ));
numPoints = zeros(dur,1);


for i = 1:dur
    
    numPoints(i) = sum( isbetween(data,tlower,tupper) );
    tlower = tlower + 1;
    tupper = tupper + 1;
    
end

bar(numPoints,1,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)