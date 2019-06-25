%% Main file for plotting data from RS Spec Analyzer.
% Comment and uncomment as needed.

%% Run this section to load the data
clear;
close all;

% Read in files
[dataIN,fileStr] = dataImport_SpecAnalyser_Rigol(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);

% Turn off warning about number of peaks. QOL change, can be removed.
warning('off','signal:findpeaks:largeMinPeakHeight');



%% Run this section to process the data and find peaks

dataOUT = dataIN;

% % Offset for any attenuation
% dataOUT(:,2,:) = dataOUT(:,2,:)+10;

% %%% Select ranges for data analysis
% highPass = 31.1e6;
% [~,shift] = min(abs(dataOUT(:,1,1)-highPass));
% if shift ~= 1
%     dataOUT(1:shift,:,:) = [];
% end
%
% lowPass = 31.8e6;
% [~,shift] = min(abs(dataOUT(:,1,1)-lowPass));
% if shift ~= length(dataOUT(:,1,1))
%     dataOUT(shift:end,:,:) = [];
% end



for ii = 1:numFilesIn
    
    
    % Average value
    meanVal(ii) = mean(dataOUT(:,2,ii));
    
    
    %     % Flatten base based on average value
    %     biasVal = 5;
    %     vals = dataOUT(:,2,ii) < (meanVal(ii) + biasVal);
    %     dataOUT(vals,2,ii) = meanVal(ii);
    
    
    % Absolute Min/Max Value
    minVal(ii) = min(dataOUT(:,2,ii));
    [maxVal(ii),pos(ii)] = max(dataOUT(:,2,ii));
    val(ii) = dataOUT(pos(ii),1,ii);
    valMHZ(ii) = val(ii)/1e6;
    
    % Find all peaks
    [pks{ii},pkPos{ii}] = findpeaks(dataOUT(:,2,ii),dataOUT(:,1,ii),...
        'MinPeakHeight',meanVal(ii)+3,'MinPeakProminence',15);
    
    
end

% %% Create averages of data sets.
%
% % Sum of the data, Shifted to zero
% dataOUT(:,1,numFilesIn+1) = dataOUT(:,1,1);
% dataOUT(:,2,numFilesIn+1) = sum(dataOUT(:,2,:),3)-sum(meanVal);
% minVal(numFilesIn+1) = min(dataOUT(:,2,numFilesIn+1));
% maxVal(numFilesIn+1) = max(dataOUT(:,2,numFilesIn+1));
%
% pks{numFilesIn+1,1} = [];
% pkPos{numFilesIn+1,1} = [];
%
% %% Analyze the data
%
% diffFreqs = pkPos{1}'/1e6;
% diffFreqs = diffFreqs(2:end)-diffFreqs(1:end-1)




%% Run this section to plot the data

numPlot = 1:size(dataOUT,3);
% numPlot = numFilesIn;

titles = {'AOFS Signal, MAXHOLD',...
    'AOFS Signal, AVERAGE',...
    'AOFS Signal, FREE RUN'};
xLabels = repmat({'Frequency (MHz)'},numFilesIn+1,1);
yLabels = repmat({'Power (dBc)'},numFilesIn+1,1);


for ii = numPlot
    figure(ii);
    ax = gca;
    
    
    plot(dataOUT(:,1,ii)/1e6,dataOUT(:,2,ii)-maxVal(ii),...
        'linewidth',4,'color',[33 54 86]/255);
    xlim([min(dataOUT(:,1,ii)/1e6) max(dataOUT(:,1,ii)/1e6)]);
    ylim([minVal(ii)-maxVal(ii) maxVal(ii)+10-maxVal(ii)]);
    
    
    % % Absolute Max Line
    % line([val val],[minVal(ii),maxVal(ii)+10],'Color',[1 0 0 .5],'LineStyle','--');
    % format shortG
    % text(82e6,-85,['Beat Signal:  ', num2str(valMHZ),' Mhz, ', num2str(maxVal), ' dBm'],'FontSize',14,'Color',[1 0 0]);
    
    % % Specific Line
    % line([valSpec valSpec],[minVal(ii),maxVal(ii)+10],'Color',[.1333 .545 .1333 .5],'LineStyle','--');
    % format shortG
    % text(82e6,-90,['Beat Signal:  ', num2str(valMHZSpec),' Mhz, ', num2str(maxValSpec), ' dBm'],'FontSize',14,'Color',[.1333 .545 .1333]);
    
%     textUse = {'$f_{LO}$','$f_{REP} - f_{LO}$','$f_{REP}$'};
%     textUse = {'$f_{CEO}$','$f_{REP} - f_{CEO}$','$f_{REP}$'};
    % All lines for peaks
    for jj = 1:length(pkPos{ii})
        pos = pkPos{ii}(jj)/1e6;
        format shortG
        line([pos pos],[minVal(ii)-maxVal(ii),pks{ii}(jj)-maxVal(ii)],...
            'Color',[1 0 0 .5],...
            'LineStyle','--',...
            'LineWidth',5);
        t = text(pos,pks{ii}(jj)+3-maxVal(ii),...
            [num2str(pkPos{ii}(jj)/1e6,'%.4g'),' Mhz'],...
            'FontSize',38,...
            'Color',[0 0 0],...
            'HorizontalAlignment','center',...
            'Interpreter','latex');
        set(t,'Rotation',00);
    end
    

    
    % Make it look nice
    xlabel(xLabels{ii},'FontSize',24);
    ylabel(yLabels{ii},'FontSize',24);
    grid on
    ax.FontSize = 40;
%     title(titles{ii},'FontSize',40);

    
end
%% Special Plots %%

numPlot = 1:size(dataOUT,3);
% numPlot = numFilesIn;
% numPlot = 2;


%%%%% Bandwidth Calculation %%%%%


titles = {'AOFS Signal, MAXHOLD',...
    'AOFS Signal, AVERAGE',...
    'AOFS Signal, FREE RUN'};

xLabels = repmat({'Frequency (MHz)'},numFilesIn+1,1);
yLabels = repmat({'Power (dBm)'},numFilesIn+1,1);

for ii = numPlot
    
    figure(99-ii)
    ax = gca;
    
    scale = 1e6;
    
    % Get BW measurements
    [bw,fLo,fHi,pow] = obw(db2pow(dataOUT(:,2,ii)),dataOUT(:,1,ii));
    
    plot(dataOUT(:,1,ii)/scale,dataOUT(:,2,ii),...
        'linewidth',3);
    xlim([min(dataOUT(:,1,ii)/scale) max(dataOUT(:,1,ii)/scale)]);
    ylim([minVal(ii) maxVal(ii)+10]);
    
    
    
    % fLo and fHi lines
    color = ax.ColorOrder(6,:);
    l = line([fLo/scale fLo/scale],[minVal(ii),maxVal(ii)+10],...
        'Color',color,'LineWidth',2);
    l.Color(4) = 0.6;
    l = line([fHi/scale fHi/scale],[minVal(ii),maxVal(ii)+10],...
        'Color',color,'LineWidth',2);
    l.Color(4) = 0.6;
    p = patch([fLo/scale fHi/scale fHi/scale fLo/scale],...
        [minVal(ii) minVal(ii) maxVal(ii)+10 maxVal(ii)+10],...
        color,'FaceAlpha',0.4,'EdgeColor','k','LineStyle','none');
    
    % All lines for peaks
    t = text((fHi+bw/5)/scale,max(pks{ii})+1,...
        ['Bandwidth: ',num2str(bw/1e3,'%.4g'),' kHz'],...
        'FontSize',18,...
        'Color',[0 0 0],...
        'HorizontalAlignment','left');
    
    
    ax.FontSize = 24;
    
    
    % Make it look nice
    title(titles{ii},'FontSize',30);
    xlabel(xLabels{ii},'FontSize',24);
    ylabel(yLabels{ii},'FontSize',24);
    
    
end







