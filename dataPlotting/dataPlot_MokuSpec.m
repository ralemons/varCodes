% General file to plot a number of different possible aspects. Comment and
% uncomment as needed.

%% Run this section to load the data
clear;
close all;

% Read in files
[dataIN,fileStr] = dataImport_Moku(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);

% Turn off warning about number of peaks. QOL change, can be removed.
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Run this section to process the data and find peaks

dataOUT = dataIN;
numChans = size(dataOUT,2)-1;

%
% %%% Select ranges for data analysis
% highPass = 70e6;
% [~,shift] = min(abs(dataOUT(:,1,1)-highPass));
% if shift ~= 1
%     dataOUT(1:shift,:,:) = [];
% end
%
% lowPass = 90e6;
% [~,shift] = min(abs(dataOUT(:,1,1)-lowPass));
% if shift ~= length(dataOUT(:,1,1))
%     dataOUT(shift:end,:,:) = [];
% end


%%% Perform operations on each file
for ii = 1:numFilesIn
    
    % % Offset for any attenuation
    % dataOUT(:,2,ii) = dataOUT(:,2,ii)+10;
    % dataOUT(:,3,ii) = dataOUT(:,3,ii)+10;
    
    for jj = 1:numChans
        
        % Average value
        meanVal(ii,jj) = mean(dataOUT(:,1+jj,ii));
        
        %     % Flatten base based on average value
        %     biasVal = 5;
        %     vals = dataOUT(:,1+jj,ii) < (meanVal(ii,jj) + biasVal);
        %     dataOUT(vals,1+jj,ii) = meanVal(ii,jj);
        
        
        % Absolute Min/Max Value
        minVal(ii,jj) = min(dataOUT(:,1+jj,ii));
        [maxVal(ii,jj),pos(ii,jj)] = max(dataOUT(:,1+jj,ii));
        val(ii,jj) = dataOUT(pos(ii,jj),1+jj,ii);
        valMHZ(ii,jj) = val(ii,jj)/1e6;
        
        % Find all peaks channel 1
        [pks{ii,jj},pkPos{ii,jj}] = findpeaks(dataOUT(:,1+jj,ii),dataOUT(:,1,ii),...
            'MinPeakHeight',meanVal(ii,jj)+3,'MinPeakProminence',30);
        
    end
end

% %% Analyze frequency comb peaks
%
% diffFreqs = pkPos{1}'/1e6;
% diffFreqs = diffFreqs(2:end)-diffFreqs(1:end-1)


%% Run this section to plot the data


chansPlot = 1:numChans;
filesPlot = 1:numFilesIn;

titles = {'AOFS Signal',...
    'CEP Measurement in OOL F2F'};
xLabels = repmat({'Frequency (MHz)'},numFilesIn+1,1);
yLabels = repmat({'Power (dBm)'},numFilesIn+1,1);

colorsPlot = {[33 54 86]/255,[124 33 0]/255};


for ii = filesPlot
    
    figure(ii);
    ax = gca;
    
    for jj = chansPlot
        
        hold on;
            p(jj) = plot(dataOUT(:,1,ii)/1e6,dataOUT(:,1+jj,ii),'linewidth',2,'color',colorsPlot{jj});
        hold off;
        %     xticks(highPass/1e6:.05:lowPass/1e6);
        xlim([min(dataOUT(:,1,ii)/1e6) max(dataOUT(:,1,ii)/1e6)]);
        ylim([min(minVal(ii,:)) max(maxVal(ii,:))+25]);
        
        
%         % Absolute Max Line
%         line([val val],[minVal(ii),maxVal(ii)+10],'Color',[1 0 0 .5],'LineStyle','--');
%         format shortG
%         text(82e6,-85,['Beat Signal:  ', num2str(valMHZ),' Mhz, ', num2str(maxVal), ' dBm'],'FontSize',14,'Color',[1 0 0]);
%         
%         % Specific Line
%         line([valSpec valSpec],[minVal(ii),maxVal(ii)+10],'Color',[.1333 .545 .1333 .5],'LineStyle','--');
%         format shortG
%         text(82e6,-90,['Beat Signal:  ', num2str(valMHZSpec),' Mhz, ', num2str(maxValSpec), ' dBm'],'FontSize',14,'Color',[.1333 .545 .1333]);
        
        
        % All lines for each channel
        for kk = 1:length(pkPos{ii,jj})
            pos = pkPos{ii,jj}(kk)/1e6;
            format shortG
            line([pos pos],[min(minVal(ii,:)),pks{ii,jj}(kk)],...
                'Color',[1 0 0 .8],...
                'LineStyle','--');
            t = text(pos,pks{ii,jj}(kk)+1,...
                [num2str(pkPos{ii,jj}(kk)/1e6,'%.4g'),' MHz'],...
                'FontSize',28,...
                'Color',[0 0 0],...
                'HorizontalAlignment','left');
            set(t,'Rotation',45);
        end
        
    end
    
    
    % Make it look nice
%     title(titles{ii},'FontSize',40);
    title('Moku Spectrum Analyzer, Crosstalk','FontSize',40);
    legend([p(1) p(2)],{'Channel 1','Channel 2'},'Location','northeast');
    xlabel(xLabels{ii},'FontSize',34);
    ylabel(yLabels{ii},'FontSize',34);
    ax.FontSize = 40;
    
end