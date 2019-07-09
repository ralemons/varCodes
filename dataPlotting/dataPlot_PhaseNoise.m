% General file to plot a number of different possible aspects. Comment and
% uncomment as needed.

%% Run this section to load the data
clear;
close all;

% Read in files
[dataIN,headers,fileStr] = dataImport_PhaseNoise(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);

% Turn off warning about number of peaks. QOL change, can be removed.
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Run this section to process the data and find peaks

dataOUT = dataIN;

optFreq = 1.934e14; % Optical Frequency in Hz

% fileOrder = {'Only CEP','ULM w/ LOCSET','ULM w/o LOCSET'};
fileOrder = {'Clock','Locked'};




%%%% Frequency Range %%%%%

%%% Select ranges for data analysis
highPass = .01;
[~,shift] = min(abs(dataOUT(:,1,1)-highPass));
if shift ~= 1
    dataOUT(1:shift,:,:) = [];
end

lowPass = 3e6;
[~,shift] = min(abs(dataOUT(:,1,1)-lowPass));
if shift ~= length(dataOUT(:,1,1))
    dataOUT(shift:end,:,:) = [];
end



for ii = 1:numFilesIn
    
    loopBW(ii) = str2double(headers{88,2,ii});
    carrier.freq(ii) = str2double(headers{72,2,ii});
    carrier.volt(ii) = str2double(headers{74,2,ii});
    phaseMod(ii) = str2double(headers{78,2,ii});
    rmsJitter(ii) = str2double(headers{80,2,ii});
    intPhaseNoise.dBc(ii) = str2double(headers{81,2,ii});
    intPhaseNoise.rad(ii) =...
        sqrt( 2 * trapz(dataOUT(:,1,ii) , (10 .^ (dataOUT(:,2,ii) ./ 10) ) ) );
    jitter.rf(ii) = (1/(2*pi*carrier.freq(1))) * intPhaseNoise.rad(ii);
    jitter.opt(ii) = (1/(2*pi*optFreq)) * intPhaseNoise.rad(ii);
    
    minVal.dBc(ii) = min(dataOUT(:,2,ii));
    maxVal.dBc(ii) = max(dataOUT(:,2,ii));
    
    
end

% jitter = (1/(2*pi*carrier.freq(1))) * sqrt( 2 * trapz(dataOUT(:,1,1) , (10 .^ (dataOUT(:,2,1) ./ 10) ) ) );
% intNoise = 10 * log10( trapz(dataOUT(:,1,1) , (10 .^ (dataOUT(:,2,1) ./10) ) ) );



%%%%% Create dBm data %%%%%
for ii = 1:numFilesIn
    
    dataOUT(:,3,ii) = dataOUT(:,2,ii) - carrier.volt(ii);
    minVal.dBm(ii) = min(dataOUT(:,3,ii));
    maxVal.dBm(ii) = max(dataOUT(:,3,ii));
    
end

%%%%% Create the rad/sqrt(Hz) data %%%%%
for ii = 1:numFilesIn
    
    dataOUT(:,4,ii) = sqrt( 2 * (10 .^ (dataOUT(:,2,ii) ./ 10) ) );
    %     dataOUT(:,4,ii) = 2 * (10 .^ (dataOUT(:,2,ii) ./ 10) );
    minVal.rad(ii) = min(dataOUT(:,4,ii));
    maxVal.rad(ii) = max(dataOUT(:,4,ii));
    
end


%%%%% Create the integrated data %%%%%
for ii = 1:numFilesIn
    
    dataOUT(:,5,ii) = imag( flip( sqrt( 2 *...
        cumtrapz( flip(dataOUT(:,1,ii)) , (10 .^ ( flip(dataOUT(:,2,ii)) ./ 10) ) )...
        ) ) );
%     dataOUT(:,5,ii) = flip(sqrt( 2 *...
%         cumtrapz( dataOUT(:,1,ii) , (10 .^ ( dataOUT(:,2,ii) ./ 10) ) )...
%         ) );
    minVal.int(ii) = min(dataOUT(:,5,ii));
    maxVal.int(ii) = max(dataOUT(:,5,ii));
    
end



for ii = 1:numFilesIn
    
    disp([fileOrder{ii},' RMS jitter is: ',num2str(jitter.opt(ii)/1e-18),' as'])
    
end



%% Run this section to plot the data

numPlot = numFilesIn;
% numPlot = size(dataOUT,3);

%%%%% Individual Plots %%%%%

titles = {['Phase Noise, ',fileOrder{1}],...
    ['Phase Noise, ',fileOrder{2}],...
    ['Phase Noise, ',fileOrder{3}]};

xLabels = repmat({'Frequency (Hz)'},numFilesIn+1,1);
yLabels = repmat({'Power (dBc)'},numFilesIn+1,1);


for ii = 1:numPlot
    figure(ii);
    ax = gca;
    
    
    semilogx(dataOUT(:,1,ii),dataOUT(:,2,ii),'LineWidth',2);
    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([minVal.dBc(ii)-10 maxVal.dBc(ii)+10]);
    
    % Loop BW Line
    line([loopBW(ii) loopBW(ii)],[minVal.dBc(ii)-10 maxVal.dBc(ii)+10],...
        'Color',[0 0 0 .5],'LineStyle','--');
    format shortG
    text(loopBW(ii)*1.1,maxVal.dBc(ii)-10,...
        ['Loop BW: ', num2str(loopBW(ii)),' Hz'],...
        'FontSize',14,'Color',[0 0 0]);
    
    ax.FontSize = 24;
    
    % Make it look nice
    title(titles{ii},'FontSize',30);
    xlabel(xLabels{ii},'FontSize',24);
    ylabel(yLabels{ii},'FontSize',24);
    
    
end

%% Special Plots %%

numPlot = 1:numFilesIn;
% numPlot = 1:size(dataOUT,3);
% numPlot = 2;


% %%%%% dBc All Plot %%%%%
% titles = {'PSD of Phase Noise, All'};
% xLabels = repmat({'Frequency (Hz)'},numFilesIn+1,1);
% yLabels = repmat({'Power (dBc/Hz)'},numFilesIn+1,1);
% 
% 
% figure(99);
% ax = gca;
% 
% for ii = numPlot
%     semilogx(dataOUT(:,1,ii),dataOUT(:,2,ii),'linewidth',2);
%     hold on
%     xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
%     ylim([min(minVal.dBc)-10 max(maxVal.dBc)+10]);
% end
% 
% 
% hold off
% 
% % Make it look nice
% legend({'Locked','Clock','Unlocked'},'Location','northeast');
% ax.FontSize = 24;
% 
% title(titles{1},'FontSize',30,'Interpreter','latex');
% xlabel(xLabels{1},'FontSize',24,'Interpreter','latex');
% ylabel(yLabels{1},'FontSize',24,'Interpreter','latex');



% %%%%% rad/sqrt(Hz) w/o IPN All Plot %%%%%
% figure(98);
% ax = gca;
% 
% titles = {'PSD of Phase Noise, All'};
% xLabels = repmat({'Frequency (Hz)'},numFilesIn+1,1);
% yLabels = {'PND (rad $$\mathrm{Hz}^{-1/2}$$)'};
% 
% for ii = numPlot
%     loglog(dataOUT(:,1,ii),dataOUT(:,4,ii),'linewidth',2);
%     hold on
% end
% 
% xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
% ylim([10^-8 10^0]);
% 
% hold off
% 
% % Make it look nice
% legend(fileOrder,'Location','northeast');
% 
% % Labels
% ylabel(yLabels{1},'FontSize',24,'Interpreter','latex');
% xlabel(xLabels{1},'FontSize',24,'Interpreter','latex');
% 
% ax.FontSize = 34;
% title(titles{1},'FontSize',40,'Interpreter','latex');
% 
% 
% % Force numbers on axis
% xTicks = 10.^(0:2:6)';
% ax.XTick = xTicks;
% ax.XTickLabel = cellstr(num2str(round(log10(xTicks)), '10^%d'));
% 
% yTicks = 10.^(-8:1:0)';
% ax.YTick = yTicks;
% ax.YTickLabel = cellstr(num2str(round(log10(yTicks)), '10^%d'));




% %%%%% rad/sqrt(Hz) w/ IPN All Plot %%%%%
% figure(97);
% ax = gca;
% 
% lineTypes = {'-','-','-'};
% lineColors = {[120 137 163]/255,[33 54 86]/255,[0 0 0]/255};
% 
% titles = {'PSD of Phase Noise, All'};
% xLabels = repmat({'Frequency (Hz)'},numFilesIn+1,1);
% yLabels = {'PND (rad $$\mathrm{Hz}^{-1/2}$$)',...
%     'Integrated PND (mrad)'};
% 
% yyaxis left
% 
% for ii = numPlot
% %     loglog(dataOUT(:,1,ii),dataOUT(:,4,ii),lineTypes{ii},'linewidth',4);
%     loglog(dataOUT(:,1,ii),dataOUT(:,4,ii),lineTypes{ii},...
%         'color',lineColors{ii},...
%         'linewidth',4);
%     hold on
% end
% 
% ylabel(yLabels{1},'FontSize',24,'Interpreter','latex');
% 
% 
% yyaxis right
% 
% for ii = 2
% % for ii = numPlot
%     semilogx(dataOUT(:,1,ii),dataOUT(:,5,ii)*1e3,lineTypes{ii},'linewidth',4);
% end
% 
% 
% xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
% ylim([0 max(intPhaseNoise.rad(1:2)*1e3)*1.25]);
% % ylim([0 max(intPhaseNoise.rad(1:2)*1e3)*.0005]);
% ylabel(yLabels{2},'FontSize',24,'Interpreter','latex');
% 
% 
% hold off
% 
% 
% % Make it look nice
% % legend([fileOrder,'Integrated PND (Clock)','Integrated PND (Locked)'],'Location','northeast');
% % legend([fileOrder,'Integrated PND'],'Location','northeast');
% legend({'RF Reference','$$f_{OOL}$$','Integrated PND'},'Location','northeast','Interpreter','latex');
% % legend(fileOrder,'Location','northeast','Interpreter','latex');
% grid on
% 
% xlabel(xLabels{1},'FontSize',24,'Interpreter','latex');
% 
% ax.FontSize = 40;
% % title(titles{1},'FontSize',40,'Interpreter','latex');


%% Fancy Table with the Results %%


Range = repmat({'1 Hz to 3 MHz'},1,3);
% Range = repmat({'1 Hz to 1 MHz'},1,3);

jitPrint = { [num2str(jitter.opt(1)/1e-18),' as (Opt)'] ...
    [num2str(jitter.opt(2)/1e-18),' as (Opt)'] ...
    [num2str(jitter.opt(3)/1e-18),' as (Opt)']};


T = table(...
    fileOrder',Range',intPhaseNoise.dBc',intPhaseNoise.rad',jitPrint',...
    'VariableNames',{'State','Range','IPN_dBc','IPN_rad','JitterRMS'});
for ii = 1:3
    fprintf('\n')
end
disp(T)

%%
writetable(T,'D:\CEP\Tabulated\PhaseNoise_ULM.txt',...
    'Delimiter',',',...
    'QuoteStrings',true);






