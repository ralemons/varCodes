%% Load in data

clear

saveData = load('~/Documents/GitHub/varCodes/savedVars/nonLinVarAng_Ideal.mat');
N = size(saveData.var,2);

%% Figure 1, Independant Var vs Theta-0 Build-Up

figure(1)

p = numSubPlots(N);

yMax = max([saveData.Y{:}],[],'all')*1.1;
yMin = min([saveData.Y{:}],[],'all')*0.9;

for ii = 1:N
    subplot(p(1),p(2),ii)
    
    plot(saveData.T{ii},saveData.Y{ii})
    title([saveData.indVar,': ' num2str(saveData.var(ii))])
    ylim([yMin yMax])
    
end


%% Figure 2, Independant Var vs Sound Data

figure(2)

p = numSubPlots(N);

yMax = max([saveData.sigAmp],[],'all')*1.1;
yMin = min([saveData.sigAmp],[],'all')*0.9;


for ii = 1:N
    subplot(p(1),p(2),ii)
    
    plot(saveData.sigAmp(:,:,ii))
    title([saveData.indVar,': ' num2str(saveData.var(ii))])
    ylim([yMin yMax])
    
end


%%

figure(3)

p = numSubPlots(N);

yMax = -Inf; yMin = Inf;

for ii = 1:N
    subplot(p(1),p(2),ii)
    
    tmp = abs(fft(saveData.sigAmp(:,:,ii))).^2;
    tmp(1,:) = [];
        
    yMax = max(max(tmp,[],'all')*1.1,yMax);
    yMin = min(min(tmp,[],'all')*0.9,yMin);
    
    plot(tmp)
    title([saveData.indVar,': ' num2str(saveData.var(ii))])
    
end

for ii = 1:N
    subplot(p(1),p(2),ii)
    ylim([yMin yMax])
end


%%
function p = numSubPlots(n)

p(2) = ceil(sqrt(n));
p(1) = ceil(n/p(2));

end