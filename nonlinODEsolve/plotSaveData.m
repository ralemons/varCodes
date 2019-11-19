%% Load in data

clear
loadVar = 1;

if loadVar == 1
    
    load('~/Documents/GitHub/varCodes/savedVars/nonLinVarAmp.mat');
    N = size(saveData.energy,2);
    
else


    saveData = load('~/Documents/GitHub/varCodes/savedVars/nonLinVarAng.mat');
    N = size(saveData.var,2);


end

%% Figure 1, Independant Var vs Theta-0 Build-Up

figure(1)

if loadVar == 1
    
    
    p = numSubPlots(N);
    
    for ii = 1:N
        subplot(p(1),p(2),ii)
        
        plot(saveData.T{ii},saveData.Y{ii})
        title(['energy: ', num2str(saveData.energy(ii,1)/10^-9), 'nJ' ])
        
    end
    
    
else
    
    
    p = numSubPlots(N);
    
    for ii = 1:N
        subplot(p(1),p(2),ii)
        
        plot(saveData.T{ii},saveData.Y{ii})
        title([saveData.indVar,': ' num2str(saveData.var(ii))])
        
    end
    
    
end


%% Figure 2, Independant Var vs Sound Data

figure(2)

if loadVar == 1
    
    
    p = numSubPlots(N);
    
    for ii = 1:N
        subplot(p(1),p(2),ii)
        
        plot(saveData.songAmp(:,ii))
        title(['energy: ', num2str(saveData.energy(ii,1)/10^-9), 'nJ' ])
        
    end
    
    
else
    
    
    p = numSubPlots(N);
    
    for ii = 1:N
        subplot(p(1),p(2),ii)
        
        plot(saveData.songAmp(:,ii))
        title([saveData.indVar,': ' num2str(saveData.var(ii))])
        
    end
    
    
end


%%

figure(3)

if loadVar == 1
    
    
    p = numSubPlots(N);
    
    for ii = 1:N
        subplot(p(1),p(2),ii)
        
        tmp = abs(fft(saveData.songAmp(:,ii))).^2;
        tmp(1) = [];        
        x = 1:1:4000;
        
        
        plot(x,tmp(x))
        title(['energy: ', num2str(saveData.energy(ii,1)/10^-9), 'nJ' ])
        
    end
    
    
else
    
    
    p = numSubPlots(N);
    
    for ii = 1:N
        subplot(p(1),p(2),ii)
        
        tmp = abs(fft(saveData.songAmp(:,ii))).^2;
        tmp(1) = [];        
        x = 1:1:4000;
        
        
        plot(x,tmp(x))
        title([saveData.indVar,': ' num2str(saveData.var(ii))])
        
    end
    
    
end



%%
function p = numSubPlots(n)

p(2) = ceil(sqrt(n));
p(1) = ceil(n/p(2));

end