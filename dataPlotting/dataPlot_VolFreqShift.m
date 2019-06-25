voltages = [1.16,.66,0,-.66,-1.16,-1.66,-1.86,...
    -2.16,-2.66,-3.16,-3.66,-4.16,-4.36,-4.96]';
frequencies = [85.70,86.55,87.116,87.683,87.683,...
    87.12,87.16,86.83,85.98,84.85,83.43,81.16,80.6,77.48]';
disp([voltages,frequencies])

figure(1);
ax1 = gca;

plot(voltages,frequencies,'b-o')

xlim([-5,1.6]);
ylim([75,90]);
xticks(fliplr(voltages'));
xtickangle(45)
yticks([75:2.5:90]);
ax1.FontSize = 16;

title('Change in CEO Frequency wrt Pump Voltage','FontSize',30);
xlabel('Pump Voltage','FontSize',24);
ylabel('CEO Frequency','FontSize',24);