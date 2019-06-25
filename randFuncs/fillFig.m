function [] = fillFig(uFillW,uFillH)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-uFillW;
ax_height = outerpos(4) - ti(2) - ti(4)-uFillH;
ax.Position = [left bottom ax_width ax_height];


end