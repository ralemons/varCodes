function out1 = bboShift1(lam1,lam2,lam3,theta,phi1,phi2,l)
%BBOSHIFT1
%    OUT1 = BBOSHIFT1(LAM1,LAM2,LAM3,THETA,PHI1,PHI2,L)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    05-Dec-2019 10:58:57

t2 = lam1.^2;
t3 = lam2.^2;
t6 = (theta.*pi)./1.8e2;
t4 = cos(t6);
t5 = lam3.^2;
t7 = sin(t6);
t8 = t5.*1.0e12;
out1 = -phi2-(l.*pi.*sqrt(1.0./(t4.^2./(t5.*-1.354e10+1.878e-2./(t8-1.822e-2)+2.7359)+t7.^2./(t5.*-1.516e10+1.224e-2./(t8-1.667e-2)+2.3753))).*2.0)./lam3+(l.*pi.*sqrt(t2.*-1.354e10+1.878e-2./(t2.*1.0e12-1.822e-2)+2.7359).*2.0)./lam1+(l.*pi.*sqrt(t3.*-1.354e10+1.878e-2./(t3.*1.0e12-1.822e-2)+2.7359).*2.0)./lam2;