%Import the data from the MadgeTech Enviromental Monitor
%   Has an optional input of the filename. If nothing is input then it will
%   pop up an interactive box to choose file.
function [output,varargout] = dataImport_EnvMon(fileStr)


if nargin == 0
    [fileName,filePath] = uigetfile('D:\EnviroMon\*.*');
    fileStr = fullfile(filePath,fileName);
end

% The files that the comes from the MadgeTech cloud has a super stupid 14
% line header that has a different number of columns each line. This idiot
% move on their part means we have to work around it. So we get no useful
% header info because they designed to look pretty in Excel.
output = readtable(fileStr,'HeaderLines',16);

varargout{1} = fileStr;


end