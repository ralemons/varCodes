%Import the data from the fancy Optical Spectrum Analyser
%   Has an optional input of the filename. If nothing is input then it will
%   pop up an interactive box to choose file.
function [output,varargout] = dataImport_OOS(fileStr)


if nargin == 0
    [fileName,filePath] = uigetfile('D:\LCLS-II Injector\TXTs\*.*','MultiSelect','on');
    fileStr = fullfile(filePath,fileName);
end

singFile = whos('fileStr');

if strcmpi(singFile.class,'cell')
    for ii = 1:length(fileName)
        temp = table2array(readtable(fileStr{ii}));
        output(:,:,ii) = temp;
    end
elseif strcmpi(singFile.class,'char')
    output = readtable(fileStr,'HeaderLines',34);
    output = table2array(output);
end


varargout{1} = fileStr;


end