%Import the data from the simpler Spectrum Analyser
%   Has an optional input of the filename. If nothing is input then it will
%   pop up an interactive box to choose file.
function [output,varargout] = dataImport_SpecAnalyser_RS(fileStr)

if nargin == 0
    [fileName,filePath] = uigetfile('D:\CEP\DATs\*.*','MultiSelect','on');
    fileStr = fullfile(filePath,fileName);
end

singFile = whos('fileStr');

if strcmpi(singFile.class,'cell')
    for ii = 1:length(fileName)
        temp = readtable(fileStr{ii});
        temp = table2array(temp(27:end,1:2));
        output(:,:,ii) = str2double(strrep(temp,',','.'));
    end
elseif strcmpi(singFile.class,'char')
    temp = readtable(fileStr);
    temp = table2array(temp(27:end,1:2));
    output = str2double(strrep(temp,',','.'));
end

varargout{1} = fileStr;


end