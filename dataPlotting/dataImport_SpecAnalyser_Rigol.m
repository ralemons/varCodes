%Import the data from the simpler Spectrum Analyser
%   Has an optional input of the filename. If nothing is input then it will
%   pop up an interactive box to choose file.
function [output,varargout] = dataImport_SpecAnalyser(fileStr)

if nargin == 0
    [fileName,filePath] = uigetfile('D:\CEP\CSVs\*.*','MultiSelect','on');
    fileStr = fullfile(filePath,fileName);
end

singFile = whos('fileStr');

if strcmpi(singFile.class,'cell')
    for ii = 1:length(fileName)
        temp = readtable(fileStr{ii});
        output(:,:,ii) = table2array(temp(:,[1 3]));
    end
elseif strcmpi(singFile.class,'char')
    temp = readtable(fileStr);
    output = table2array(temp(:,[1 3]));
end

varargout{1} = fileStr;


end