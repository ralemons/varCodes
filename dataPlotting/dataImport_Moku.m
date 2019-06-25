%Import the data from the Moku Lab Spectrum Analyzer
%   Has an optional input of the filename. If nothing is input then it will
%   pop up an interactive box to choose file.
function [output,varargout] = dataImport_Moku(fileStr)


if nargin == 0
    [fileName,filePath] = uigetfile('D:\CEP\CSVs\*.*','MultiSelect','on');
    fileStr = fullfile(filePath,fileName);
end

singFile = whos('fileStr');

if strcmpi(singFile.class,'cell')
    for ii = 1:length(fileStr)
        temp = importdata(fileStr{ii});
        output(:,:,ii) = temp.data;
    end
elseif strcmpi(singFile.class,'char')
    temp = importdata(fileStr);
    output = temp.data;
end


varargout{1} = fileStr;

end