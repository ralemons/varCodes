%Import the data from the simpler Spectrum Analyser
%   Has an optional input of the filename. If nothing is input then it will
%   pop up an interactive box to choose file.
function [output,header,varargout] = dataImport_PhaseNoise(fileStr)

if nargin == 0
    [fileName,filePath] = uigetfile('D:\CEP\DATs\*.*','MultiSelect','on');
    fileStr = fullfile(filePath,fileName);
end

singFile = whos('fileStr');

if strcmpi(singFile.class,'cell')
    for ii = 1:length(fileName)
        temp = readtable(fileStr{ii});
        firstNum = find(~cellfun(@isempty,regexp(temp{:,1},'^\d','match','once')),1);
        header(:,:,ii) = table2array(temp(1:firstNum-1,:));
        temp = table2array(temp(firstNum:end,1:2));
        output(:,:,ii) = str2double(strrep(temp,',','.'));
    end
elseif strcmpi(singFile.class,'char')
    temp = readtable(fileStr);
    firstNum = find(~cellfun(@isempty,regexp(temp{:,1},'^\d','match','once')),1);
    header = table2array(temp(1:firstNum-1,:));
    temp = table2array(temp(firstNum:end,1:2));
    output = str2double(strrep(temp,',','.'));
end

varargout{1} = fileStr;


end