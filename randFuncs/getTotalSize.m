function varargout = getTotalSize(obj)

props = properties(obj);
totSize = 0;

if isempty(props)
    totSize = whos('obj');
    totSize = totSize.bytes;
    
else
    for ii  = 1:length(props)
        s = obj.(props{ii});
        if isobject(s)
            s = getTotalSize(s);
            totSize = totSize + s;
        else
            s = whos('s');
            totSize = totSize + s.bytes;
        end
    end
end




if nargout == 0
    
    fprintf(1, '>> %.4f MiB\n', totSize/1024^2);
    
elseif nargout == 1
    
    varargout{1} = totSize;
    
end


end