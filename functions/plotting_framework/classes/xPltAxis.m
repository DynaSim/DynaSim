
classdef xPltAxis
    % Metadata information. This is inserted into the "meta" property of
    % xPlt class

    properties
        name = ''         % (optional) 1x1 - string naming the axis in question (e.g. time, voltage, cell number, etc)
        values = []       % (optional) 1xN - can be numerical matrix or cell array of strings describing axis in question
        astruct = struct   % (optional) 1xN structure array of additional information
    end
    
    methods
        function getaxisinfo(xpa)
            temp = [xpa.name, ' -> '];
            Nvals = length(xpa.values);
            for i = 1:Nvals-1
                
                temp = [temp,xpa.getvaluestring(i),', '];
                
            end
            temp = [temp,xpa.getvaluestring(Nvals),'\n'];
            
        
            fprintf(temp)
        end
        
        function out = getvalue_noncell(xpa,i)
            % Looks at entry xpa.value(i) and returns its output as a
            % numeric, regardless of whether it is actually cell array
            if iscell(xpa.values)
                out = xpa.values{i};
            else
                out = xpa.values(i);
            end
        end
        
        function out = getvaluestring(xpa,i)
            % Looks at entry xpa.value(i) and returns its output as a
            % string, regardless of what data type it actually is (string,
            % cell, numeric, etc).
            out = xpa.getvalue_noncell(i);
            if isnumeric(out)
                out = num2str(out);
            end
        end

        
    end
end
