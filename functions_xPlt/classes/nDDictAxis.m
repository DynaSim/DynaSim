
classdef nDDictAxis
    % Metadata information. This is inserted into the "meta" property of
    % xPlt class

    properties
        name = ''         % (optional) 1x1 - string naming the axis in question (e.g. time, voltage, cell number, etc)
        values = []       % (optional) 1xN - can be numerical matrix or cell array of strings describing axis in question
        astruct = struct   % (optional) 1xN structure array of additional information
    end
    
    methods
        function out = getaxisinfo(obj)
            max_values_to_display = 10;
            
            temp = [obj.name, ' -> '];
            Nvals = length(obj.values);
            Nvals = min(Nvals,max_values_to_display);          % Limit to displaying 10 values
            for i = 1:Nvals-1
                temp = [temp,obj.getvaluestring(i),', '];
            end
            
            if length(obj.values) > max_values_to_display
                temp = [temp,obj.getvaluestring(Nvals),', ...'];
            else
                temp = [temp,obj.getvaluestring(Nvals)];
            end
            
            if nargout > 0
                out = temp;
            else
                fprintf([temp, '\n'])
            end
        end
        
        function out = getvaluenoncell(obj,i)
            % Looks at entry obj.value(i) and returns its output as a
            % numeric, regardless of whether it is actually cell array
            
            if length(i) > 1; error('i must be singleton'); end
            if iscell(obj.values)
                out = obj.values{i};
            else
                out = obj.values(i);
            end
        end
        
        function out = getvaluestring(obj,i)
            % Looks at entry obj.value(i) and returns its output as a
            % string, regardless of what data type it actually is (string,
            % cell, numeric, etc).
            if length(i) > 1; error('i must be singleton'); end
            out = obj.getvaluenoncell(i);
            if isnumeric(out)
                out = num2str(out);
            end
        end
        
        function out = getvaluescellstring(obj)
            % Looks at entry obj.value(i) and returns its output as a
            % cell array of strings
            out = cell(1,length(obj.values));
            for i = 1:length(obj.values)
                out{i} = num2str(obj.getvaluenoncell(i));
            end
            
            
        end

        
    end
end
