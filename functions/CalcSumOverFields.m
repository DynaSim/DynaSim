


function data = CalcSumOverFields(data, fields, output_field_name)
    % Purpose: Creates a new field in data that is the sum of a bunch of
    % other fields (specified by the "fields" cell array). Stores the output
    % in "output_field_name." Useful for adding multiple ionic currents
    % together.
    % Inputs:
    %   data - DynaSim data structure (see CheckData)
    %   fields - cell array of field names. These will be summed and stored
    %       in the field called output_field_name.
    %   output_field_name - string containing the desired output field
    %   name.
    % Outputs:
    %   data: data structure containing summed data

    
    if nargin < 3
        % Default field name
        output_field_name = 'summed';
        
        % Add default prefix
        ind = strfind(fields{1},'_');
        prefix = fields{1}(1:ind);
        output_field_name = strcat(prefix,output_field_name);
    end
    
    % Add prefix if necessary
    if isempty(strfind(output_field_name,'_'))
        warning('No population prefix specified (see documentation). Adding a default prefix.');
        ind = strfind(fields{1},'_');
        prefix = fields{1}(1:ind);
        output_field_name = strcat(prefix,output_field_name);
        
    end
    
    data = CheckData(data);
    % note: calling CheckData() at beginning enables analysis function to
    % accept data matrix [time x cells] in addition to DynaSim data structure.

    for i = 1:length(data)
        dat = data(i);
        sumfields = zeros(size(dat.(fields{1})));
        for j = 1:length(fields)
            sumfields = sumfields + dat.(fields{j});
        end
        data(i).(output_field_name) = sumfields;
        data(i).labels(end:end+1) = {output_field_name, data(i).labels{end}};   % Store new label as 2nd last entry (so time is always last; not sure if this is important...
    end

end