

function [dt_formatted, axn_formatted] = previewTable(data_table,ax_names)

    % All data, formatted into cells
    dt_formatted = data_table;
    for i = 1:length(dt_formatted)
        if ~iscell(dt_formatted{i}); dt_formatted{i} = num2cell(dt_formatted{i});
        end
    end
    dt_formatted = horzcat(dt_formatted{:});
    
    % Axis names
    axn_formatted = {'data', ax_names{:}};
    
    % Build a divider
    divider = repmat({'------'},1,length(axn_formatted));
    
    display(vertcat(axn_formatted,divider,dt_formatted));
%     display(dt_formatted)
    

end
    