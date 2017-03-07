

function [dt_formatted, axn_formatted] = previewTable(data_table,col_names)

    % All data, formatted into cells
    dt_formatted = data_table;
    for i = 1:length(dt_formatted)
        if ~iscell(dt_formatted{i}); dt_formatted{i} = num2cell(dt_formatted{i});
        end
    end
    dt_formatted = horzcat(dt_formatted{:});
    
    % Build a divider
    divider = repmat({'------'},1,length(col_names));
    
    % Display the concatenated table
    display(vertcat(col_names,divider,dt_formatted));
    

end
    