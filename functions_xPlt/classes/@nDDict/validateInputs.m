


function outputformats = validateInputs(xp,field,field_type)
    % For internal use in importing data in to nDDict.

    switch field_type
        case 'data'
            if iscell(field)
                outputformats = 'cell';          
            elseif isnumeric(field)
                outputformats = 'numeric';
        %     elseif isstring(X)
        %         outputformat = 'string';
            else
                error('Input obj.data must be cell or numeric');
            end


        case 'axis'
            outputformats = cell(1,length(field));
            for i = 1:length(field)
                outputformats{i} = validate_axis(field{i});
                if strcmp(outputformats{i},'unknown')
                    error(['axislabels' num2str(i) ' must be of type numeric or cell array of character vectors']);
                end
            end
        otherwise
            error('Unrecognized input foramt');

    end
    
    
    function ax_type = validate_axis(axlinear)

        if isnumeric(axlinear)
            ax_type = 'numeric';
    %     elseif isstring(axlinear)
    %         ax_type = 'string';
        elseif iscell(axlinear)
            if is_cellarray_of_chars(axlinear)
                ax_type = 'cell_of_chars';
            else
                ax_type = 'unknown';
            end
        else
            ax_type = 'unknown';
        end

    end


    function out = is_cellarray_of_chars(c)
        out = all(cellfun(@ischar,c));
    end


end
