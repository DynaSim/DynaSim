function [outputs, output_variables] = dsGetOutputList(spec)

    model = dsGenerateModel(spec);
    output_variables = cat(2, 'time', model.state_variables);
    
    if ~isempty(model.monitors)
    
        output_variables = cat(2, output_variables, fieldnames(model.monitors)');
    
    end
    
    if ~isempty(model.fixed_variables)
        
        fields=fieldnames(model.fixed_variables)';
        output_variables=cat(2,output_variables,fields);
        outputs=cell(1,length(output_variables));
    
    end

end

