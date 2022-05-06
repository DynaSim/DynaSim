function dsParamsModifier(tempfuncname, map)

    fileID = fopen(tempfuncname, 'w');
    fprintf(fileID, 'function dlTempFuncParamsChanger(dlPath)\n\n');
    fprintf(fileID, '\tp = load([dlPath, ''/params.mat'']);\n\n');
    n = size(map);
    
    labels = map.keys();
    values = map.values();
    
    for i = 1:n
    
        fprintf(fileID, '\tif sum(strcmpi(fieldnames(p.p), ''%s''))\n', labels{1, i});
        if strcmpi(class(values{1, i}), 'double')
             
            m = max(size(values{1, i}));
            if m == 1
                
                fprintf(fileID, '\t\tp.p.%s = %d;\n', labels{1, i}, values{1, i});
            
            else
                
                x = values{1, i};
                m = size(x, 1);
                l = size(x, 2);
                
                fprintf(fileID, '\t\tp.p.%s = [', labels{1, i});
                for j = 1:m
                    
                    if j > 1                    
                        fprintf(fileID, ';');
                    end
                    
                    for k = 1:l
                        fprintf(fileID, ' %d', x(j, k));
                    end
                    
                end
                fprintf(fileID, '];\n');
                
            end
            
        elseif ischar(values{1, i})
            
            fprintf(fileID, '\t\tp.p.%s = ''%s'';\n', labels{1, i}, values{1, i});
        
        else
            
            fprintf(fileID, '\t\tp.p.%s = %g;\n', labels{1, i}, values{1, i});
        
        end
        
        fprintf(fileID, '\telse\n');
        fprintf(fileID, '\t\tfprintf("Parameter or variable ''%s'' not found in params.mat file. Check if you are refering to a correct variable.\\n");\n', labels{1, i});
        fprintf(fileID, '\tend\n\n');
        
    end
    
    fprintf(fileID, '\tsave([dlPath, ''/params.mat''], ''-struct'', ''p'');\n\nend');
    fclose(fileID);

end