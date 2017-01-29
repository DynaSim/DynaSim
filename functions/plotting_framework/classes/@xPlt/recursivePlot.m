
function varargout = recursivePlot(xp,function_handles,dimensions)
    

    sz = size(xp);
    if length(function_handles) > 1
        
        xp2 = xPlt;
        selection_curr = cell(1,length(sz));
        % Just cardcode in the various cases for dimensionality. It's
        % unlikely they will ever want to go above showing 4D in a single
        % plot.
        switch length(dimensions{1})        
            case 1                          % 1D
                dim1 = dimensions{1}(1);
                for i = 1:sz(dim1)
                    selection_curr{dim1} = i;
                    xp2.data{i} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end));
                end
                
            case 2                          % 2D
                dim1 = dimensions{1}(1);
                dim2 = dimensions{1}(2);
                for i = 1:sz(dim1)
                    for j = 1:sz(dim2)
                        selection_curr{dim1} = i;
                        selection_curr{dim2} = j;
                        xp2.data{i} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end));
                    end
                end
                
            case 3                          % 3D
                dim1 = dimensions{1}(1);
                dim2 = dimensions{1}(2);
                dim3 = dimensions{1}(3);
                for i = 1:sz(dim1)
                    for j = 1:sz(dim2)
                        for k = 1:sz(dim3)
                            selection_curr{dim1} = i;
                            selection_curr{dim2} = j;
                            selection_curr{dim3} = k;
                            xp2.data{i} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end));
                        end
                    end
                end
                
            case 4                          % 4D
                dim1 = dimensions{1}(1);
                dim2 = dimensions{1}(2);
                dim3 = dimensions{1}(3);
                dim4 = dimensions{1}(4);
                for i = 1:sz(dim1)
                    for j = 1:sz(dim2)
                        for k = 1:sz(dim3)
                            for l = 1:sz(dim4)
                                selection_curr{dim1} = i;
                                selection_curr{dim2} = j;
                                selection_curr{dim3} = k;
                                selection_curr{dim4} = l;
                                xp2.data{i} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end));
                            end
                        end
                    end
                end
                
                
        end
        
        % Update axes
        for i = 1:length(dimensions{1})
            dim_curr = dimensions{1}(i);
            xp2.axis(i).values = xp.axis(dim_curr).values;
            xp2.axis(i).name = xp.axis(dim_curr).name;
        end
        
        [varargout{1:nargout}] = function_handles{1}(xp2);
        
    else
        [varargout{1:nargout}] = function_handles{1}(xp);
    end
end
            
        