
function varargout = recursivePlot(xp,function_handles,dimensions)
    
    if length(function_handles) > 1
        switch dimensions{1}
            case 1
                xp2 = xPlt;
                for i = 1:size(xp.data,1)
                    xp_temp = xPlt(xp.data(i,:,:,:,:,:));
                    xp2.data{i} = @() recursivePlot(squeeze(xp_temp),function_handles(2:end),dimensions(2:end));
                end
                [varargout{1:nargout}] = function_handles{1}(xp2);
            case 2

        end
    else
        [varargout{1:nargout}] = function_handles{1}(xp);
    end
end
            
        