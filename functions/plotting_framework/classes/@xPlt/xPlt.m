
classdef xPlt

    properties
        data        % Storing the actual data
        meta        % Cell array of xPltMeta classes for each axis
        info        % xPltMeta classes
    end
    
    methods
        function [xp2, popnames] = inds(xp,xp_inds)
            % xp - Class xPlt
            % xp_inds - Class of xPltInds
            
            xp2 = xp;
            
            all_true_flag = all(xp_inds.ind_selected == true);
            
            % Pull out relevant info
            im = xp_inds.ind_meta;
            is = xp_inds.ind_selected;
            
            % Permute the chosen meta index to the front
            ndims = length(xp.meta);
            xp2 = xp2.permute([im, 1:im-1, im+1:ndims]);
            
            % Extract chosen data from xp.data
            sz = size(xp2.data);
            xp2.data = reshape(xp2.data,[sz(1),sum(sz(2:end))]);    % Compress everything into a 2x2 matrix
            xp2.data = xp2.data(is,:);
            xp2.data = reshape(xp2.data, [size(xp2.data,1), sz(2:end)]);
            
            % Extract the chosen data from xp.meta
            xp2.meta{1} = xp2.meta{1}.inds(is);
            
            % Permute back to the original order
            xp2 = xp2.permute([2:im, 1, im+1:ndims]);

%             if all_true_flag
%                 % Everything is true so can remove the dimension entirely
%                 xp2.data = xp2.data(
%             else
%             end
        end
        
        function xp = permute(xp,order)
            xp.data = permute(xp.data,order);
            xp.meta = xp.meta(order);
        end
    end
end


function output = inheritObj(output,input)
    % Merges contents of input into output.
    C = metaclass(input);
    P = C.Properties;
    for k = 1:length(P)
        if ~P{k}.Dependent
            output.(P{k}.Name) = input.(P{k}.Name);
        end
    end
end


