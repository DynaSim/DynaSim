function filename = filepartsNameExt(file)
%% filepartsNameExt
% Purpose: Return filename with extension from path

[~, name, ext] = fileparts(file);

filename = [name, ext];

end