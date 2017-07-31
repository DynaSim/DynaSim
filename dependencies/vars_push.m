
% Pushes all variables in the current workspace to the strcuture s
function s = vars_push()

% Grab workspace
w = evalin('caller', 'whos');
names = {w.name};
s = struct;
for i = 1:numel(w)
    s.(names{i}) = evalin('caller', names{i});
end
