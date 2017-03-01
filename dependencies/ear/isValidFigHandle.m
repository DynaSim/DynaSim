function validBool = isValidFigHandle(h)
%% isValidFigHandle
% Purpose: check if input exists and is a valid figure handle.
%
% Usage:
%   validBool = isValidFigHandle(handle)
%   validBool = isValidFigHandle(strForHandle)
%
% Input: string with handle name, or handle itself
%
% Note: will work with string containing any callable matlab indexing/object

try
  if ischar(h) %string input
    validBool = evalin('caller', ['isvalid(' h ')']);
  else
    validBool = isvalid(h);
  end
catch
  validBool = 0;
end

end