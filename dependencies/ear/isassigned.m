function returnBool = isassigned(varStr)
%% isassigned
% Author: Erik Roberts
%
% Purpose: Check if variable exists and is not empty
%
% Usage: returnBool = isassigned(varStr)

  returnBool = evalin('caller', ['exist(''' varStr ''', ''var'')']) && ~isempty(evalin('caller', varStr));
  
end