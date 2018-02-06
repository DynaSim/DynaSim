function totalGB = memoryUsageCallerGB()
  vars = evalin('caller','whos');
  totalBytes = sum([vars.bytes]);
  
  totalGB = totalBytes/1073741824;
end