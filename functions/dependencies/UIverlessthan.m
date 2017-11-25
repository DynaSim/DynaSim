function result = UIverlessthan (ref_version);
UIversion = version;
if strcmp(reportUI,'matlab')
  UIversion = UIversion(1:5);
end
UIversion = hex2dec(UIversion(1:2:end));
ref_version = hex2dec(ref_version(1:2:end));
result = UIversion < ref_version;
