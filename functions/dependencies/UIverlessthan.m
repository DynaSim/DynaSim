function result = UIverlessthan (ref_version);
  UIversion = version;
  [UIversion,UIversion2] = versionReplaceDots(UIversion);
  UIversion = hex2dec(UIversion);
  [ref_version,ref_version2] = versionReplaceDots(ref_version);
  ref_version = hex2dec(ref_version);
  result = UIversion < ref_version;
  if ~result
    UIversion2 = hex2dec(UIversion2);
    ref_version2 = hex2dec(ref_version2);
    result = UIversion2 < ref_version2;
  end
end

function [versionNoDots,versionNoDots2] = versionReplaceDots(versionDots)
  versionDots_parts = strsplit(versionDots,'.');
  versionNoDots2 = '0';
  switch numel(versionDots_parts)
    case 1, versionNoDots = [[versionDots_parts{:}],'00'];
    case 2, versionNoDots = [[versionDots_parts{:}],'0'];
    case 3, versionNoDots = [versionDots_parts{:}];
    case 4,
      versionNoDots = [versionDots_parts{1:3}];
      tmp = strsplit(versionDots_parts{4},' ');
      versionNoDots2 = [tmp{1}];
    otherwise error('version wrong format');
  end
end
