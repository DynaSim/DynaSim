function mlBool = ismatlab()
  % returns logical of whether running in matlab environment
  %
  % Author: Erik Roberts

  mlBool = strcmp(reportUI,'matlab');
end