function displayError(err)
%DISPLAYERROR - Wrapper to print MATLAB error
  fprintf('Error: %s\n',err.message);
  for i=1:length(err.stack)
    fprintf('\t in %s (line %g)\n',err.stack(i).name,err.stack(i).line);
  end
