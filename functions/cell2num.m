function m = cell2num(c)
%
% This function takes a cell matrix of strings and converts
% each cell into a number.

  sz = size(c);
  m = zeros(sz);
  for i=1:prod(size(c))
    m(i) = sscanf(c{i},'%f');
  end;