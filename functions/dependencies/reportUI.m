function ui = reportUI
% function that reports ui, that is whether Octave or matlab is in use
if exist('OCTAVE_VERSION', 'builtin')
  ui = 'octave';
else
  ui = 'matlab';
end
