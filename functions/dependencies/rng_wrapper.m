function out = rng_wrapper(in,varargin)
% Purpose: call the Matlab or Octave specific unique function

if strcmp(reportUI,'matlab')
  if nargin==0
    out = rng;
  elseif isempty(varargin)
    out = rng(in);
  else
    out = rng(in,varargin{:});
  end
else
  if nargin==0
    out = rng_octave;
  elseif isempty(varargin)
    out = rng_octave(in);
  else
    out = rng_octave(in,varargin{:});
  end
end
