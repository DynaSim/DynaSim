function out = rng_wrapper(in,varargin)
% Purpose: call the Matlab or Octave specific unique function
%
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA

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
    [rand_seed, randn_seed] = rng_octave;
  elseif isempty(varargin)
    [rand_seed, randn_seed] = rng_octave(in);
  else
    [rand_seed, randn_seed] = rng_octave(in,varargin{:});
  end
  out = [rand_seed, randn_seed];
end
