function out = rng_wrapper(in,varargin)
% Purpose: call the Matlab or Octave specific unique function
%
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA

if strcmp(reportUI,'matlab')
  if nargin > 0 && isempty(varargin)
    rng(in);
  elseif nargin > 0
    rng(in,varargin{:});
  end
  out = rng;
else
  if nargin > 0 && isempty(varargin)
    rng_octave(in);
  elseif nargin > 0
    rng_octave(in,varargin{:});
  end
  out = rng_octave;
end
