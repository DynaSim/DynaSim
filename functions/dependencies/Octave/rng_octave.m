function out = rng_octave(x)
%rng_octave(x) establishes the seed and reports seed and state for the random number generator of rand, rande, randg, randn and randp
%x must be:
%- a numeric value that sets the seed value of the rng
%- 'shuffle': it seeds rng based on current time
%- 'default': it sets the seed = 0 (actually 5489 for Matlab compatibility)
%if x is empty or missing, rng returns the seed of the random number generator

  errmsg = 'only a single numeric or limited string argument is supported, type help rng_octave for details';

  if nargin == 0
    x = [];
  elseif nargin > 1
    error(errmsg);
  elseif ischar(x)
    if strcmp(x,'shuffle')
      tnow = clock;
      x = 1e6*tnow(end);
    elseif strcmp(x,'default')
      x = 0;
    else
      error(errmsg);
    end
  elseif ~isnumeric(x)
    error(errmsg);
  end

  if isempty(x)
    rand_seed = uint32(rand('seed'));
    rande_seed = uint32(rande('seed'));
    randg_seed = uint32(randg('seed'));
    randn_seed = uint32(randn('seed'));
    randp_seed = uint32(randp('seed'));
    if isequal(rand_seed,rande_seed,randg_seed,randn_seed,randp_seed)
      out.Seed = rand_seed;
    else
      out.rand_seed = rand_seed;
      out.rande_seed = rande_seed;
      out.randg_seed = randg_seed;
      out.randn_seed = randn_seed;
      out.randp_seed = randp_seed;
    end
    rand_state = uint32(rand('state'));
    rande_state = uint32(rande('state'));
    randg_state = uint32(randg('state'));
    randn_state = uint32(randn('state'));
    randp_state = uint32(randp('state'));
    if isequal(rand_state,rande_state,randg_state,randn_state,randp_state)
      out.State = rand_state;
    else
      out.rand_state = rand_state;
      out.rande_state = rande_state;
      out.randg_state = randg_state;
      out.randn_state = randn_state;
      out.randp_state = randp_state;
    end
  else
    out.Seed = uint32(x);
    out.State = twister_state(out.Seed);
    rand('seed',out.Seed);
    rande('seed',out.Seed);
    randg('seed',out.Seed);
    randn('seed',out.Seed);
    randp('seed',out.Seed);
    rand('state',out.State);
    rande('state',out.State);
    randg('state',out.State);
    randn('state',out.State);
    randp('state',out.State);
  end
  out.Type = 'twister';

  % states in Octave and Matlab are different, this fixes the difference for rand
  % current limitations:
  % - unfortunately it does not seem to work for matching randn
  % - in Matlab, the same rng applies to rand and randn, whereas in Octave each one evolves in its way. This means that in Matlab rand, randn and randn, rand give different series, whereas in Octave they are the same
  function state = twister_state(seed)
    if nargin < 1 || seed == 0
      seed = 5489; % for Matlab compatibility, 0 and 'default' actually corresponds to seed 5489 in Matlab
    end
    state = uint32(zeros(625,1));
    state(1) = seed;
    for N = 1:623
        %% initialize_generator
        % bit-xor (right shift by 30 bits)
        uint64(1812433253)*uint64(bitxor(state(N),bitshift(state(N),-30)))+N; % has to be uint64, otherwise in 4th iteration hit maximum of uint32!
        state(N+1) = uint32(bitand(ans,uint64(intmax('uint32')))); % untempered numbers
    end
    state(end) = 1; % Matlab displays 624 here, but for compatibility, this has to be set to 1
  end
end
