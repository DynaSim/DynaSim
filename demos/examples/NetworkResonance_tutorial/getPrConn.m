function K=getPrConn(Npre,Npost,prob,seed,zero_diag)
% get probabilistic connectivity from Npre to Npost w/ Pr[src->dst]=prob
% every target cell receives inputs from exactly prob*Npre random cells in 
% the source population. default prob = .6 (60%)
% 
% usage:
% mods={'E1->E2','netcon','getPrConn(Npre,Npost,prob)'};
% spec=ApplyModifications(base,mods);

if nargin<5, zero_diag=0; end
if nargin<4, seed=getfield(rng,'Seed'); end
if nargin<3, prob=.6; end
if nargin<2, Npost=1; end
if nargin<1, Npre=1; end

% set random seed
rng(seed);

% create probabilistic connectivity matrix
num_sources=ceil(Npre*prob);
K=zeros(Npre,Npost); % E->E, [N_pre x N_post]
for i=1:Npost
  source_list=randperm(Npre); % guarantees each source has at most 1 connection to a given target cell
  if zero_diag==1
    % exclude connections to self
    viable=~ismember(source_list,i);
  else
    % allow all sources
    viable=ones(size(source_list))==1; % logical type for consistency
  end
  viable_list=source_list(viable);
  sources=viable_list(1:min(num_sources,numel(viable_list)));
  %sources=source_list(1:num_sources);
  K(sources,i)=1;
end

% normalize kernels by number of presynaptic connections (num_sources)
K=K./repmat(max(1,sum(K,1)),[size(K,1) 1]);
% K=K/max(1,num_sources);
