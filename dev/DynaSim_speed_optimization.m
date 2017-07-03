%% optimize conditionals with spike times

thresh=.5;
Npop=8e3;
spike_buffer=5; % 1, 100
tabs=0;
T=0:.01:200; t=0;
L=length(T);
V=rand(L,Npop);
tspike=rand(spike_buffer,Npop);

tic
for n=2:L, test=(V(n,:)>=thresh&V(n-1,:)<thresh); end
toc

tic
for n=2:L, test=(any(t<tspike+tabs,1)); end
toc

%{

f = @(V,tspike_pre)gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(V-ESYN)
g = @(V)0
Isyn = @(V)iff(any(t<tspike_pre+max_int,1),V,f,g)

where
function result = iff(condition,stateVar,f,g)
global tspike_pre
if condition
  result = f(stateVar,tspike_pre);
else
  result = g(stateVar);
end

%}


%% Compare sparse vs full matrix operations

p=.05; % connection probability
Xf=(rand(1000,1000)<p).*rand(1000,1000);
Xs=sparse(Xf);
af=rand(1,size(Xf,1));
as=sparse(af);

tic; Xf*Xf; toc1=toc;
tic; af*Xf; toc2=toc;
tic; as*Xf; toc3=toc;
tic; (af*Xf).*af; toc4=toc;
[toc1 toc2 toc3 toc4]

tic; Xs*Xs; toc1=toc;
tic; af*Xs; toc2=toc;
tic; as*Xs; toc3=toc;
tic; (af*Xs).*af; toc4=toc;
[toc1 toc2 toc3 toc4]
