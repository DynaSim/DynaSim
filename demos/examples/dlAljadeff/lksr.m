function r = lksr(x,y,h,N)
% lksr Laplacian kernel smoothing regression
%
% r = lksr(x,y) returns the Laplacian kernel regression in structure r such that
%   r.f(r.x) = y(x) + e
% The bandwidth and number of samples are also stored in r.h and r.n
% respectively.
%
% r = lksr(x,y,h) performs the regression using the specified bandwidth, h.
%
% r = lksr(x,y,h,n) calculates the regression in n points (default same as input).
%
% Without output, lksr(x,y) or lksr(x,y,h) will display the regression plot.
%
% Algorithm
% The kernel regression is a non-parametric approach to estimate the
% conditional expectation of a random variable:
%
% E(Y|X) = f(X)
%
% where f is a non-parametric function. Based on the kernel density
% estimation, this code implements the Nadaraya-Watson kernel regression
% using the Laplacian kernel as follows:
%
% f(x) = sum(kerf((x-X)/h).*Y)/sum(kerf((x-X)/h))
%
% See also gksr, eksr, gkde, ksdensity

% Example 1: smooth curve with noise
%{
x = 1:100;
y = sin(x/10)+(x/50).^2;
yn = y + 0.2*randn(1,100);
r = lksr(x,yn);
plot(x,y,'b-',x,yn,'co',r.x,r.f,'r--','linewidth',2)
legend('true','data','regression','location','northwest');
title('Laplacian kernel regression')
%}
% Example 2: with missing data
%{
x = sort(rand(1,100)*99)+1;
y = sin(x/10)+(x/50).^2;
y(round(rand(1,20)*100)) = NaN;
yn = y + 0.2*randn(1,100);
r = lksr(x,yn);
plot(x,y,'b-',x,yn,'co',r.x,r.f,'r--','linewidth',2)
legend('true','data','regression','location','northwest');
title('Laplacian kernel regression with 20% missing data')
%}

% By Yi Cao at Cranfield University on 12 March 2008.
%

% Check input and output
narginchk(2,4);
nargoutchk(0,1);
if numel(x) ~= numel(y)
  error('x and y are in different sizes.');
end

x = x(:);
y = y(:);
% clean missing or invalid data points
inv = (x~=x) | (y~=y);
x(inv) = [];
y(inv) = [];

% Default parameters
if nargin < 4
  % N = 100;
  N = length(x);
elseif ~isscalar(N)
  error('N must be a scalar.')
end
r.n = length(x);
if nargin < 3
  % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
  hx = median(abs(x-median(x)))/0.6745*(4/3/r.n)^0.2;
  hy = median(abs(y-median(y)))/0.6745*(4/3/r.n)^0.2;
  if ~hx && ~hy
    h = sqrt(hy*hx);
  elseif ~hy
    h = hy;
  else
    h = hx;
  end
end
r.h = h;

% Laplacian kernel function
kerf = @(z)exp(-abs(z));

r.x = linspace(min(x),max(x),N);
r.f = zeros(1,N);

% for k = 1:N
%   z = kerf((r.x(k)-x)/h);
%   % figure, plot(z)
%   r.f(k) = sum(z.*y)/sum(z);
% end

% vectorized
z = kerf((r.x-x)./h');
r.f = sum(z.*y,1)./sum(z,1);

% Plot
if ~nargout
  figure
  plot(r.x,r.f,'r',x,y,'bo')
  hold on
  ylabel('f(x)')
  xlabel('x')
  title('Kernel Smoothing Regression');
end
