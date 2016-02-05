function str=dynasim_strrep(str,oldstr,newstr,lpad,rpad)
% Purpose: replace full words by new character strings, ignoring matches
% that appear as sub-strings. note: built-in strrep replaces all matches.
% Examples:
% dynasim_strrep('(v)*(-av)','v','pop1_v')
% dynasim_strrep('v-v^2+vav','v','pop1_v')
% dynasim_strrep('v-v-v','v','pop1_v')
% dynasim_strrep('v-v-v^2','v','pop1_v')
% dynasim_strrep('(v-v-v^2)','v','pop1_v')
% dynasim_strrep('E-pop1_V+1','pop1_V','pop1_V(n-1)')
% dynasim_strrep('v=1; u=u+d','u','u(n,test)')

if nargin<4, lpad=''; end
if nargin<5, rpad=''; end
  if isempty(str)
    return;
  end
  pat=['([^\w\.]{1})' oldstr '([^\w]{1})']; % in the middle
    % note: exclude .oldstr for case where prefix has already been prepended and oldstr appears >1x in str
  rep=['$1' lpad newstr rpad '$2'];
  str=regexprep(str,pat,rep);
  % check for neighboring occurrence that wasn't substituted (e.g., (v-v^2) -> (pop1_v-v^2))
  % note: this is only a possible issue for strings "in the middle"
  test=['([^\w\.]{1})' oldstr '([^\w(' newstr ')]{1})'];
  if ~isempty(regexp(str,test,'once'))%~isempty(regexp(str,test,'match'))
    % substitute remaining occurrences
    str=regexprep(str,test,rep);
  end
  pat=['([^\w\.]{1})' oldstr '$'];    % at the end
  rep=['$1' lpad newstr rpad];
  str=regexprep(str,pat,rep);
  pat=['^' oldstr '([^\w]{1})'];      % at the beginning
  rep=[lpad newstr rpad '$1'];
  str=regexprep(str,pat,rep);
  pat=['^' oldstr '$'];               % all there is
  rep=[lpad newstr rpad];
  str=regexprep(str,pat,rep);
