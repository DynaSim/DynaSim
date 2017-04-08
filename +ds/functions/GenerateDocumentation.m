function GenerateDocumentation
%GENERATEDOCUMENTATION - Build DynaSim function documentation
%
%  This is a simple call to `m2html()` to make it easy for anyone with DynaSim
%  on their path to build the documentation.
%
%  HOWEVER, there are caveats:
%  1. This MUST be run from the root DynaSim directory, meaning the current
%     directory of your MATLAB instance must be in your "/path/to/dynasim".
%    - This is because `m2html` must be clear on what "functions" folder you are
%      referring to, and because `function` and `functions` are reserved words in
%      MATLAB, you can't just simply locate the folder globally using
%      `which(functions)`. There may be a way around this, but you'll have to
%      step very carefully around the reserved words.
%  2. This is a custom version of m2html, in that I replaced all deprecated
%     'error(nargin(...' uses, which were giving warnings that slowed down the
%     program, into modern 'narginchk' uses, so this is NOT identical to the
%     downloadable m2html.

m2html('mfiles',{'functions','functions_xPlt'},...
       'htmldir','docs',...
       'recursive','on',...
       'global','on',...
       'template','frame',...
       'index','menu',...
       'graph','on');
