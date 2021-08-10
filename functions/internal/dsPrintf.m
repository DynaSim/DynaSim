function dsPrintf(options, varargin)
%% dsPrintf
% printf if options.verbose_flag

if options.verbose_flag
  fprintf(varargin{:});
end

end