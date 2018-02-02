function dsVprintf(options, varargin)

if options.verbose_flag
  fprintf(varargin{:});
end

end