function dsUnitSave_figHandles( handles, outputDir )
%SAVE_FIGHANDLES saves fig handles to outputDir

for iH = 1:length(handles)
  thisFigName = sprintf('fig%i',iH);
  savefig(handles(iH), fullfile(outputDir, thisFigName));
    % NOTE: compareFigFiles won't work with compact mode for savefig
end

end
