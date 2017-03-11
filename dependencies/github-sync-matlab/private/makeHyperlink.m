function hyperlink = makeHyperlink(text, codeToRun)
codeToRun = ['matlab: ' codeToRun];
hyperlink = sprintf('<a href="%s">%s</a>', codeToRun, text);
return
