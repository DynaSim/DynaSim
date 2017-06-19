function editModelFiles(files)

[~,eqnfiles]= dsLocateModelFiles(files);

for file = eqnfiles(:)'
  edit(file{1});
end

end
