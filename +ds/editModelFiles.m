function editModelFiles(files)

[~,eqnfiles]= ds.locateModelFiles(files);

for file = eqnfiles(:)'
  edit(file{1});
end

end
