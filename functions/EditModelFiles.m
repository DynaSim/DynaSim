function EditModelFiles(files)

[~,eqnfiles]= LocateModelFiles(files);

for file = eqnfiles(:)'
  edit(file{1});
end

end