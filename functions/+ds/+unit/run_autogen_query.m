function result = run_autogen_query(query)
%% Run queried autogen tests
if ~exist('query','var')
  ds.unit.inputQuery
else
  evalin('base', ['query=''' query ''';'])
end

fprintf('Running autogen tests with query: ''%s''\n\n', query)
result = runtests('ds.unit.test_autogen_query');
end