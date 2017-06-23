function result = run_autogen_query(query)
%% Run queried autogen tests
if ~exist('query','var')
  dsUnit.inputQuery
else
  evalin('base', ['query=''' query ''';'])
end

fprintf('Running autogen tests with query: ''%s''\n\n', query)
result = runtests('dsUnit.test_autogen_query', 'UseParallel', true);
end