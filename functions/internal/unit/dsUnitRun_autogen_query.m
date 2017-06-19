function result = dsUnitRun_autogen_query(query)
%% Run queried autogen tests
if ~exist('query','var')
  dsUnitInputQuery
else
  evalin('base', ['query=''' query ''';'])
end

fprintf('Running autogen tests with query: ''%s''\n\n', query)
result = runtests('dsUnitTest_autogen_query', 'UseParallel', true);
end