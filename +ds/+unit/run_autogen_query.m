%% Run queried autogen tests
query = 'applyModifications';

fprintf('Running autogen tests with query: ''%s''\n\n', query)
result = runtests('ds.unit.test_autogen_query');