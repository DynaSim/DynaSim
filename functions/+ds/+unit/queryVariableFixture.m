classdef queryVariableFixture < matlab.unittest.fixtures.Fixture
  
  properties (SetAccess=private)
    query
  end
  
  methods
    function fixture = queryVariableFixture(query)
      fixture.query = query;
    end
    
    function setup(fixture)
    end
  end
  
end