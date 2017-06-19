classdef dsUnitTest_strrep2 < matlab.unittest.TestCase
  
  methods (Test,  TestTags = {'utility'})
    function testDocExamples(testCase)
      testCase.verifyEqual('(pop1_v)*(-av)', dsStrrep2('(v)*(-av)','v','pop1_v') )
      testCase.verifyEqual('pop1_v-pop1_v^2+vav', dsStrrep2('v-v^2+vav','v','pop1_v') )
      testCase.verifyEqual('pop1_v-pop1_v-pop1_v', dsStrrep2('v-v-v','v','pop1_v') )
      testCase.verifyEqual('pop1_v-pop1_v-pop1_v^2', dsStrrep2('v-v-v^2','v','pop1_v') )
      testCase.verifyEqual('(pop1_v-pop1_v-pop1_v^2)', dsStrrep2('(v-v-v^2)','v','pop1_v') )
      testCase.verifyEqual('E-pop1_V(n-1)+1', dsStrrep2('E-pop1_V+1','pop1_V','pop1_V(n-1)') )
      testCase.verifyEqual('v=1; u(n,test)=u(n,test)+d', dsStrrep2('v=1; u=u+d','u','u(n,test)') )
    end
    
    function testSubString(testCase)
      testCase.verifyEqual('new.new', dsStrrep2('old.old','old','new') )
      testCase.verifyEqual('new.new.new', dsStrrep2('old.old.old','old','new') )
    end
  end
  
end
