classdef dsUnitTest_strrep < matlab.unittest.TestCase
  
  methods (Test,  TestTags = {'utility'})
    function testDocExamples(testCase)
      testCase.verifyTrue(strcmp('(pop1_v)*(-av)', dsStrrep('(v)*(-av)','v','pop1_v') ))
      testCase.verifyTrue(strcmp('pop1_v-pop1_v^2+vav', dsStrrep('v-v^2+vav','v','pop1_v') ))
      testCase.verifyTrue(strcmp('pop1_v-pop1_v-pop1_v', dsStrrep('v-v-v','v','pop1_v') ))
      testCase.verifyTrue(strcmp('pop1_v-pop1_v-pop1_v^2', dsStrrep('v-v-v^2','v','pop1_v') ))
      testCase.verifyTrue(strcmp('(pop1_v-pop1_v-pop1_v^2)', dsStrrep('(v-v-v^2)','v','pop1_v') ))
      testCase.verifyTrue(strcmp('E-pop1_V(n-1)+1', dsStrrep('E-pop1_V+1','pop1_V','pop1_V(n-1)') ))
      testCase.verifyTrue(strcmp('v=1; u(n,test)=u(n,test)+d', dsStrrep('v=1; u=u+d','u','u(n,test)') ))
    end
    
    function testSubString(testCase)
      testCase.verifyTrue(strcmp('new.old', dsStrrep('old.old','old','new') ))
      testCase.verifyTrue(strcmp('new.old.old', dsStrrep('old.old.old','old','new') ))
    end
  end
  
end
