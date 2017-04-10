classdef test_strrep < matlab.unittest.TestCase
  
  methods (Test,  TestTags = {'utility'})
    function testDocExamples(testCase)
      testCase.verifyTrue(strcmp('(pop1_v)*(-av)', ds.strrep('(v)*(-av)','v','pop1_v') ))
      testCase.verifyTrue(strcmp('pop1_v-pop1_v^2+vav', ds.strrep('v-v^2+vav','v','pop1_v') ))
      testCase.verifyTrue(strcmp('pop1_v-pop1_v-pop1_v', ds.strrep('v-v-v','v','pop1_v') ))
      testCase.verifyTrue(strcmp('pop1_v-pop1_v-pop1_v^2', ds.strrep('v-v-v^2','v','pop1_v') ))
      testCase.verifyTrue(strcmp('(pop1_v-pop1_v-pop1_v^2)', ds.strrep('(v-v-v^2)','v','pop1_v') ))
      testCase.verifyTrue(strcmp('E-pop1_V(n-1)+1', ds.strrep('E-pop1_V+1','pop1_V','pop1_V(n-1)') ))
      testCase.verifyTrue(strcmp('v=1; u(n,test)=u(n,test)+d', ds.strrep('v=1; u=u+d','u','u(n,test)') ))
    end
    
    function testSubString(testCase)
      testCase.verifyTrue(strcmp('new.old', ds.strrep('old.old','old','new') ))
      testCase.verifyTrue(strcmp('new.old.old', ds.strrep('old.old.old','old','new') ))
    end
  end
  
end