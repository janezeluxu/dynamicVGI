ucoths = main(1,0,1e-3,-1,0,[0,1],41,[0,1],10,'ro','strong','yes');
errcoths= ErrorEstimate(ucoths,60,41,-1,1e-3);
[u1,errui1] = main(3,50,1e-3,-1,0,[0,1],41,[0,1],10,'r*','strong','no');
errorVGI1 = ErrorEstimate(u1,60,41,-1,1e-3);