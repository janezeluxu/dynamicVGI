% weak and no upwind
clear all
ucoth = main(1,0,1e-3,-1,[0,1],6,[0,1],2,'bo','weak','yes');
[u,erru,iteration] = main(3,500,1e-3,-1,[0,1],6,[0,1],2,'r*','weak','yes');
[error] = ErrorEstimate(u,60,6,-1,1e-3);
fprintf('Test weak boundary condition with no upwind weighting, meshSize is 5 \n');
fprintf('Convergence error: %d\n',erru);
fprintf('Used number of iterations: %d\n',iteration);
fprintf('local VGI Solution error to analytical: %d\n',error);

clear all
ucoth = main(1,0,1e-3,-1,[0,1],11,[0,1],2,'bo','weak','no');
[u,erru,iteration] = main(3,500,1e-3,-1,[0,1],11,[0,1],2,'r*','weak','yes');
[error] = ErrorEstimate(u,100,11,-1,1e-3);
fprintf('Test weak boundary condition with no upwind weighting, meshSize is 10 \n');
fprintf('Convergence error: %d\n',erru);
fprintf('Used number of iterations: %d\n',iteration);
fprintf('local VGI Solution error to analytical: %d\n',error);

clear all
ucoth = main(1,0,1e-3,-0.5,[0,1],21,[0,1],3,'b*','weak','no');
[u,erru,iteration] = main(3,500,1e-3,-0.5,[0,1],21,[0,1],3,'r*','weak','yes');
[errorVGI] = ErrorEstimate(u,60,21,-1,1e-3);
[errorcoth] = ErrorEstimate(ucoth,60,21,-1,1e-3);
fprintf('Test weak boundary condition with no upwind weighting, meshSize is 20 \n');
fprintf('Convergence error: %d\n',erru);
fprintf('Used number of iterations: %d\n',iteration);
fprintf('local VGI Solution error to analytical: %d\n',errorVGI);
fprintf('coth Solution error to analytical: %d\n',errorcoth);

clear all
ucoth = main(2,0,1e-3,-1,[0,1],41,[0,1],4,'bo','weak','no');
[u,erru,iteration] = main(3,500,1e-3,-1,[0,1],41,[0,1],4,'r*','weak','yes');
[errorVGI] = ErrorEstimate(u,60,41,-1,1e-3);
[errorcoth] = ErrorEstimate(ucoth,60,41,-1,1e-3);
fprintf('Test weak boundary condition with no upwind weighting, meshSize is 40 \n');
fprintf('Convergence error: %d\n',erru);
fprintf('Used number of iterations: %d\n',iteration);
fprintf('local VGI Solution error to analytical: %d\n',errorVGI);
fprintf('coth Solution error to analytical: %d\n',errorcoth);
% 
clear all
ucoth = main(1,0,1e-3,-1,[0,1],81,[0,1],5,'b*','weak','no');
%ucompress = main(4,0,1e-3,-1,[0,1],81,[0,1],5,'k*','weak','no');
[u,erru,iteration] = main(3,50,1e-3,-1,[0,1],81,[0,1],5,'r*','weak','yes');
[errorVGI] = ErrorEstimate(u,60,81,-1,1e-3);
fprintf('Test weak boundary condition with no upwind weighting, meshSize is 80 \n');
fprintf('Convergence error: %d\n',erru);
fprintf('Used number of iterations: %d\n',iteration);
fprintf('local VGI Solution error to analytical: %d\n',errorVGI);
% % 
clear all
ucoth = main(1,0,1e-3,-1,[0,1],161,[0,1],6,'bo','weak','no');
[u,erru,iteration] = main(3,500,1e-3,-1,[0,1],161,[0,1],6,'r*','weak','yes');
[error] = ErrorEstimate(u,60,161,-1,1e-3);
fprintf('Test weak boundary condition with no upwind weighting, meshSize is 160 \n');
fprintf('Convergence error: %d\n',erru);
fprintf('Used number of iterations: %d\n',iteration);
fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],321,[0,1],7,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],321,[0,1],7,'r*','weak','yes');
% [errorVGI] = ErrorEstimate(u,60,321,-1,1e-3);
% [errorcoth] = ErrorEstimate(ucoth,60,321,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 320 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',errorVGI);
% fprintf('coth Solution error to analytical: %d\n',errorcoth);
% 
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],641,[0,1],2,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],641,[0,1],2,'r*','weak','yes');
% [error] = ErrorEstimate(u,60,641,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 640 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],1281,[0,1],2,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],1281,[0,1],2,'r*','weak','yes');
% [error] = ErrorEstimate(u,60,1281,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 1280 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],2561,[0,1],2,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],2561,[0,1],2,'r*','weak','yes');
% [error] = ErrorEstimate(u,60,2561,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 2560 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],5121,[0,1],2,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],5121,[0,1],2,'r*','weak','yes');
% [error] = ErrorEstimate(u,60,5121,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 5120 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],10241,[0,1],2,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],10241,[0,1],2,'r*','weak','yes');
% [error] = ErrorEstimate(u,60,10241,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 10240 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],20481,[0,1],2,'bo','weak','no');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],20481,[0,1],2,'r*','weak','yes');
% [error] = ErrorEstimate(u,60,20481,-1,1e-3);
% fprintf('Test weak boundary condition with no upwind weighting, meshSize is 20480 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);

%% strong and upwind
% clear all
% ucoth = main(2,0,1e-3,-1,[0,1],6,[0,1],2,'bo','strong','no');
% [u,erru,iteration] = main(3,10,1e-3,-1,[0,1],6,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,6,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 10 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],11,[0,1],3,'bo','strong','yes');
% [u,erru,iteration] = main(3,500,1e-3,-1,[0,1],11,[0,1],3,'r*','strong','yes');
% [error] = ErrorEstimate(u,60,11,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 10 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],21,[0,1],4,'bo','strong','yes');
% [u,erru,iteration] = main(3,500,1e-3,-1,[0,1],21,[0,1],4,'r*','strong','yes');
% [error] = ErrorEstimate(u,60,21,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 20 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],41,[0,1],5,'bo','strong','yes');
% [u,erru,iteration] = main(3,500,1e-3,-1,[0,1],41,[0,1],5,'r*','strong','yes');
% [error] = ErrorEstimate(u,60,41,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 40 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],81,[0,1],6,'bo','strong','yes');
% [u,erru,iteration] = main(3,500,1e-3,-1,[0,1],81,[0,1],6,'c*','strong','yes');
% [error] = ErrorEstimate(u,60,81,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 80 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],161,[0,1],7,'bo','strong','yes');
% [u,erru,iteration] = main(3,500,1e-3,-1,[0,1],161,[0,1],7,'r*','strong','yes');
% [error] = ErrorEstimate(u,60,161,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 160 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],321,[0,1],8,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],321,[0,1],8,'r*','strong','yes');
% [error] = ErrorEstimate(u,60,321,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 320 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],641,[0,1],2,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],641,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,641,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 640 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],1281,[0,1],2,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],1281,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,1281,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 1280 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],2561,[0,1],2,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],2561,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,2561,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 2560 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],5121,[0,1],2,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],5121,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,5121,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 5120 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],10241,[0,1],2,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],10241,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,10241,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 10240 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);
% 
% clear all
% ucoth = main(1,0,1e-3,-1,[0,1],20481,[0,1],2,'bo','strong','yes');
% [u,erru,iteration] = main(3,50,1e-3,-1,[0,1],20481,[0,1],2,'r*','strong','no');
% [error] = ErrorEstimate(u,60,20481,-1,1e-3);
% fprintf('Test strong boundary condition with upwind weighting, meshSize is 20480 \n');
% fprintf('Convergence error: %d\n',erru);
% fprintf('Used number of iterations: %d\n',iteration);
% fprintf('local VGI Solution error to analytical: %d\n',error);