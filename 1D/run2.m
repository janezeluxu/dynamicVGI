% clear all
nPoints = 41;
l = ones(nPoints-1,1)*0;
bcValue = [0,1];
c2 = 2;
    %ucoth = main(1,0,1e-3,-1,l,[0,0],41,[0,1],21,'cs-','strong','LS','pathTubec1');
   [u3,upnorm3,~] = main(0,20,1e-1,1,l,bcValue,41,[0,1],3,'k-','LS',c2);
   
   [u3,upnorm3,~] = main(2,20,1e-1,1,l,bcValue,41,[0,1],3,'bs-','LS',c2);
   
   [u3,upnorm3,iteration] = main(3,20,1e-1,1,l,bcValue,41,[0,1],3,'r*-','MULTI',c2);
   
   iteration
   %title('Strong BC Converged uh no c1 upwind')
   %legend('Analytical Solution','Static Tau','Dynamic Tau');%'no Stronly imposed uH upwind c1new','Stronly imposed uH upwind c1new')
   
   
%    [u3,upnorm3,iteration] = main(3,20,1e-3,-1,l,bcValue,41,[0,1],31,'bo-','strong','LS','pathLinec1');
%    iteration
%    [u3,upnorm3,iteration] = main(3,20,1e-3,-1,l,bcValue,41,[0,1],31,'ro-','strong','MULTI','pathLinec1');
%    iteration
%    title('Strong BC Converged uh upwind c1new')
%    legend('no Stronly imposed uH','Stronly imposed uH');
%    
%    [u3,upnorm3,iteration] = main(3,20,1e-3,-1,l,bcValue,41,[0,1],32,'bs-','strong','LS','pathLinec1');
%    iteration
%    [u3,upnorm3,iteration] = main(3,20,1e-3,-1,l,bcValue,41,[0,1],32,'rs-','strong','MULTI','pathLinec1');
%    iteration
%    title('Strong BC Converged uh')
%    legend('no Stronly imposed uH','Stronly imposed uH');
   
   %'no Stronly imposed uH upwind c1new','Stronly imposed uH upwind c1new')
   
   
   %[u3,upnorm3,iteration] = main(3,20,1e-3,-1,l,bcValue,41,[0,1],31,'b*-','weak','LS','pathLinec1');
   %iteration
   %title('weak BC Converged uh')
   %figure (21)
   %iterations = 1:1:20;
   %plot(iterations,upnorm1,'r*-')
   %hold on
   %plot(iterations,upnorm2,'b*-')
   %hold on
   %plot(iterations,upnorm3,'ro-')
   %hold on
   %plot(iterations,upnorm4,'bo-')

l = ones(nPoints-1,1)*1;
bcValue = [0,0];
for Iter = 0:5
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'b*-','weak','LS','pathLinec1');
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'k*-','weak','LS','pathTubec1');
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'r*-','weak','LS','pathLineL1M1');
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'c*-','weak','LS','pathTubeL1M1');

   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'b*-','strong','MULTI','pathLinec1');
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'k*-','strong','MULTI','pathTubec1');
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'r*-','strong','MULTI','pathLineL1M1');
   %[u3,upnorm3] = main(3,Iter,1e-3,-1,l,bcValue,41,[0,1],Iter+1,'c*-','strong','MULTI','pathTubeL1M1');
    
   %legend('pathLinec1','pathTubec1','pathLineL1M1','pathTubeL1M1')
    %errorVGI1(Iter+1) = ErrorEstimate(u1,60,41,-1,1e-3);
    %errorVGI2(Iter+1) = ErrorEstimate(u2,60,41,-1,1e-3);
    %errorVGI3(Iter+1) = ErrorEstimate(u3,60,6,-1,1e-3);
    %errorVGI4(Iter+1) = ErrorEstimate(u4,60,6,-1,1e-3);
end
% 
% ucoths = main(1,0,1e-3,-1,[0,1],41,[0,1],10,'ro','strong','yes');
% errcoths= ErrorEstimate(ucoths,60,41,-1,1e-3);
% ucothw = main(1,0,1e-3,-1,0,[0,1],41,[0,1],10,'ro','weak','yes');
% errcothw= ErrorEstimate(ucothw,60,41,-1,1e-3);
% figure(40)
% plot(0:1:Iter,errorVGI1,'r*-')
%hold on
%plot(0:1:Iter,errorVGI2,'b*-')
% %hold on
% plot(0:1:Iter,errorVGI3,'b*-')
% hold on
% plot(0:1:Iter,errorVGI4,'bo-')
 %hold on
 %plot(0:1:Iter,ones(Iter+1,1)*errcoths,'r-')
% %hold on
% hold on
% plot(0:1:Iter,ones(Iter+1,1)*errcothw,'b-')
% hold on