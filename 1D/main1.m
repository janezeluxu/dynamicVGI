function [u,upnorm] = main1(itau,maxIter,kappa,a,BCval,nPoints,domain,FigNum,lineType,weakOrStrong,Projection)
global h;
global q;

nElements = nPoints-1;
h = zeros(1,nElements);
x = zeros(1,nPoints);

for i = 1:nElements
    h(i) = (domain(2)-domain(1))/(nPoints-1);
end

for i = 1:nPoints
    x(i) = domain(1)+ sum(h(1:i-1));
end

c1 = ones(nPoints,1)*0.0;
c2 = 1;
upnorm = 0;
dt = 3/7;
if itau == 3
    for iteration = 1:maxIter
        [u] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,B,c2);
        [c1new,U] = solvec1LS( u,x,nPoints,kappa,a,c2);
        
        c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
                
        upnorm = upHmorn(u,1,nPoints,a,c1,c2);
    end
        
end

figure (FigNum)
plot(x,u,lineType)
hold on
%plot(x,U,'ro')
%hold on
%plot(x1,ux1,'co')
%xlim([0 0.5])
r = a(1)/kappa(1,1);
b = exp(-r);
[ue] = (exp(r/2*x)+1)./(exp(r/2)+1).*(exp(r/2*x)-1)./(exp(r/2)-1); 

Pe = norm(a)*h(1)/(2*norm(kappa));
fprintf('Mesh Peclet Number: %f\n',Pe);

end