function [u,upnorm,iteration] = main(itau,maxIter,kappa,a,f,BCval,nPoints,domain,FigNum,lineType,Projection,c2)
global h;
global q;
B = 0;
p=1;
q = 1.0;
nElements = nPoints-1;
h = zeros(1,nElements);
x = zeros(1,nPoints);

if q == 1
    h1 = (domain(2)-domain(1))/(nPoints-1);
else
    h1 = (1-q)/(1-q^(nElements));
end

for i = 1:nElements
    %h(i) = (domain(2)-domain(1))/(nPoints-1);
    h(i) = h1*q^(i-1);
end

for i = 1:nPoints
    x(i) = domain(1)+ sum(h(1:i-1));
end
x1 = [x(1),x(3:41)];

c1 = ones(nPoints,1)*0.0;

l1 = ones(nPoints,1)*0;
m1 = ones(nPoints,1)*0;

[u] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,f,c2);   
%[c10,~,l1new,m1new] = solvec1LS( u0,x,nPoints,kappa,a,c2);

if strcmp(Projection, 'LS') == 1
    [l1new,m1new,~] = solvec1LS( u,x,nPoints,kappa,a,c2);
elseif strcmp(Projection, 'MULTI') == 1
    [l1new,m1new,~] = solvec1Project( u,x,nPoints,kappa,a,c2,BCval);
else
    fprintf('Should specify projection');
end
        
for node = 3:nPoints-2
    l1new(node) = (l1new(node-1)+l1new(node)+l1new(node+1));
    m1new(node) = (m1new(node-1)+m1new(node)+m1new(node+1));
end
c10 = zeros(nPoints,1);
for node = 2:nPoints-1
    c10(node) = (m1new(node)/l1new(node))/(a^2);
    if c10(node)<0|| l1new(node) ==0
        c10(node) = 0;
    end
end
         
upnorm_u0 = upHmorn(u,1,nPoints,a,c10,c2);
%u = 0;
%ux1 = solveUstrong(nPoints-1,x1,BCval,c1,itau,kappa,a,B,c2); 
%U = 0;
uTol = 1e-4;
upnorm = 0;
dt = 3/7;
iteration = 0;
if itau == 3
    for iteration = 1:maxIter
        %iteration
        [u] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,f,c2); 

        if strcmp(Projection, 'LS') == 1
            [l1new,m1new,c1new] = solvec1LS( u,x,nPoints,kappa,a,c2);
        elseif strcmp(Projection, 'MULTI') == 1
            [l1new,m1new,c1new] = solvec1Project( u,x,nPoints,kappa,a,c2,BCval);
        else
            fprintf('Should specify projection');
        end

        upnorm_u = upHmorn(u,1,nPoints,a,c1,c2);
        c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
        %upwind
        %c1(1)=c1(2);
        c1(nPoints) = c1(nPoints-1);
        
        %c1(1) = 2*c1(2)-c1(3);
        %c1(nPoints) = 2*c1(nPoints-1)-c1(nPoints-2);
        [unew] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,f,c2);
        
        normu = norm(u-unew)/norm(u);
        upnorm_unew = upHmorn(unew,1,nPoints,a,c1,c2);
        up = norm(upnorm_u-upnorm_unew)/norm(upnorm_u0);
        if  (normu<uTol )%||up <uTol) 
            break
        end
        
    end
        
end
c1;
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

%figure (20)
%iterations = 1:1:iteration;
%plot(iterations,upnorm,'b*-')
%xlim([0 0.1])
%ylim([0 1])
%loglog(x,ue,'b-')
Pe = norm(a)*h(1)/(2*norm(kappa));
%fprintf('Mesh Peclet Number: %f\n',Pe);
title(horzcat('Pe_{cell} = ',num2str(Pe)),'FontSize',30);
%title('Element Number = 80, Coarse  Grid');
%legend('Coth Classical Tau','Compressible Tau','Weak BC local VGI')
%legend('8 Element','16 Element','32 Element','64 Element','128 Element')
hold on
end

