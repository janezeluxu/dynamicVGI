function [l1,m1,c1] = solvec1Project( u,x,nPoints,kappa,a,c2,bcValue)

global h 
c1 = zeros(nPoints,1);
l1 = zeros(nPoints,1);
m1 = zeros(nPoints,1);
U = zeros(nPoints,1);
for node = 2:nPoints-1
    H = h(node-1)+h(node);
    
    uh = [u(node-1);u(node);u(node+1)];
    CordMatrix = [1,x(node-1);1,x(node);1,x(node+1)];
    
    [a1, a2] = LeastSquare( CordMatrix, uh,node,nPoints,bcValue);
    
    dUx = a2;
    uHa = a1+a2*0.5*(x(node-1)+x(node));
    uHb = a1+a2*0.5*(x(node+1)+x(node));
    uH2 = a1+a2*x(node);
    U(node) = uH2;
    
    uha = 0.5*(u(node-1)+u(node));
    uhb = 0.5*(u(node)+u(node+1));
    duha = -u(node-1)+u(node);
    duhb = -u(node)+u(node+1);
    
    A = duha*h(node-1)^(c2)+duhb*h(node)^(c2);
    B = H^(c2)*dUx*H;
    
    C = -a*H*uH2;
    %C = -a*(uHa*h(node-1)+uHb*h(node));
    D = kappa*dUx*H;
    E = a*(-h(node-1)*uha)+a*(-h(node)*uhb);
    F = kappa*(duha+duhb);
    
    l1(node) = (A-B);
    m1(node) = (C+D-(E+F));
end    
    


  for node = 2:nPoints-1
     c1(node) = (m1(node)/l1(node))/(a^2);
     if c1(node)<0 || abs(l1(node)) <=1e-20
         c1(node) = 0;
     end
 end
   
end
function [a1, a2] = LeastSquare( CordMatrix, uh,node,nPoints,bcValue)
if node ==2
    value = bcValue(1);
    A = CordMatrix\uh;
    a2 = A(2);
    a1 = value-a2*CordMatrix(1,2);
elseif node == nPoints-1
    value = bcValue(2);
    A = CordMatrix\uh;
    a2 = A(2);
    a1 = value-a2*CordMatrix(3,2);
else
    A = CordMatrix\uh;
    a1 = A(1);
    a2 = A(2);
end
end