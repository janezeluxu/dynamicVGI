function [ c1] = solvec1Weight( u,x,nPoints,kappa,a,c2)

global h 
c1 = zeros(nPoints,1);
U = zeros(nPoints,1);
for node = 2:nPoints-1
%for node = 7:7
    Pe = norm(a)*h(node-1)/(2*norm(kappa));
    H = h(node-1)+h(node);
    
    dfaia = h(node-1);
    dfaib = h(node);
    dfai2 = H;

    %fprintf('Mesh Peclet Number: %f\n',Pe);
    uh = [Pe*u(node-1);u(node+1)];
    CordMatrix = [1,x(node-1);1,x(node+1)]; 
    
    %uh = [u(node-1);u(node);u(node+1)];
    %CordMatrix = [1,x(node-1);1,x(node);1,x(node+1)]; 
    
    A = CordMatrix\uh;
    a1 = A(1);
    a2 = A(2);
            
    dUx = a2;
    uHa = a1+a2*0.5*(x(node-1)+x(node));
    uHb = a1+a2*0.5*(x(node+1)+x(node));
    uH2 = a1+a2*0.5*(x(node)+x(node));
    
    uha = 0.5*(u(node-1)+u(node));
    uhb = 0.5*(u(node)+u(node+1));
    duha = -u(node-1)+u(node);
    duhb = -u(node)+u(node+1);
    
    A = dfaia*duha*h(node-1)^(c2-1)+dfaib*duhb*h(node)^(c2-1);
    B = H^(c2)*dUx*dfai2;
    
    C = -a*dfai2*uH2; 
    %C = -a*(uHa*dfaia+uHb*dfaib);
    D = kappa*dUx*dfai2;
    E = a*(-dfaia*uha)+a*(-dfaib*uhb);
    %F = kappa*dfaia*duha+kappa*dfaib*duhb;
    F = kappa*((dfaia/h(node-1))*duha+(dfaib/h(node))*duhb);
    
    l1 = A-B;
    m1 = C+D-(E+F);
    
    c1(node) = (m1/l1)/(a^2);
    
    if c1(node)<0
        c1(node) = 0;
    end
end
end