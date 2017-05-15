function [u] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,l,c2)
global h;

%source
nElement = nPoints-1;
%l = zeros(nElement,1)*0.2;

IEN(:,1) = 1:nPoints-1;
IEN(:,2) = 2:nPoints;
K = zeros(nPoints,nPoints);
F = zeros(nPoints,1);

for element = 1:nElement
    J = x(element+1)-x(element);
    InvJ = 1/J;
    %element
    c1ele = 0.5*(c1(element+1)+c1(element));
    tau = tauFunc(element,c1ele,itau,kappa,a,InvJ,c2);
    %cc1 = tau/h(element)
    Stab = a*tau;
    k = [InvJ*(kappa+a*Stab)-0.5*a, -InvJ*(kappa+a*Stab)+0.5*a;-InvJ*(kappa+a*Stab)-0.5*a,InvJ*(kappa+a*Stab)+0.5*a];
    f = [l(element)*(0.5*J-Stab);l(element)*(0.5*J+Stab)];
    
    %[k,f] = ElementStiff(p,a,kappa,B,Stab,J,InvJ,l,element);
    
    %assemble    
    K(IEN(element,1),IEN(element,1)) = k(1,1)+K(IEN(element,1),IEN(element,1));
    K(IEN(element,1),IEN(element,2)) = k(1,2)+K(IEN(element,1),IEN(element,2));
    K(IEN(element,2),IEN(element,1)) = k(2,1)+K(IEN(element,2),IEN(element,1));
    K(IEN(element,2),IEN(element,2)) = k(2,2)+K(IEN(element,2),IEN(element,2));
    
    F(IEN(element,1)) = f(1)+F(IEN(element,1));
    F(IEN(element,2)) = f(2)+ F(IEN(element,2));
end

%% apply Strong boundary condition
K(1,:) = 0;
K(1,1) = 1;
K(nPoints,:) = 0;
K(nPoints,nPoints) = 1;

F(1) = BCval(1);
F(nPoints) = BCval(2);

u = K\F;


end

function [tau] = tauFunc(element,c1element,itau,kappa,a,InvJ,c2)
global h;
if(itau == 1)
    alpha = norm(a)*h(element)/(2*norm(kappa));
    zi = coth(alpha)-1/alpha;
    tau = h(element)*zi/(2*norm(a));
elseif(itau == 2)     
    tau1 = h(element)/2*a(1);
    tau2 = h(element)^2/(4*kappa(1,1));
    tau = 1/sqrt(tau1^-2 + 9*tau2^-2);
elseif(itau == 3)                    
    tau = c1element*(h(element)^c2); 
elseif(itau == 4)
    tau = (2*InvJ*norm(a)+2*InvJ^2*kappa)^-1;
elseif(itau == 0)
    tau = 0;
end

end