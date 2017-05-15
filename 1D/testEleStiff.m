function [] = testEleStiff()
a = -1;
tau = 0.099;
InvJ = 5;
kappa =1e-3; 
N = 1;
[K,Kweak] = ele_kf(a,tau,InvJ,kappa,N);
[k,kweak] = ele_k(a,tau,InvJ,kappa);

errot  = K - k;
erroweak  = Kweak - kweak;
L2 = (sum(sum(errot.^2)))^0.5;
L2weak = (sum(sum(erroweak.^2)))^0.5;
if  (L2 < 1E-15)&&(L2weak<1E-15)
   disp( 'test error estimate Successed!!!!')
else
   disp( 'test error estimate Failed!!!!')
end

end

function [K,Kweak] = ele_kf(a,tau,InvJ,kappa,N)
Nxi = [-1;1];
[xi, w] = GaussQuad(N, 1);
nssl = length(Nxi);
K = zeros(nssl,nssl);
Kweak = zeros(nssl,nssl);
F = zeros(nssl,1);
for i = 1:N
    Na = [1-xi(i);xi(i)];
    K = K+InvJ*kappa*Nxi*Nxi'*w(i)+a*Na*Nxi'*w(i)+tau*a^2*InvJ*Nxi*Nxi';
    Kweak = Kweak+InvJ*kappa*Nxi*Nxi'*w(i)-a*Nxi*Na'*w(i)+tau*a^2*InvJ*Nxi*Nxi';
end
end

function [k,kweak] = ele_k(a,tau,InvJ,kappa)
Stab = a*tau;
element = 5;
l = zeros(element,1);
J = 1/InvJ;
k = [InvJ*(kappa+a*Stab)-0.5*a, -InvJ*(kappa+a*Stab)+0.5*a;-InvJ*(kappa+a*Stab)-0.5*a,InvJ*(kappa+a*Stab)+0.5*a];
f = [l(element)*(0.5*J-Stab);l(element)*(0.5*J+Stab)];
kweak = (a*Stab+kappa)*InvJ*[1,-1;-1,1]-a*[-0.5,-0.5;0.5,0.5];
end