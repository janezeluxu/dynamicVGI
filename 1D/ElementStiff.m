function [k,f] = ElementStiff(a,kappa,B,Stab,J,InvJ,l,element)
p=1;
n = nIntergerPoints(p,0);
[x,w] = GaussQuad(n,1);
[ShapeFunc,DivSF] = ShapeFuncDiv(p,x);
k = zeros(p+1,p+1);
%f = zeros(p+1,1);
for i = 1:length(w)
    Nax = DivSF(:,i)*InvJ;
    Nbx = Nax';
    Na = ShapeFunc(:,i);
    Nb = Na';
    k = k+Nax*kappa*Nbx*J*w(i)-Nax*a*Nb*J*w(i)+Na*B*Nb*J*w(i)+a*Stab*Nax*Nbx*J*w(i);
    %f = f+
end
f = [l(element)*(0.5*J-Stab);l(element)*(0.5*J+Stab)];
end

function [ShapeFunc,DivSF] = ShapeFuncDiv(p,x)
L = length(x);
ShapeFunc = zeros(p+1,L);
DivSF = zeros(p+1,L);
for v = 0:p
    coeff = nchoosek(p,v);
    ShapeFunc(v+1,:) = coeff*x.^v.*(1-x).^(p-v);
    DivSF(v+1,:) = coeff*v*x.^(v-1).*(1-x).^(p-v)-coeff*x.^(v).*(p-v).*(1-x).^(p-v-1);
end
end