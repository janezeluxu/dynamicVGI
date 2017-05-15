function [ue] = Exact(x,a,kappa)
r = a(1)/kappa(1,1);
b = exp(-r);
[ue] = (exp(r/2*x)+1)./(exp(r/2)+1).*(exp(r/2*x)-1)./(exp(r/2)-1);
end