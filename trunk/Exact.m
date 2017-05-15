function uex = Exact(x,y,force,a,kappa)
uex = zeros(1,length(x));
if force ==0
    %[uex] =(1-x.^2);
    alpha = a(1)/kappa(1,1);
    b = exp(-alpha);
    [uex] = (exp(-alpha*(1-x))-b)./(1-b);
elseif force ==1
    [uex] = (x.^2-1).*(y.^2-1);
elseif force ==2
    [uex] = (x.^2-1).*(y.^2-1);
elseif force == 3
    [uex] = SteadyHeatAnalytical(x,y);
elseif force ==4
    %uex = (1-x.^4);
    [uex] = sin(x).*sin(y);
elseif force ==5
    for i = 1:length(x)
        uex(i) = PieceWise1(x(i),y(i));
    end
elseif force ==6
    for i = 1:length(x)
        uex(i) = PieceWise2(x(i),y(i));
    end
elseif force ==7
    for i = 1:length(x)
        uex(i) = PieceWise3(x(i),y(i));
    end
elseif force ==8
    for i = 1:length(x)
        uex(i) = PieceWise4(x(i),y(i));
    end
end
end

function u = PieceWise1(x,y)
if x<= -0.6 
    u = (x+1);%*(1-y^2);%*(y+1);
elseif x>0.6 
    u = (1-x);%*(1-y^2);%*(1-y);
else
    u = (-(5/6)*x^2+0.7);%*(1-y^2);
end
%u = (x^2-1);
end

function u = PieceWise2(x,y)
if x<= -0.6 
    u = (x+1);%*(1-y^2);%*(y+1);
elseif x>0.6 
    u = (1-x);%*(1-y^2);%*(1-y);
else
    a = -1/(4*0.6^3);
    b = 0.4-a*0.6^4;
    u = (a*x^4+b);%*(1-y^2);
end
%u = (x^2-1);
end

function u = PieceWise3(x,y)
if x<= -0.6 
    u = (x+1)*(1-y^2);
elseif x>0.6 
    u = (1-x)*(1-y^2);
else
    u = (-(5/6)*x^2+0.7)*(1-y^2);
end
%u = (x^2-1);
end

function u = PieceWise4(x,y)
if (x> -0.6 && x< 0.6 && y> -0.6 && y< 0.6)
    u = (-x^2+0.6^2)*(0.6^2-y^2);
else
    u = 0;
end
end


