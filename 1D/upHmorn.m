function [unorm] = upHmorn(u,nInt,Grid_Size,a,c1,c2)
global h
[xi,weight] = GaussQuad(nInt, 1);
[uxi, uxixi] = uxiuxixi(u,Grid_Size);
nElement = Grid_Size-1;
unorm = 0;

%uxixi = zeros(nElement,1);
for ele = 1:nElement
    J = h(ele);
    c1ele = 0.5*(c1(ele+1)+c1(ele));
    tau = c1ele*(h(ele)^c2);
    for np = 1:length(weight)
        %upx = (J*uxixi(ele)-J^2*uxi(ele))/(J^3);
        upx = (uxi(ele)*J);
        unorm = unorm+a^2*tau^2*upx^2*J*weight(np);
    end
end

unorm = unorm^0.5;
end



function [uxi, uxixi] = uxiuxixi(u,Grid_Size)
nElement = Grid_Size-1;
uxi = zeros(nElement,1);
uxiR = zeros(Grid_Size,1);
uxixi = zeros(nElement,1);
for ele = 1:nElement
    uxi(ele) = -u(ele)+u(ele+1);
end

uxiR(1)=uxi(1);
uxiR(nElement+1) = uxi(nElement);
for ele = 2:nElement
    uxiR(ele)=(uxi(ele-1)+uxi(ele))/2;
    uxixi(ele) = -uxiR(ele)+uxiR(ele+1);  
end
uxixi(1) = -uxiR(1)+uxiR(2); 
end