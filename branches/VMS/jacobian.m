function [ JInverse, detJ,gij,x,y] = jacobian(simplex)

global nsd; 

%calculate jacobian inverse and det J in lambda space, can put
%this function into private access
elementID = simplex.EleID;
IENstruct = simplex.IENstructure;
vertexData = simplex.VertexData;

[~, vIDs,~] = IENstruct(elementID,:).vertexIEN;
x1 = vertexData{vIDs(1),2}(1);
y1 = vertexData{vIDs(1),2}(2);
x2 = vertexData{vIDs(2),2}(1);
y2 = vertexData{vIDs(2),2}(2);
x3 = vertexData{vIDs(3),2}(1);
y3 = vertexData{vIDs(3),2}(2);

x = [x1,x2,x3];
y = [y1,y2,y3];
A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));

JInverse = (0.5/A)*[(y3-y1),(x1-x3);(y1-y2),(x2-x1);];

gradx = [(-x1+x2), (-x1+x3);(-y1+y2), (-y1+y3)];
detJ = 0.5*(gradx(1,1)*gradx(2,2)-gradx(1,2)*gradx(2,1));

gij = zeros(nsd,nsd); 
dxi1(1) = 0.5/A*(y2-y3);
dxi1(2) = 0.5/A*(x3-x2);
dxi2(1) = 0.5/A*(y3-y1);
dxi2(2) = 0.5/A*(x1-x3); 

for i = 1:nsd
    for j = 1:nsd
        gij(i,j) = dxi1(i)*dxi1(j) + dxi2(i)*dxi2(j);
    end
end
    
end