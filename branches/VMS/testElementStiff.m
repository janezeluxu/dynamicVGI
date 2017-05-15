function [] = testElementStiff()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
xPoints = 2;
yPoints =2;
xDomain = [-1,1];
yDomain = [-1,1];
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
OrderList = [1];

kappa = [4.7,0;0,5.1];
force = 1;

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[IBC] = mesh.BoundaryConditionDOF(edgeData);

ien = IEN(xPoints,yPoints,meshData,edgeData,faceData);
[vertexIEN,edgeIEN,faceIEN] = ien.Construct_IEN();

TotalDOF = max([vertexIEN{:,2},edgeIEN{:,:},faceIEN{:,2}]);
IENstruct = struct('vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN);

p = OrderList(1);
n = p+1;
qPoints = TriGaussPoints(n);
a = SimplexShapeFunc(qPoints,p);
[ShapeFuncTable, divShapeFuncTable ] = a.uniformPShapeFunctionTable();
                
[IENall,pAll] = elementIEN(1,IENstruct);
nInt = max(pAll)+1;
simplex = SimplexStiffMatrix(1,IENstruct,IENall,nInt,vertexData,kappa,force,ShapeFuncTable,divShapeFuncTable);

[ JInverse, detJ] = simplex.jacobian()
[elementK,elementF] = simplex.eleStiffMatrix()

%p = 1 test
k = [9.3492, -3.9108, -5.4383; -3.9108, 2.7192,1.1917; -5.4383, 1.1917, 4.2467];
OnlineK = k/1.3;
% error = elementK-OnlineK;
% if error<1E-5
%     disp('Pass StiffnessMatrix test')
% else
%     disp('StiffnessMatrix test Failed')
% end

%p = 3 test
k = analytical(p);
f = analytical(f);
end


function [IENall,pAll] = elementIEN(ele,IENstruct)
[~, IENvertex,~] = IENstruct(ele,:).vertexIEN;

[edge1, edge2,edge3] = IENstruct(ele,:).edgeIEN;

if edge1(1) >1 %there is edge modes, pEdge>1
    IENedge1 = edge1(2:end);
else
    IENedge1 = [];
end

if edge2(1) >1 %there is edge modes, pEdge>1
    IENedge2 = edge2(2:end);
else
    IENedge2 = [];
end

if edge3(1) >1 %there is edge modes, pEdge>1
    IENedge3 = edge3(2:end);
else
    IENedge3 = [];
end

[~, face,~] = IENstruct(ele,:).faceIEN;

if face(1) > 2 %there is face modes, pFace>2
    IENface = face(2:end);
else
    IENface = [];
end

IENall = [IENvertex,IENedge1,IENedge2,IENedge3,IENface];
pAll = [edge1(1),edge2(1),edge3(1),face(1)];
end

function [k,f] = analytical(p)


end


