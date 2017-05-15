function [] = testMain()
xPoints = 5;
yPoints = 5;
xDomain = [-1,1];
yDomain = [-1,1];
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
OrderList = [1];

kappa = [1,0;0,1];
force = 3;
tic

[ u ] = testGlobalKF( xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList,kappa,force );
max(u)
toc

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
[edgeData] = mesh.edgeMesh(meshData);
[x,y,value] = DiagnalPlot(u,xPoints,edgeData,vertexData);
%  
%  for i = 1:length(x)
%     s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;
%  end
%  figure (2)
%  plot(s,value,'bo')
%  hold on
%  [uex] = SteadyHeatAnalytical(x,y);
%  max(uex);
%  error = sum(uex-value');
%  
%  %plot analytical solution
%  x = sort(x);
%  y = sort(y,'descend');
%  for i = 1:length(x)
%     s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;
%  end
%  [uex] = SteadyHeatAnalytical(x,y);
%   plot(s,uex,'k-')
  
end

function [ u ] = testGlobalKF( xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList,kappa,force )
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

[kGlobal, fGlobal] = globalKF(OrderList,IENstruct,IBC,TotalDOF,...
                                vertexData,kappa,force,xPoints,yPoints);
[K,F] = KFfullMatrx(OrderList,IENstruct,IBC,TotalDOF,...
                                vertexData,kappa,force,xPoints,yPoints);

u = K\F;
end

function [K,F] = KFfullMatrx(p,IENstruct,IENB,TotalDOF,...
                                vertexData,kappa,q,M,N)
% Assemble of element level matrix   

%create uniform p shape function tables, p is a list
ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    n = nIntergerPoints(p(j));
    qPoints = TriGaussPoints(n);
    simplexsf = SimplexShapeFunc(qPoints,p(j));
    [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();    
    ShapeFuncTable{p(j)} =  ShapeFunc;
    divSFtable{p(j)} = divShapeFunc;
end

Grid_size = (M-1)*(N-1)*2;
K = zeros(TotalDOF,TotalDOF);
F = zeros(TotalDOF,1);
for ele = 1:Grid_size

%for ele = 1:1    
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions    
    [IENall,pAll] = elementIEN(ele,IENstruct);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll));
    
    if range(pAll) == 0
        %disp('uniform p element');
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
    else
        %disp('nonuniform p element');
        qPoints = TriGaussPoints(n);
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, divShapeFunc ] = simplexsf.variablePShapeFunctionTable();
    
    end
    
    %find element Shape Functions, insert stiffness matrix and force vectors
    %into global matrixs(spars)
   
    simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,kappa,q,ShapeFunc,divShapeFunc);

    [elementK,elementF] = simplex.eleStiffMatrix();

    %construct K,F directly
    ele
    for i = 1:nssl
        rowTemp = IENall(i);
        F(rowTemp) = elementF(i) + F(rowTemp);
        for j = 1:nssl           
            columnTemp = IENall(j);
            K(rowTemp,columnTemp) = K(rowTemp,columnTemp) + elementK(i,j);
            
        end
    end
    K77 = K(7,7)
    %add boundary conditions
    for i = 1: length (IENB)
        F(IENB(i)) = 0;
        K(IENB(i),:) = 0;
        K(:,IENB(i)) = 0;
        K(IENB(i),IENB(i)) = 1;
    end
    
end
                            
end


function [nInt] = nIntergerPoints(p)
nInt = p+3;
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

