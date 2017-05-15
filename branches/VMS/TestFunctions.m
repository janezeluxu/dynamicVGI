function [] = TestFunctions()
p=4;
n = p+1;
qPoints = TriGaussPoints(n);
a = SimplexShapeFunc(qPoints,p);

xPoints = 10;
yPoints = 10;
xDomain = [-1,1];
yDomain = [-1,1];
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
OrderList = [4];

kappa = [1,0;0,1];
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

[IENall,pAll] = elementIEN(1,IENstruct);
b = SimplexShapeFunc(qPoints,IENall,pAll);

for i = 1:size(qPoints,1)
    IntPoint = qPoints(i,:);
    %IntPoint = [0,1];
    [errorSF, errorDiv] = testShapeFunctions(a,IntPoint,p);
    [errorSF2, errorDiv2] = testVariablePShapeFunctions(a,b,IntPoint);
    maxdivError = max(sum(errorDiv2));
    if  (maxdivError < 1E-13) && (max(errorSF2) < 1E-13) 
        disp( 'testUniformShapeFunction Successed!!!!')
    else
        disp( 'testUniformShapeFunction Failed!!!!')
    end
end


[errSF, errDiv] = testShapeFunctionsTable(a,qPoints,p);
if  (errSF < 1E-13) && (max(errDiv) < 1E-13)
    disp( 'testUniformShapeFuncTable Successed!!!!')
else
    disp( 'testUniformShapeFuncTable Failed!!!!')
end


% for ele = 1:1
%  [errorK,errorF] = testElementStiff(ele,p);
%      if  (max(sum(errorK)) < 1E-15) && (max(errorF) < 1E-15) 
%         disp( 'test element StiffnessMatrix Successed!!!!')
%     else
%         disp( 'test element StiffnessMatrix Failed!!!!')
%     end
% end
 end

function [errSF, errDiv] = testShapeFunctionsTable(a,qPoints,p)

[ShapeFuncTable, divShapeFuncTable ] = a.uniformPShapeFunctionTable();
nQuadPoints = size(qPoints,1);
ShapeFuncAna = zeros((p+1)*(p+2)/2,nQuadPoints);
divShapeFuncAna = zeros((p+1)*(p+2)/2,2,nQuadPoints);
for i = 1: nQuadPoints
    IPoint = qPoints(i,:);
    syms lam2 lam3
    [SFAnalytical,divAnalytical] = simplexAnalytical(a,p,lam2, lam3); 
    lambda2 = IPoint(1);
    lambda3 = IPoint(2);

    SFAnalytical = subs(subs(SFAnalytical,lam2,lambda2),lam3,lambda3);
    divAnalytical = subs(subs(divAnalytical,lam2,lambda2),lam3,lambda3);
 
    ShapeFuncAna(:,i) =  SFAnalytical;
    divShapeFuncAna(:,:,i) = divAnalytical;
end

errSF = max(sum(ShapeFuncTable - ShapeFuncAna));
errDiv = max(sum(divShapeFuncTable - divShapeFuncAna));

end

function [errorSF, errorDiv] = testShapeFunctions(a,IntPoint,p)
[ShapeFunc, divShapeFunc] = a.uniformPShapeFunction(IntPoint);
 syms lam2 lam3
 [SFAnalytical,divAnalytical] = simplexAnalytical(a,p,lam2, lam3);
 
 lambda2 = IntPoint(1);
 lambda3 = IntPoint(2);
 

 %lambda2 = 0;
 %lambda3 = 0;
 
 SFAnalytical = subs(subs(SFAnalytical,lam2,lambda2),lam3,lambda3);
 divAnalytical = subs(subs(divAnalytical,lam2,lambda2),lam3,lambda3);
 
 errorSF = ShapeFunc - double(SFAnalytical);
 errorDiv = divShapeFunc - double(divAnalytical);
end

function [errorSF, errorDiv] = testVariablePShapeFunctions(a,b,IntPoint)
[ ShapeFunc, divShapeFunc ] = b.variablePShapeFunction(IntPoint);
[ShapeFuncUniform, divShapeFuncUniform] = a.uniformPShapeFunction(IntPoint);
errorSF = ShapeFunc - ShapeFuncUniform;
errorDiv = divShapeFunc - divShapeFuncUniform;
end

function [errorK,errorF] = testElementStiff(elementNumber,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
xPoints = 10;
yPoints =10;
xDomain = [-1,1];
yDomain = [-1,1];
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
%OrderList = [4];

kappa = [1,0;0,1];
force = 1;

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,p);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[IBC] = mesh.BoundaryConditionDOF(edgeData);

ien = IEN(xPoints,yPoints,meshData,edgeData,faceData);
[vertexIEN,edgeIEN,faceIEN] = ien.Construct_IEN();

TotalDOF = max([vertexIEN{:,2},edgeIEN{:,:},faceIEN{:,2}]);
IENstruct = struct('vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN);

%p = OrderList(1);
n = p+1;
qPoints = TriGaussPoints(n);
a = SimplexShapeFunc(qPoints,p);
[ShapeFuncTable, divShapeFuncTable ] = a.uniformPShapeFunctionTable();
                
[IENall,pAll] = elementIEN(elementNumber,IENstruct);
nInt = max(pAll)+1;
simplex = SimplexStiffMatrix(elementNumber,IENstruct,IENall,nInt,vertexData,kappa,force,ShapeFuncTable,divShapeFuncTable);

[ ~, ~,x,y] = simplex.jacobian();
[elementK,elementF] = simplex.eleStiffMatrix();
x1 = x(1); x2 = x(2); x3 = x(3); y1 = y(1); y2 = y(2); y3 = y(3);
[k,f] = analytical(a,p,x1,x2,x3,y1,y2,y3,kappa,force);
errorK = elementK-double(k);
errorF = elementF-double(f');

%k = [9.3492, -3.9108, -5.4383; -3.9108, 2.7192,1.1917; -5.4383, 1.1917, 4.2467];
%OnlineK = k/1.3;

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


function [k,f] = analytical(a,p,x1,x2,x3,y1,y2,y3,kappa,q)

[lexicoOrder] = testLexiOrder(a,p);
syms x y 
A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
lam1 = 1-lam2-lam3;

nsSF = size(lexicoOrder,1);
ShapeFunc = cell(nsSF,1);
divShapeFunc = cell(nsSF,2);

for i = 1: nsSF
    a1 = lexicoOrder(i,1);
    a2 = lexicoOrder(i,2);
    a3 = lexicoOrder(i,3);
    
    b1 = lam1.^a1;
    b2 = lam2.^a2;
    b3 = lam3.^a3;
    
    p = a1+a2+a3;
    coeff = factorial(p)/(factorial(a1)*factorial(a2)*factorial(a3));
    ShapeFunc{i} = coeff*b1*b2*b3;
    divShapeFunc{i,1} = diff(ShapeFunc(i),x);
    divShapeFunc{i,2} = diff(ShapeFunc(i),y);
end
a = x2+y2;
for i = 1:nsSF
    fe = ShapeFunc{i}*q;
    f(i) = int(int(fe,x,x1,a-y),y,y1,y3);
    for j = 1:nsSF
        ke = [divShapeFunc{i,1},divShapeFunc{i,2}]*kappa*[divShapeFunc{j,1};divShapeFunc{j,2}];
        k(i,j) =  int(int(ke,x,x1,a-y),y,y1,y3);
    end
end
end

function [ShapeFunc,divShapeFunc] = simplexAnalytical(a,p,lam2, lam3)
%lam2 = IntPoint(1);
%lam3 = IntPoint(2);

[lexicoOrder] = testLexiOrder(a,p);

%syms lam2 lam3
lam1 = 1-lam2-lam3;

nsSF = size(lexicoOrder,1);
ShapeFunc = cell(nsSF,1);
divShapeFunc = cell(nsSF,2);

for i = 1: nsSF
    a1 = lexicoOrder(i,1);
    a2 = lexicoOrder(i,2);
    a3 = lexicoOrder(i,3);
    
    b1 = lam1.^a1;
    b2 = lam2.^a2;
    b3 = lam3.^a3;
    
    p = a1+a2+a3;
    coeff = factorial(p)/(factorial(a1)*factorial(a2)*factorial(a3));
    ShapeFunc{i} = coeff*b1*b2*b3;
    divShapeFunc{i,1} = diff(ShapeFunc(i),lam2);
    divShapeFunc{i,2} = diff(ShapeFunc(i),lam3);
end


end

function lexicoOrder = testLexiOrder(a,p)
vertexAlpha = a.getLexico(p,'vertex',0);

edgeAlpha1 = a.getLexico(p,'edge',1);
edgeAlpha2 = a.getLexico(p,'edge',2);
edgeAlpha3 = a.getLexico(p,'edge',3);
faceAlpha = a.getLexico(p,'face',0);
lexicoOrder = [vertexAlpha;edgeAlpha1;edgeAlpha2;edgeAlpha3;faceAlpha];

end