function [] = TestFunctions()

%number of points in x
xPoints = 5;
%number of points in y
yPoints = 5;

%x domain coordinates
xDomain = [-1,1];
%y domain coordinates
yDomain = [-1,1];
%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [-1,-0.5];
y1 = [-1,-1];
x2 = [1,0.5];
y2 = [1,1];
OrderList = [2,4];

%define kappa and source term
kappa = [1,0;0,1];

force = 3;

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

for ele = 2: 2
    [errork,errorf] = testVariableStiff(ele,OrderList,IENstruct,vertexData,kappa,force);
    L2K = (sum(sum(errork.^2)))^0.5
    L2F = (sum(sum(errorf.^2)))^0.5
    if  (L2K < 1E-13) && (L2F < 1E-14) 
        disp( 'test element StiffnessMatrix Successed!!!!')
    else
        disp( 'test element StiffnessMatrix Failed!!!!')
    end
    
   [fh]  = testVariableShapeFunc(ele,IENstruct);
end

% p=1;
% n = p+1;
% qPoints = TriGaussPoints(n);
% a = SimplexShapeFunc(qPoints,p);
% [IENall,pAll] = elementIEN(1,IENstruct);
% b = SimplexShapeFunc(qPoints,IENall,pAll);
% 
% for i = 1:size(qPoints,1)
%     IntPoint = qPoints(i,:);
%     %IntPoint = [0,1];
%     [errorSF, errorDiv] = testShapeFunctions(a,IntPoint,p);
%     %[errorSF2, errorDiv2] = testVariablePShapeFunctions(a,b,IntPoint);
%     L2SF = (sum(sum(errorSF.^2)))^0.5;
%     L2Div = (sum(sum(errorDiv.^2)))^0.5;
%     %maxdivError = max(sum(errorDiv2));
%     if  (L2SF < 1E-14) && (L2Div < 1E-15) 
%         disp( 'testUniformShapeFunction Successed!!!!')
%     else
%         disp( 'testUniformShapeFunction Failed!!!!')
%     end
% end
% 
% 
% [errSF, errDiv] = testShapeFunctionsTable(a,qPoints,p);
% L2SF = (sum(sum(errSF.^2)))^0.5;
% L2Div = (sum(sum(errDiv.^2)))^0.5;
% if  (L2SF < 1E-15) && (L2Div  < 1E-15)
%     disp( 'testUniformShapeFuncTable Successed!!!!')
% else
%     disp( 'testUniformShapeFuncTable Failed!!!!')
% end
% 
% 
% for ele = 1:1
%  [errorK,errorF,vertexData] = testElementStiff(ele,p);
%  L2K = (sum(sum(errorK.^2)))^0.5;
%  L2F = ((sum(errorF.^2)))^0.5;
% 
%      if  (L2K < 1E-13) && (L2F < 1E-15) 
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
 %[lexicoOrder] = testLexiOrder(a,p);
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

function [errorK,errorF,vertexData] = testElementStiff(elementNumber,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
xPoints = 6;
yPoints =6;
xDomain = [-1,1];
yDomain = [-1,1];
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
%OrderList = [4];

kappa = [1,0;0,1];
force = 3;

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
n = p+3;
qPoints = TriGaussPoints(n);
a = SimplexShapeFunc(qPoints,p);
[ShapeFuncTable, divShapeFuncTable ] = a.uniformPShapeFunctionTable();
                
[IENall,pAll] = elementIEN(elementNumber,IENstruct);
nInt = max(pAll)+3;
simplex = SimplexStiffMatrix(elementNumber,IENstruct,IENall,nInt,vertexData,kappa,force,ShapeFuncTable,divShapeFuncTable);

[ JInverse, ~,x,y] = simplex.jacobian();
[elementK,elementF] = simplex.eleStiffMatrix();
    
elementF;
x1 = x(1); x2 = x(2); x3 = x(3); y1 = y(1); y2 = y(2); y3 = y(3);
[k,f] = analytical(a,p,x1,x2,x3,y1,y2,y3,kappa,1);

double(f);
f = double(f);
errorK = elementK-double(k);
errorF = elementF-double(f');

k = [9.3492, -3.9108, -5.4383; -3.9108, 2.7192,1.1917; -5.4383, 1.1917, 4.2467];
OnlineK = k/1.3;
OnlineK;

end


function [errork,errorf] = testVariableStiff(ele,p,IENstruct,vertexData,kappa,q)
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

[IENall,pAll] = elementIEN(ele,IENstruct);
n = nIntergerPoints(max(pAll));

if range(pAll) == 0
    %disp('uniform p element');
    ShapeFunc = ShapeFuncTable{pAll(1)};
    divShapeFunc = divSFtable{pAll(1)};
else
    ele
    disp('nonuniform p element');
    qPoints = TriGaussPoints(n);
    simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
    [ShapeFunc, divShapeFunc ] = simplexsf.variablePShapeFunctionTable();
    
end

simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,kappa,q,ShapeFunc,divShapeFunc);

[elementK,elementF] = simplex.eleStiffMatrix();

[~, ~,xCord,yCord] = simplex.jacobian();

q = 1;
[K,F] = variableKFAnalytic(IENall,pAll, xCord,yCord,1,kappa);
k = double(K);

errork = elementK-double(K);
errorf = elementF-double(F');

end

function [k,f] = variableKFAnalytic(IENall,pAll, xCord,yCord,q,kappa)
[ShapeFunc,divShapeFunc] = variableSFAnalytic(IENall,pAll, xCord,yCord);
x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);
syms x y
nsSF = length(IENall);
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

function [fb] = testVariableShapeFunc(ele,IENstruct)
[IENall,pAll] = elementIEN(ele,IENstruct);
xCord = [0,1,0]; yCord = [0,0,1];
[ShapeFunc,divShapeFunc] = variableSFAnalytic(IENall,pAll, xCord,yCord);
ShapeFunc{1}
ShapeFunc{2}
ShapeFunc{3}
ShapeFunc{4}
ShapeFunc{5}
ShapeFunc{6}
ShapeFunc{7}
ShapeFunc{8}

% sumSF = 0;
% sumDiv1 = 0;
% sumDiv2 = 0;
q=0;
nsSF = length(IENall);

syms x y 
x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);
a = x2+y2;

% A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
% lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
% lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
% lam1 = 1-lam2-lam3;

%ShapeFunc = {lam1,lam2,lam3}

% for i = 1:nsSF
%     sumSF = sumSF+ShapeFunc{i};
%     sumDiv1  = sumDiv1+divShapeFunc{i,1};
%     sumDiv2  = sumDiv2+divShapeFunc{i,2};
% end
% sumSF
% sumDiv1
% sumDiv2

Source = 1;

for i = 1:nsSF
    fe = ShapeFunc{i}*Source;
    f(i) = int(int(fe,x,x1,a-y),y,y1,y3);
    for j = 1:nsSF
        Me = ShapeFunc{i}*ShapeFunc{j};
        M(i,j) =  int(int(Me,x,x1,a-y),y,y1,y3);
    end
end
f = double(f)';
M = double(M);
fb = M\f

for i = 1:nsSF
fbE(i) = 1;
end
xA = [0,1,0,1/2,1/2,0,0];
yA = [0,0,1,0,1/2,1/3,2/3];
%fbE = xA+yA;

 fbE'-fb
end

function [Source] = MMS(x,y,q)
            if q ==0
                Source = 1;
            elseif q == 1
                Source = -2*(x.^2+y.^2-2);
            elseif q == 2
                Source = x*y;
            elseif q ==3
                Source = (x^2-1)*(y^2-1);
            elseif q ==4
                Source = 0;
            end
end

function [ShapeFunc,divShapeFunc] = variableSFAnalytic(IENall,pAll, xCord,yCord)
%[lexicoOrder] = testLexiOrder(a,p);

x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);

syms x y 
A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
lam1 = 1-lam2-lam3;

pE1 = pAll(1);
pE2 = pAll(2);
pE3 = pAll(3);
pF = pAll(4);

nssl = length(IENall);
ShapeFunc = cell(nssl,1);
divShapeFunc = cell(nssl,2);

n = nIntergerPoints(max(pAll));
qPoints = TriGaussPoints(n);
simplex = SimplexShapeFunc(qPoints,IENall,pAll);

if pE1>1
    edgeAlpha1 = simplex.getLexico(pE1,'edge',1);
    [ SFedge1, divEdge1 ] = SFAnalytical(edgeAlpha1,xCord,yCord);
   
else
    SFedge1 = [];
    divEdge1 = [];
end
if pE2>1
    edgeAlpha2 = simplex.getLexico(pE2,'edge',2);
    [ SFedge2, divEdge2 ] = SFAnalytical(edgeAlpha2,xCord,yCord);
   
else
    SFedge2 = [];
    divEdge2 = [];
end
if pE3>1
    edgeAlpha3 = simplex.getLexico(pE3,'edge',3);
    [ SFedge3, divEdge3 ] = SFAnalytical(edgeAlpha3,xCord,yCord);
else
    SFedge3 = [];
    divEdge3 = [];
end
if pF>2
    faceAlpha = simplex.getLexico(pF,'face',1);
    [ SFFace, divFace ] = SFAnalytical(faceAlpha,xCord,yCord);
else
    SFFace = [];
    divFace = [];
end


lam = [lam1,lam2,lam3];
dlam = [diff(lam1,x),diff(lam1,y);diff(lam2,x),diff(lam2,y);diff(lam3,x),diff(lam3,y)];

% %corrected vertex mode
for j = 1:3
    sumE1 = 0; sumE2 = 0;sumE3 =0;sumF = 0;
    divSum11 = 0;divSum21=0;divSum31=0;divsumF1=0;
    divSum12 = 0;divSum22=0;divSum32=0;divsumF2=0;
    if pE1>1
        for i = 1:pE1-1
            coeff = edgeAlpha1(i,j)/pE1;
            %SFedge1{1}
            sumE1 = sumE1+coeff*SFedge1{i};
            divSum11 = divSum11+coeff*divEdge1{i,1};
            divSum12 = divSum12+coeff*divEdge1{i,2};
        end
    end
    if pE2>1
        for i = 1:pE2-1
            coeff = edgeAlpha2(i,j)/pE2;
            sumE2 = sumE2+coeff*SFedge2{i};
            divSum21 = divSum21+coeff*divEdge2{i,1};
            divSum22 = divSum22+coeff*divEdge2{i,2};
        end
    end
    if pE3>1
        for i = 1:pE3-1
            coeff = edgeAlpha3(i,j)/pE3;
            sumE3 = sumE3+coeff*SFedge3{i};
            divSum31 = divSum31+coeff*divEdge3{i,1};
            divSum32 = divSum32+coeff*divEdge3{i,2};
        end
    end
    if pF>2
        numberFaceMode = (pF-1)*(pF-2)/2;
        for i = 1:numberFaceMode
            coeff = faceAlpha(i,j)/pF
            sumF = sumF+coeff*SFFace{i};
            divsumF1 = divsumF1 + coeff*divFace{i,1};
            divsumF2 = divsumF2 + coeff*divFace{i,2};
        end
    end
    ShapeFunc{j} = lam(j) - sumE1 - sumE2 - sumE3 - sumF;
    divShapeFunc{j,1} = dlam(j,1)-divSum11-divSum21-divSum31-divsumF1;
    divShapeFunc{j,2} = dlam(j,2)-divSum12-divSum22-divSum32-divsumF2;
end

if nssl>3
ShapeFunc(4:nssl) = [SFedge1;SFedge2;SFedge3;SFFace];
divShapeFunc(4:nssl,:) = [divEdge1;divEdge2;divEdge3;divFace];
end
end


function [nInt] = nIntergerPoints(p)
%define intergral rules
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

function [ShapeFunc,divShapeFunc] = SFAnalytical(lexicoOrder,xCord, yCord)

x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);

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


end

function [ShapeFunc,divShapeFunc] = simplexAnalytical(a,p,lam2, lam3)

[lexicoOrder] = testLexiOrder(a,p);

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