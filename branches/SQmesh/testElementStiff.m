function [] = testElementStiff()
%This program solve the 2D diffution equation for variable or unifrom order
%simplex elements, the mesh is structured triangles.
maxIter = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nsd; 
global stabflag;
global itau; 
global kappa; 
global direction
global force
global BCType;
global a;
global TotalDOF;
global BC; 
global tol; 
global Grid_size;
global node_Size;
global gamma;
global h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = [1,0;0,1]*1e-3;
force = 0;
a = [-1 -1];
nin = [1;0];
nout = [-1;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 2; 
xPoints = 3;   % number of points in x
yPoints = 3;   % number of points in x
xDomain = [0,1];   % limit of domain in x
yDomain = [0,1];   % limit of domain in x

%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [0];
y1 = [0];
x2 = [1];
y2 = [1];
OrderList = [1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-7; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCbot = 0;
BCtop = 0;  
BCleft = 0;
BCright = 1; 
BC = [BCbot BCtop BCleft BCright]; 
direction = 'X';
BCType = {'Strong','Strong','Strong','Strong'};
gamma = 1;
global CbI;
CbI = 4;
if(nsd==1)
    h = range(xDomain)/(xPoints-1);    
else
    h = sqrt(0.5*range(xDomain)/(xPoints-1) * range(yDomain)/(yPoints-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stabflag = 1;   % Add stabilization (1: VMS, 2: ...)
itau = 2;  
Grid_size = (xPoints-1)*(yPoints-1);
node_Size = xPoints*yPoints;%+(xPoints-1)*(yPoints-1);
nInt = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the mesh data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TotalDOF,vertexData,~,~,IBC,BCval,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList); 
 
[IENall,pAll] = elementIEN(1,IENstruct);
nInt = max(pAll)+1;
%[ShapeFuncTable,divSFtable] = ShapeTable(1,1);
[ShapeFuncTable,divSFtable] = ShapeTable(0,pAll(1));

ShapeFunc = ShapeFuncTable{pAll(1)};
divShapeFunc = divSFtable{pAll(1)};
                
c1Ele = 0;
%simplex = SimplexStiffMatrix(1,IENstruct,IENall,nInt,vertexData,kappa,force,ShapeFuncTable,divShapeFuncTable);
simplex = SquareStiffMatrix(1,IENstruct,IENall,2,vertexData,ShapeFunc,divShapeFunc,c1Ele);
[ JInverse, detJ] = simplex.jacobian();
[elementK,elementF] = simplex.eleStiffMatrix();
%BCvalEle = 1;
%[InletK,InletF] = simplex.WeakBCInletStiffMatrix([1;0],BCvalEle);
%[OutletK,OutletF] = simplex.WeakBCOutlettiffMatrix([-1;0],BCvalEle);
%InletK
% error = elementK-OnlineK;
% if error<1E-5
%     disp('Pass StiffnessMatrix test')
% else
%     disp('StiffnessMatrix test Failed')
% end

p = OrderList(1);
[k,f] = analyticalInlet(JInverse,nin,0);
double(k)
%f = analytical(f);
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

function [elementK,elementF] = analyticalInlet(JInverse,nin,BCval)
global kappa; 
global a;
global h;
global CbI;
syms lam1 
NaGlobal =  [lam1;0;1-lam1];
NbGlobal = [lam1;0;1-lam1];

gradNaGlobal = [-1,-1;1,0;0,1]*JInverse;
gamma = 1;
elementF = (-gamma*(kappa*gradNaGlobal')'*nin*h)*BCval;
elementF = elementF-a*nin*NaGlobal*BCval*h;
elementF = elementF+(CbI*norm(kappa)/h)*(NaGlobal)*h*BCval;

elementK =  int( -NaGlobal*((kappa*gradNaGlobal')'*h*nin)'-gamma*(kappa*gradNaGlobal')'*nin*NbGlobal'*h...
 +(CbI*norm(kappa)/h)*(NaGlobal*NbGlobal')*h,lam1,0,1);



end


