function[TotalDOF,vertexData,edgeData,faceData,IBC,BCval,uHBCE,uHBCV,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,varargin)

global nsd;
global direction;
global BCType;
if length(varargin) ==1
    OrderList = varargin{1};
    mesh = geo_Mesh_Square(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
elseif length(varargin) ==3
    PV = varargin{1};
    PH = varargin{2};
    PT = varargin{3};
    mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,PV,PH,PT);
end

%mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh(meshData);
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[IBC1,IBC2,IBC3,IBC4,V,BCval1,BCval2,BCval3,BCval4,VBCval] = mesh.BoundaryCondition(edgeData,vertexData);
%IBC1,IBC2,IBC3,IBC4,V,BCval1,BCval2,BCval3,BCval4,VBCval

uHBCV(1,:) = stronguH(vertexData,V,[IBC1,IBC2,IBC3,IBC4]);
uHBCV(2,:) = [uHBCV(1,1)-xPoints,uHBCV(1,2)-xPoints,V(3:end)];

uHBCE1 = stronguH(vertexData,[IBC1],[IBC2,IBC3,IBC4,V]);
VAL1 = uHBCE1-xPoints-1;
uHBCE2 = stronguH(vertexData,[IBC2],[IBC1,IBC3,IBC4,V]);
VAL2 = uHBCE2+xPoints+1;
uHBCE3 = stronguH(vertexData,[IBC3],[IBC1,IBC2,IBC4,V]);
VAL3 = uHBCE3-xPoints-1;
uHBCE4 = stronguH(vertexData,[IBC4],[IBC1,IBC2,IBC3,V]);
VAL4 = uHBCE4+xPoints+1;
uHBC = [];
uHBCval = [];
if (nsd ==2 && strcmp(BCType{1},'Strong')==1)   
    uHBC = [uHBC,uHBCE1];
    uHBCval = [uHBCval,VAL1];
end
if (nsd ==2 && strcmp(BCType{2},'Strong')==1)   
    uHBC = [uHBC,uHBCE2];
    uHBCval = [uHBCval,VAL2];
end
if (nsd ==2 && strcmp(BCType{3},'Strong')==1)   
        uHBC = [uHBC,uHBCE3];
    uHBCval = [uHBCval,VAL3];
end
if (nsd ==2 && strcmp(BCType{4},'Strong')==1)   
        uHBC = [uHBC,uHBCE4];
    uHBCval = [uHBCval,VAL4];
end
uHBCE(1,:) = uHBC;
uHBCE(2,:) = uHBCval;

% uHBCE1 = setdiff(uHBCE1,uHBCV(1,:));
% uHBCE2 = setdiff(uHBCE2,uHBCV(1,:));
% uHBCE3 = setdiff(uHBCE3,uHBCV(1,:));
% uHBCE4 = setdiff(uHBCE4,uHBCV(1,:));

if ( nsd ==1 && strcmp(direction,'X')==1)
    %[IBC1,IBC2,IBC3,IBC4,~,~,BCval3,BCval4] = mesh.BoundaryCondition1D(edgeData);
    IBC3 = [IBC3,V(1),V(3)];
    IBC4 = [IBC4,V(2),V(4)];
    BCval3 = [BCval3,VBCval(1),VBCval(3)];
    BCval4 = [BCval4,VBCval(2),VBCval(4)];
    IBC = {IBC1;IBC2;IBC3;IBC4};
    BCval = {BCval1,BCval2,BCval3,BCval4};
    disp ('1D in X');
elseif (nsd ==1 && strcmp(direction,'Y')==1)
    %[IBC1,IBC2,IBC3,IBC4,BCval1,BCval2,~,~] = mesh.BoundaryCondition1D(edgeData);
    IBC1 = [IBC1,V(1),V(2)];
    IBC2 = [IBC2,V(3),V(4)];
    BCval1 = [BCval1,VBCval(1),VBCval(2)];
    BCval2 = [BCval2,VBCval(3),VBCval(4)];
    IBC = {IBC1;IBC2;IBC3;IBC4};
    BCval = {BCval1,BCval2,BCval3,BCval4};
    disp ('1D in Y');
elseif (nsd ==2 && strcmp(BCType{1},'Strong')==1||strcmp(BCType{1},'flux')==1)
    IBC1 = [IBC1,V(1)];
    IBC2 = [IBC2,V(3),V(4)];
    IBC4 = [IBC4,V(2)];
    
    BCval1 = [BCval1,VBCval(1)];
    BCval2 = [BCval2,VBCval(3),VBCval(4)];
    BCval4 = [BCval4,VBCval(2)];
    
    IBC = {IBC1;IBC2;IBC3;IBC4};
    BCval = {BCval1,BCval2,BCval3,BCval4};
elseif (nsd ==2 && (strcmp(BCType{1},'Outlet')==1 ||strcmp(BCType{1},'Inlet')==1))
    IBC1 = [IBC1,V(1),V(2)];
    IBC2 = [IBC2,V(3),V(4)];
    IBC3 = [IBC3,V(1),V(3)];
    IBC4 = [IBC4,V(2),V(4)];
    BCval1 = [BCval1,VBCval(1),VBCval(2)];
    BCval2 = [BCval2,VBCval(3),VBCval(4)];
    BCval3 = [BCval3,VBCval(1),VBCval(3)];
    BCval4 = [BCval4,VBCval(2),VBCval(4)];
    IBC = {IBC1;IBC2;IBC3;IBC4};
    BCval = {BCval1,BCval2,BCval3,BCval4};
   %disp ('2D!!!');
end

ien = IEN(xPoints,yPoints,meshData,edgeData,faceData);
[vertexIEN,edgeIEN,faceIEN] = ien.Construct_IEN();

TotalDOF = max([vertexIEN{:,2},edgeIEN{:,:},faceIEN{:,2}]);
IENstruct = struct('vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN);

% Build periodic
iper = zeros(TotalDOF,1); 
for i = 1:TotalDOF
    iper(i) = i;
    if(strcmp(BCType{1},'PeriodicY')==1)  % Only set up for 1D right now
        for j = 1:length(IBC1)
            iper(IBC1(j)) = IBC2(j);
        end
        
    elseif(strcmp(BCType{3},'PeriodicX')==1)  % Only set up for 1D right now
        for j = 1:length(IBC3)
            iper(IBC3(j)) = IBC4(j);
        end
    end
end
end
 
function uHBC = stronguH(vertexData,BC,V)
SurroundNode = [];
for i = 1:length(BC)
    SurroundNode = [SurroundNode,vertexData{BC(i),4}];
end
a = unique(SurroundNode);
uHBC =setxor(a,BC);
uHBC =setdiff(uHBC,V);
end