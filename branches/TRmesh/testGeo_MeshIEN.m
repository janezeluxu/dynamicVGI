%function [] = testGeo_MeshIEN( )
function [ meshData,vertexData,edgeDatatest,faceData,vertexIEN,edgeIEN,faceIEN] = testGeo_MeshIEN( )
% test geo_Mesh class, output all mesh and IEN cells
%number of points in x

xPoints = 6;
%number of points in y
yPoints = 6;
%x domain coordinates
xDomain = [-1,1];
%y domain coordinates
yDomain = [-1,1];
%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
% 
% %define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [-1,-0.6];
y1 = [-1,-1];
x2 = [1,0.6];
y2 = [1,1];

pV = [2,2];
pH = [2,2];
pT = [2,2];
OrderList = [1,1];

[meshData,vertexData,edgeDatatest,faceData] = testMesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
 
[vertexIEN,edgeIEN,faceIEN] = IENIBC(xPoints,yPoints,meshData,edgeDatatest,faceData);
end

function [meshData,vertexData,edgeData,faceData] = testMesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList)

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
BCtop = 0; 
BCbot = 0; 
BCleft = 0;
BCright = 1; 
global BC
BC = [BCtop BCbot BCleft BCright]; 

[meshData] = mesh.MeshData();
meshData{1}(1,:)
%[edgeData] = mesh.edgeAllMesh(meshData);
[vertexData] = mesh.vertexMesh(meshData);
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
%[IBC1,IBC2,IBC3,IBC4,V,BCval1,BCval2,BCval3,BCval4,VBCval] = mesh.BoundaryCondition(edgeData);
 %%IBC1,IBC2,IBC3,IBC4,V,BCval1,BCval2,BCval3,BCval4,VBCval
end

function [vertexIEN,edgeIEN,faceIEN] = IENIBC(xPoints,yPoints,Mesh,Edge,Face)
global Grid_size
Grid_size = (xPoints-1)*(yPoints-1)*2;
ienbc = IEN(xPoints,yPoints,Mesh,Edge,Face);
[vertexIEN,edgeIEN,faceIEN] = ienbc.Construct_IEN();
end
