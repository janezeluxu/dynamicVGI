%function [] = testGeo_MeshIEN( )
function [ meshData,vertexData,edgeDatatest,faceData,IBC,vertexIEN,edgeIEN,faceIEN] = testGeo_MeshIEN( )
% test geo_Mesh class, output all mesh and IEN cells
%number of points in x
xPoints = 11;
%number of points in y
yPoints = 11;
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
OrderList = [1,2];

[meshData,vertexData,edgeDatatest,faceData,IBC] = testMesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,pV,pH,pT);
 
[vertexIEN,edgeIEN,faceIEN] = IENIBC(xPoints,yPoints,meshData,edgeDatatest,faceData);
end

function [meshData,vertexData,edgeData,faceData,IBC] = testMesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,pV,pH,pT)

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,pV,pH,pT);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
%[edgeData] = mesh.edgeAllMesh(meshData);
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[IBC] = mesh.BoundaryConditionDOF(edgeData);
 
end

function [vertexIEN,edgeIEN,faceIEN] = IENIBC(xPoints,yPoints,Mesh,Edge,Face)
ienbc = IEN(xPoints,yPoints,Mesh,Edge,Face);
[vertexIEN,edgeIEN,faceIEN] = ienbc.Construct_IEN();
end
