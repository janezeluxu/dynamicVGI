function [ meshData,vertexData,edgeData,faceData,IBC,vertexIEN,edgeIEN,faceIEN] = testGeo_MeshIEN( )
% test geo_Mesh class, output all mesh and IEN cells
%number of points in x
xPoints = 10;
%number of points in y
yPoints = 10;
%x domain coordinates
xDomain = [-1,1];
%y domain coordinates
yDomain = [-1,1];
%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [-1,-0.5];
y1 = [-1,-0.5];
x2 = [1,0.5];
y2 = [1,0.5];
OrderList = [2,4];

[meshData,vertexData,edgeData,faceData,IBC] = testMesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
 
[vertexIEN,edgeIEN,faceIEN] = IENIBC(xPoints,yPoints,meshData,edgeData,faceData);
end

function [meshData,vertexData,edgeData,faceData,IBC] = testMesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList)

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[IBC] = mesh.BoundaryConditionDOF(edgeData);
 
end

function [vertexIEN,edgeIEN,faceIEN] = IENIBC(xPoints,yPoints,Mesh,Edge,Face)
ienbc = IEN(xPoints,yPoints,Mesh,Edge,Face);
[vertexIEN,edgeIEN,faceIEN] = ienbc.Construct_IEN();
end
