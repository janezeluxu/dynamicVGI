function[TotalDOF,vertexData,IBC,BCval,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList)

global nsd;

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[IBC,BCval] = mesh.BoundaryConditionDOF(edgeData);

ien = IEN(xPoints,yPoints,meshData,edgeData,faceData);
[vertexIEN,edgeIEN,faceIEN] = ien.Construct_IEN();

TotalDOF = max([vertexIEN{:,2},edgeIEN{:,:},faceIEN{:,2}]);
IENstruct = struct('vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN);

% Build periodic
iper = zeros(TotalDOF,1); 
for i = 1:TotalDOF
    if(nsd==1 && i<=(xPoints-1) && i ~=1)  % Only set up for 1D right now
        iper(i) = i+(yPoints-1)*xPoints;
    else
        iper(i) = i;
    end
end

end