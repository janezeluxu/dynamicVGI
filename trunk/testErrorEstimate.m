function [ output_args ] = testErrorEstimate( input_args )
%number of points in x
xPoints = 6;
%number of points in y
yPoints = 6;

%x domain coordinates
xDomain = [-1,1];
%y domain coordinates
yDomain = [-1,1];
%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
OrderList = [4];

%define kappa and source term
kappa = [1,0;0,1];

force = 1;

mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
[vertexData] = mesh.vertexMesh();

ien = IEN(xPoints,yPoints,meshData,edgeData,faceData);

[vertexIEN,edgeIEN,faceIEN] = ien.Construct_IEN();
TotalDOF = max([vertexIEN{:,2},edgeIEN{:,:},faceIEN{:,2}]);
IENstruct = struct('meshData',meshData,'vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN);
u = zeros(TotalDOF,1);

%loop  only one element
GridSize = (xPoints-1)*(yPoints-1)*2;
%GridSize = 1;
nInt = 10;
[Error] = errorEstimate(u,xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList,force,GridSize,nInt);
ErrorA = analyticError(force,IENstruct,vertexData,GridSize,OrderList);
errot  = Error - double(ErrorA);
L2 = (sum(sum(errot.^2)))^0.5;
if  (L2 < 1E-13)
   disp( 'test error estimate Successed!!!!')
else
   disp( 'test error estimate Failed!!!!')
end
end


function IntEx = analyticError(force,IENstruct,vertexData,GridSize)

IntEx = 0;
IntA = 0;
for i = 1:GridSize
    [~, vIDs,~] = IENstruct(i,:).vertexIEN;
    x1 = vertexData{vIDs(1),2}(1);
    y1 = vertexData{vIDs(1),2}(2);
    x2 = vertexData{vIDs(2),2}(1);
    y2 = vertexData{vIDs(2),2}(2);
    x3 = vertexData{vIDs(3),2}(1);
    y3 = vertexData{vIDs(3),2}(2);
    
%     %x1 = 0;x2 = 1;x3 = 0;y1 = 0;y2 = 0;y3 = 1;
%     xCord = [x1,x2,x3];
%     yCord = [y1,y2,y3];
%     
%     gradx = [(-x1+x2), (-x1+x3);(-y1+y2), (-y1+y3)];
%     detJ = 0.5*(gradx(1,1)*gradx(2,2)-gradx(1,2)*gradx(2,1));
%             
%     %OrderList
%     n = OrderList(1)+5;        
%     qPoints = TriGaussPoints(n);
%         
%     lam = [1-qPoints(:,1)-qPoints(:,2),qPoints(:,1:2)];
%     x = xCord*lam';
%     y = yCord*lam';
%     uA = Exact(x,y,force);
%     IntA = (uA.^2*qPoints(:,3)*detJ)+IntA;
    
    syms x y
    a = x2+y2;
    uex = Exact(x,y,force);
    IntEx = int(int(uex^2,x,x1,a-y),y,y1,y3)+IntEx;
end

end

function uex = Exact(x,y,force)
if force ==1
    [uex] = (x.^2-1).*(y.^2-1);
elseif force ==2
    [uex] = sin(x).*sin(y);
elseif force == 3
    [uex] = SteadyHeatAnalytical(x,y);
elseif force ==4
    uex = zeros(1,length(x));
    for i = 1:length(x)
        uex(i) = PieceWise(x(i),y(i));
    end
end
end


