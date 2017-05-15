function [x,y,value] = Reconstruct(u,M,N,xDomain,yDomain,x1,y1,x2,y2,OrderList,xLine,yLine,nPoints)
% %2D figure only vertex
% nVertex = size(vertexData,1);
% uVertex = u(1:nVertex);
% 
% %dx = 1/M;
% %dy = 1/N;
% 
% [X,Y] = meshgrid(1:N, 1:M);
% 
% Z = zeros(M,N);
% for i = 1:N
%     Z(:,i) = uVertex(1+(i-1)*M:i*M);
%     
% end
% %surf(X,Y,Z)

%line plot

xL = xDomain(1);     
xR = xDomain(2);    
yL = yDomain(1);    
yU = yDomain(2); 

xs = xLine(1);
xt = xLine(2);
ys = yLine(1);
yt = yLine(2);

if xs == xt %a vertical line
    dy = abs(yt-ys)/nPoints;
    y = ys:dy:yt;
    %x = [xs];
    x = ones(1,nPoints+1)*xs;
    
    elseif ys ==yt %a horizontal line
        dx = abs(xt-xs)/nPoints;
        x = xs:dx:xt;
        y = ones(1,nPoints+1)*ys;
    else %not a vertical or horizontal line
        k = (ys-yt)/(xs-xt);
        b = ys-k*xs;
        dx = abs(xt-xs)/nPoints;
        
        x = xs:dx:xt;
        y = k*x+b;
        
end
        

mesh = geo_Mesh(M,N,xDomain,yDomain,x1,y1,x2,y2,OrderList);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh();
[edgeData] = mesh.edgeMesh(meshData);
[faceData] = mesh.faceMesh(meshData,edgeData);
%[IBC] = mesh.BoundaryConditionDOF(edgeData);

ien = IEN(M,N,meshData,edgeData,faceData);
[vertexIEN,edgeIEN,faceIEN] = ien.Construct_IEN();

%TotalDOF = max([vertexIEN{:,2},edgeIEN{:,:},faceIEN{:,2}]);
IENstruct = struct('meshData',meshData,'vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN);

%find element #based on x,y
value = zeros(1,length(x));
for i = 1:length(x)
    %s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;

    ele = faceNum(x(i),y(i),xL,yL,xR,yU,M,N,IENstruct,vertexData,edgeData);
    % point is a vertex
    if length(ele)==2
        value(i) = u(ele(1));
        
    elseif length(ele)==2
        %point is on an edge
        
    else
        [~, vIDs,~] = IENstruct(ele,:).vertexIEN;
        x1 = vertexData{vIDs(1),2}(1);
        y1 = vertexData{vIDs(1),2}(2);
        x2 = vertexData{vIDs(2),2}(1);
        y2 = vertexData{vIDs(2),2}(2);
        x3 = vertexData{vIDs(3),2}(1);
        y3 = vertexData{vIDs(3),2}(2);
        A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
        
        lam2 = ((x3*y1-x1*y3)+(y3-y1)*x(i)+(x1-x3)*y(i))/(2*A);
        lam3 = ((x1*y2-x2*y1)+(y1-y2)*x(i)+(x2-x1)*y(i))/(2*A);
        quadraturePoint = [lam2,lam3];
        
        [IENall,pAll] = elementIEN(ele,IENstruct);
        if range(pAll) == 0
            for j = 1:length(IENall)
                
                %[ ShapeFunc, ~ ] = uniformPShapeFunction(pAll(1),quadraturePoint);
                %qPoints = TriGaussPoints(n);
                simplexsf = SimplexShapeFunc(quadraturePoint,pAll(1));
                [ShapeFunc, ~ ] = simplexsf.uniformPShapeFunction(quadraturePoint);     
                value(i) =  value(i)+u(IENall(j))*ShapeFunc(j);
                
            end
            
        else
            for j = 1:length(IENall)
                %[ ShapeFunc, ~ ] = variablePShapeFunction(ele,IENstruct,quadraturePoint);
                
                simplexsf = SimplexShapeFunc(quadraturePoint,IENall,pAll);
                [ShapeFunc, ~ ] = simplexsf.variablePShapeFunction(quadraturePoint);
        
                value(i) = value(i) + u(IENall(j))*ShapeFunc(j);
            end
        end
    end
end

end

%from x y cordinate find element number
function elementID = faceNum(xCord,yCord,xL,yL,xR,yU,M,N,IENstruct,vertexData,edgeData)
dx = (xR-xL)/(M-1);
dy = (yU-yL)/(N-1);

xIndL = floor((xCord-xL)/dx)+1;
xIndU = ceil((xCord-xL)/dx)+1;
yIndL = floor((yCord-yL)/dy)+1;
yIndU = ceil((yCord-yL)/dy)+1;


vID1 = (yIndL-1)*M+xIndL;
vID3 = (yIndU-1)*M+xIndL;

vID2 = (yIndL-1)*M+xIndU;
vID4 = (yIndU-1)*M+xIndU;

if vID1 == vID2== vID3 == vID4
    %only 1 point
    elementID = [vID1,NaN];
elseif vID1 == vID3 && vID2 ==vID4
    %yUpp = yLower, on the edge of xlower and x upper, get the edgeID
    %instead of elementID
    for ele = 1:size(IENstruct,1)
        [eID, vID,~] = IENstruct(ele,:).meshData;
        if ismember([vID1,vID2],vID) == 1
            elementID = eID;
        end
    end
else
    vIDLow = [vID1,vID2,vID3];
    vIDUpp = [vID4,vID3,vID2];

    for ele = 1:size(IENstruct,1)
        [eID, vID,~] = IENstruct(ele,:).meshData;
        if isequal(vIDLow,vID) == 1
            eleLow = eID;
        end
        
        if isequal(vIDUpp,vID) == 1
            eleUpp = eID;
        end
    end

    %vIDLow
    %vIDUpp
    %eleLow
    %eleUpp
    %determine if it's lower or upper element
    k = -dy/dx;
    b = dy;
    
    Coordinate = vertexData{vID1,2};
    xLC = Coordinate(1);
    yLC = Coordinate(2);
    if (yCord - yLC)<= k*(xCord - xLC)+b
        elementID = eleLow;
    else
        elementID = eleUpp;
    end
end
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
