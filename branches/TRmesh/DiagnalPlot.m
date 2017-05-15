function [x,y,value] = DiagnalPlot(u,M,vertexData)
%2D figure only vertex
nVertex = size(vertexData,1);
uVertex = u(1:nVertex);

%dx = 1/M;
%dy = 1/N;

[X,Y] = meshgrid(1:M, 1:M);

Z = zeros(M,M);
for i = 1:M
    Z(:,i) = uVertex(1+(i-1)*M:i*M);
    
end

%surf(X,Y,Z)

xIndex = 1:1:M;
%xIndex = 3;
yIndex = M:-1:1;
%yIndex = 1:1:M;
vID = zeros(1,M);

for i = 1:M
        %convert double index to 1 index, relate vertex coordinated to vID
        vID(i) = (yIndex(i)-1)*M+xIndex(i);
        xv(i) = vertexData{vID(i),2}(1);
        yv(i) = vertexData{vID(i),2}(2);
end

% edge_highOrder = [];
% edgeSize = size(edgeData,1);
% for edge = 1:edgeSize
%     count = 0;
%     for vertx = 1:length(vID)
%         if any(edgeData{edge,2} == vID(vertx))
%             count = count+1;
%         end
%     end
%     
%     if (count ==2)
%         edge_highOrder = [edge_highOrder, edge];
%     end
%        
% end
% edgeDOF = [];
% count = 0;
% for i = 1: length(edge_highOrder)
%     edgeDOFs = edgeData{edge_highOrder(i),4};
%     nEdge = length(edgeDOFs);
%     for j = 1:nEdge
%         count = count+1;
%         edgeDOF(count) = edgeDOFs(j);
%         xe(count) = (-xv(i)+xv(i+1))*j/(nEdge+1)+xv(i);
%         ye(count) = (-yv(i)+yv(i+1))*j/(nEdge+1)+yv(i);
%     end
% end
% 
% uE =u(edgeDOF);
value = [u(vID)];
x = [xv];
y = [yv];
end