classdef geo_Mesh < handle
    
    properties (Access = private)
        %number of Points in x direction
        XPoints;
        %number of points in y
        YPoints;
        %left bound x coordinate
        LeftBound;
        %right bound x coordinate
        RightBound;
        %lower bound y coordinate
        LowerBound;
        %upper bound y coordinate
        UpperBound;
        
        %each block span from [xStart,yStart] to [xStop,yStop] with order P
        xStart;
        xStop;
        yStart;
        yStop;
        % Order p list
        Order;
    end
    
    methods(Access = public)
        
        function geoMesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList)
        geoMesh.XPoints = xPoints;          
        geoMesh.YPoints = yPoints;          
        geoMesh.LeftBound = xDomain(1);     
        geoMesh.RightBound = xDomain(2);    
        geoMesh.LowerBound = yDomain(1);    
        geoMesh.UpperBound = yDomain(2);    
        geoMesh.xStart = x1;                
        geoMesh.xStop = x2;                 
        geoMesh.yStart = y1;
        geoMesh.yStop = y2;
        geoMesh.Order = OrderList;          
        
        end
        
        
        function [meshData] = MeshData(geoMesh)
            %[faceID, vertexID, edgeID]
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            Grid_size = (M-1)*(N-1)*2;
            meshData = cell(1,3);
            %meshData1 = zeros(Grid_size,1);
            %meshData2 = zeros(Grid_size,3);
            %meshData3 = zeros(Grid_size,3);
            %element number and vertex
            for i  = 1:Grid_size
                meshData{i,1} = i;
                %meshData1(i) = i;
                
                meshData{i,2} = vertex(geoMesh,i);
                
                %[meshData3(i,1),startDOF] = edgeH(geoMesh,i);
                %[meshData3(i,2),startDOF] = edgeT(geoMesh,i,startDOF);
                %meshData3(i,3) = edgeV(geoMesh,i,startDOF);
                
                [meshData{i,3}(1),startDOF] = edgeH(geoMesh,i);
                [meshData{i,3}(2),startDOF] = edgeT(geoMesh,i,startDOF);
                meshData{i,3}(3) = edgeV(geoMesh,i,startDOF);
            end
            %meshData{1} = meshData1;
            %meshData{2} = meshData2;
            %meshData{3} = meshData3;
            
        end
        
        function [vertexData] = vertexMesh(geoMesh)
            %[vID, vCoordinate]
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);
            
            %vertexData = cell(1,2);
            vertexData = cell(M*N,2);
            for vID  = 1:M*N
                vertexData{vID,1} = vID;
                yIndex = fix((vID-1)/M) ;
                xIndex = mod((vID-1),M) ;
                xCoord = xIndex*dx+xL;
                yCoord = yIndex*dy+yL;
                vertexData{vID,2} = [xCoord,yCoord];
            end
        end
        
        function [edgeData] = edgeMesh(geoMesh,meshData)
            %[edgeID, VID for that edge, edgeOrder, global_edge_dof] 
            
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            p = geoMesh.Order;
            
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            Grid_size = (M-1)*(N-1)*2;
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);
            nBlock = length(p);
            
            %total number of edges
            numberEdge = max([meshData{:,3}]);
            edgeData = cell(numberEdge,4);
            meshDataEdge = cell2mat(meshData(:,3));
            
            for i = 1:Grid_size
            edgeData{meshDataEdge(i,1),2} = [meshData{i,2}(1),meshData{i,2}(2)];
            edgeData{meshDataEdge(i,2),2} = [meshData{i,2}(2),meshData{i,2}(3)];
            edgeData{meshDataEdge(i,3),2} = [meshData{i,2}(3),meshData{i,2}(1)];
            end
            
            %edgeData{1:numberEdge,1} = 1:numberEdge;
            %edgeData{:,3} = p(1);
            
            for i = 1:numberEdge
                edgeData{i,1} = i; %edgeID
                %edgeData{i,2} = geoMesh.edgeVertex(i,meshData,meshDataEdge); %vertex connect to this edge
                edgeData{i,3} = p(1); %edge Order initialize to p(1)
                %edgeData{i,4} = -1;%edge DOF initialize to -1
            end
                    
            edgeV = cell2mat(edgeData(:,2));
            %process all high order blocks
            if nBlock>1
                for i = 2:nBlock
                    pOrder = p(i);
                    %from block information, calculate affected vertex
                    vertexList =  vertex_affected(geoMesh,dx,dy,i);
                    %from the vertexList, get the high order edges, correct edge orders
                    edgeList = geoMesh.edgeList(vertexList, edgeV);
                    for edge = 1:length(edgeList)
                        edgeData{edgeList(edge),3} = max(pOrder,edgeData{edgeList(edge),3});
                    end
                end
            end
            %global edge DOF
            edgeStartCount = max(max(meshData{Grid_size,2}));
            edgeCount = 1;
            for i = 1:numberEdge
                pEdge = edgeData{i,3};
                
                if pEdge>1
                    for nEdge = 1:pEdge-1
                        edgeData{i,4} = (edgeCount-(pEdge-2):edgeCount)+ edgeStartCount;
                        edgeCount = edgeCount+1;
                    end
                end
            end
            
        end
        
        function [faceData] = faceMesh(geoMesh,meshData,edgeData)
            %[fID, fOrder, global_face_dof]
            
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            
            Grid_size = (M-1)*(N-1)*2;
            faceData = cell(Grid_size,4);
            
            for i  = 1:Grid_size
                faceData{i,1} = i;
                edgeList = meshData{i,3};
                %faceOrder is the minimum of edge orders
                faceData{i,2} = min([edgeData{edgeList(1),3},edgeData{edgeList(2),3},edgeData{edgeList(3),3}]);
                %initialize face global_DOF to -1
                %faceData{i,3} = -1;
            end
            
            %faceGlobalDOF
            faceStartCount = max([edgeData{:,4}]);
            faceCount = 1;
            for face = 1:Grid_size
                order = faceData{face,2};
                if order > 1
                    faceModes = (order-1)*(order-2)/2;
                    for nFace = 1:faceModes
                        faceData{face,3} = (faceCount-(faceModes-1):faceCount)+ faceStartCount;
                        faceCount = faceCount+1;
                    end
                end
            end
        end
        
        function [IBC,BCval] = BoundaryConditionDOF(geoMesh,edgeData)
            global BC;
            global nsd; 
            
            %All boundary DOFs list
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            
            vertexID1 = zeros(1,M);
            vertexID2 = zeros(1,M);
            vertexID3 = zeros(1,N);
            vertexID4 = zeros(1,N);
            
            %bottom and top vertex
            if(nsd~=1)
                for y = 1: M
                    %x = 1
                    vertexID1(y) =  y;
                    %x = N
                    vertexID2(y) =   (N-1)*M+y;
                end
            else
                vertexID1 = []; 
                vertexID2 = [];                 
            end
            
            %left and right vertex
            for x = 1: N
                vertexID3(x) =   (x-1)*M+1;
                vertexID4(x) =   (x-1)*M+M;
            end
            
            if(nsd==1)
                vertexID = [vertexID3,vertexID4];
            else
                vertexID = [vertexID1,vertexID2,vertexID3(2:end-1),vertexID4(2:end-1)];                           
            end
            
            %numberEdge = length(edgeData);
            edgeV = cell2mat(edgeData(:,2));
            
            if(nsd==1)
                edge1 = [];
                edge2 = []; 
            else
                edge1 = geoMesh.edgeList( vertexID1, edgeV);
                edge2 = geoMesh.edgeList( vertexID2, edgeV);
            end            
            edge3 = geoMesh.edgeList( vertexID3, edgeV);
            edge4 = geoMesh.edgeList( vertexID4, edgeV);
            
            edgeDOF = [];
            edge = [edge1,edge2,edge3,edge4];
            for i = 1: length(edge)
                edgeDOF = [edgeDOF,edgeData{edge(i),4}];
            end
            IBC = [vertexID,edgeDOF];
            
            BCval = zeros(length(IBC),1)*NaN; 
            count = 1; 
            if(nsd~=1)
                % Bottom
                for i = 1:length(vertexID1)
                    if(i==1)
                        BCval(count) = max(BC(2),BC(3)); 
                    elseif (i==length(vertexID1))
                        BCval(count) = max(BC(2),BC(4));                     
                    else
                        BCval(count) = BC(2);
                    end
                    count = count+1; 
                end
                % Top
                for i = 1:length(vertexID2)
                    if(i==1)
                        BCval(count) = max(BC(1),BC(3)); 
                    elseif (i==length(vertexID2))
                        BCval(count) = max(BC(1),BC(4));                     
                    else
                        BCval(count) = BC(1);
                    end                
                    count = count+1; 
                end 
                % Left
                for i = 2:(length(vertexID3)-1)
                    BCval(count) = BC(3);
                    count = count+1; 
                end            
                % Right
                for i = 2:(length(vertexID4)-1)
                    BCval(count) = BC(4);
                    count = count+1; 
                end    
            else
                % Left
                for i = 1:(length(vertexID3))
                    BCval(count) = BC(3);
                    count = count+1; 
                end            
                % Right
                for i = 1:(length(vertexID4))
                    BCval(count) = BC(4);
                    count = count+1; 
                end      
            end
                  
            
            % Edge Nodes
            if(length(edgeDOF)>0)
                if(nsd==1)
                    for i = 1:length(edge1)
                        BCval(count) = BC(2);
                        count = count+1; 
                    end
                    for i = 1:length(edge2)
                        BCval(count) = BC(1);
                        count = count+1; 
                    end     
                end
                for i = 1:length(edge3)
                    BCval(count) = BC(3);
                    count = count+1; 
                end            
                for i = 1:length(edge4)
                    BCval(count) = BC(4);
                    count = count+1; 
                end         
            end
        end        
        
    end
    methods(Access = private)
        
        
        function [vertexID] = vertex(geoMesh,ele)
            %for each odd element, calculate the vertexID
            M = geoMesh.XPoints;
            if mod(ele,2) == 1
                Reminder = mod(ele,2*(M-1));
                result = fix(ele/(2*(M-1)));
                c = fix(Reminder/2);
                a = result*M+c+1;
                vertexID = [a,a+1,a+M];
            %for each even element, calculate the vertexID 
            else
                ele = ele-1;
                Reminder = mod(ele,2*(M-1));
                result = fix(ele/(2*(M-1)));
                c = fix(Reminder/2);
                a = result*M+c+2;
                vertexID = [a+M,a+M-1,a];
            end
        end
        
        
        function [edgeMesh,lastDOF] = edgeH(geoMesh,ele)
            %for each odd element, calculate the Horizontal edgeID
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            if mod(ele,2) == 1
                m = (ele+1)/2;
                edgeMesh =  m;
                lastDOF = (M-1)*N;
            %for each even element, calculate the Horizontal edgeID    
            else
                n = (ele-1) + (M-1)*2;
                m = (n+1)/2;
                edgeMesh =  m;
                lastDOF = (M-1)*N;
            end
            
            
        end
        
        function [edgeMesh,lastDOF] = edgeT(geoMesh,ele,startDOF)
            %for each odd element, calculate the tilted edgeID
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            if mod(ele,2) == 1
                m = (ele+1)/2;
                edgeMesh =  m+startDOF;
                lastDOF = (M-1)*N+(M-1)*(N-1);
            %for each even element, calculate the tilted edgeID    
            else
                m = (ele)/2;
                edgeMesh =  m+startDOF;
                lastDOF = (M-1)*N+(M-1)*(N-1);
            end
        end
        
        
        function [edgeMesh] = edgeV(geoMesh,ele,startDOF)
            %for each odd element, calculate the vertical edgeID
            M = geoMesh.XPoints;
            if mod(ele,2) == 1
                Reminder = mod(ele,(2*(M-1)));
                nRow = fix(ele/(2*(M-1)));
                c = fix(Reminder/2);
                m = nRow*(M)+c+1;
                edgeMesh =  m+startDOF;
            %for each even element, calculate the vertical edgeID    
            else
                Reminder = mod(ele-1,(2*(M-1)));
                nRow = fix((ele-1)/(2*(M-1)));
                c = fix(Reminder/2);
                m = nRow*(M)+c+2;
                edgeMesh =  m+startDOF;
            end
        end
        
        
        
        function [vertexNumber] = vertex_affected(geoMesh,dx,dy,i)
            %from high order vertex coordinate, get the vID
            x1 = geoMesh.xStart;
            y1 = geoMesh.yStart;
            x2 = geoMesh.xStop;
            y2 = geoMesh.yStop;
            xL = geoMesh.LeftBound;
            yL = geoMesh.LowerBound;
            M = geoMesh.XPoints;
            
            xS = x1(i);
            yS = y1(i);
            xT = x2(i);
            yT = y2(i);
            
            x1Index = fix((xS-xL)/dx) +1;
            y1Index = fix((yS-yL)/dy) +1;
            
            x2Index = fix((xT-xL)/dx) +1;
            y2Index = fix((yT-yL)/dy) +1;
            
            vertexNumber = zeros(1,(x2Index-x1Index)*(y2Index-y1Index));
            
            i = 0;
            for y = y1Index:y2Index
                for x = x1Index:x2Index
                    i = i+1;
                    %convert double index to 1 index
                    vertexNumber(i) = (y-1)*M+x;
                end
            end
            
        end
        



    end
    methods(Static)
        
        function [edgeVertexList] = edgeVertex(edgeID, meshData,meshDataEdge)
            %from the edgeID, find the two vertexID
            [rowID,columnID] = find((edgeID==meshDataEdge),1);
            edgeVertexList = [meshData{rowID,2}(columnID),meshData{rowID,2}(mod(columnID,3)+1)];
            
%             for ele = 1:Grid_size
%                 if meshData{ele,3}(1) == edgeID;
%                     edgeVertexList = [meshData{ele, 2}(1),meshData{ele, 2}(2)];
%                 elseif meshData{ele,3}(2) == edgeID;
%                     edgeVertexList = [meshData{ele, 2}(2),meshData{ele, 2}(3)];
%                 elseif meshData{ele,3}(3) == edgeID;
%                     edgeVertexList = [meshData{ele, 2}(3),meshData{ele, 2}(1)];
%                 end
%             end
        end
        
        
        function [edge_highOrder] = edgeList(vertex_list, edgeV)
            %from vertexList, find the neighboring edge List
            
            E1List = [];
            E2List = [];
            %edgeV = cell2mat(edgeData(:,2));
            
            for i = 1:length(vertex_list)
                E1List = [E1List;find(vertex_list(i) == edgeV(:,1))];
                E2List = [E2List;find(vertex_list(i) == edgeV(:,2))];
            end
            [~,LocE2] = ismember(E1List,E2List);
            LocE2 = LocE2(LocE2~=0);
            
            edge_highOrder = zeros(length(LocE2),1);
            for i = 1:length(LocE2)
                edge_highOrder(i) = E2List(LocE2(i));
            end
            
            edge_highOrder = sort(edge_highOrder');
            %LocE2
            
%             edge_highOrder = [];
%             for edge = 1:edgeSize
%                 count = 0;
%                 for vertx = 1:length(vertex_list)
%                     if any(edgeData{edge,2} == vertex_list(vertx))
%                         count = count+1;
%                     end
%                 end
%                 
%                 if (count ==2)
%                     edge_highOrder = [edge_highOrder, edge];
%                 end
%                 
%             end
            
        end
        
    end
end