classdef geo_Mesh_symm < handle
    
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
        PV;
        PH;
        PT

    end
    
    methods(Access = public)
        
        function geoMesh = geo_Mesh_symm(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,varargin)
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
        
        if length(varargin) ==1
            geoMesh.Order = varargin{1}; 
            geoMesh.PV = nan;
            geoMesh.PH = nan;
            geoMesh.PT = nan;
        elseif length(varargin) ==3
            geoMesh.PV = varargin{1};
            geoMesh.PH = varargin{2};
            geoMesh.PT = varargin{3};
            geoMesh.Order = nan;
        else
          disp('number of input is wrong!')
        end
        end
        
        
        function [meshData] = MeshData(geoMesh)
            %[faceID, vertexID, edgeID]
            
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            Patch_size = (M-1)*(N-1);
            meshData = cell(Patch_size*4,3);
            
            meshData1 = zeros(Patch_size,1);
            meshData2 = zeros(Patch_size,5);
            meshData3 = zeros(Patch_size,8);
            
            %element number and vertex
            for i  = 1:Patch_size
                meshData1(i) = i;  
                meshData2(i,:) = vertex(geoMesh,i);
                meshData3(i,:) = edge(geoMesh,i);
                
                meshData{(i-1)*4+1,1} = (i-1)*4+1;
                meshData{(i-1)*4+2,1} = (i-1)*4+2;
                meshData{(i-1)*4+3,1} = (i-1)*4+3;
                meshData{(i-1)*4+4,1} = (i-1)*4+4;
                
                meshData{(i-1)*4+1,2} = [meshData2(i,1),meshData2(i,2),meshData2(i,5)];
                meshData{(i-1)*4+1,3} =  [meshData3(i,6),meshData3(i,5),meshData3(i,1)];
                meshData{(i-1)*4+2,2} = [meshData2(i,2),meshData2(i,3),meshData2(i,5)];
                meshData{(i-1)*4+2,3} =  [   meshData3(i,7),meshData3(i,6),meshData3(i,2)];
                meshData{(i-1)*4+3,2}  = [meshData2(i,5),meshData2(i,3),meshData2(i,4)];
                meshData{(i-1)*4+3,3}  =   [  meshData3(i,3),meshData3(i,8),meshData3(i,7)];
                meshData{(i-1)*4+4,2}  = [meshData2(i,1),meshData2(i,5),meshData2(i,4)];
                meshData{(i-1)*4+4,3}  =   [  meshData3(i,8),meshData3(i,4),meshData3(i,5)];
            end
        end
        %meshData{1}(i,:)
        function [vertexData] = vertexMesh(geoMesh,meshData)
            %[vID, vCoordinate, element around vertex, vertex around vertex]
            %global Grid_size
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            Grid_size = (M-1)*(N-1)*4;
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);
            
            %vertexData = cell(1,2);
            %edgeVertex = cell2mat(edgeData(:,2));
            vertexData = cell(M*N+(M-1)*(N-1),2);
            %vertex1 = zeros(M*N+(M-1)*(N-1),1);
            %vertex2 = zeros(M*N+(M-1)*(N-1),2);
            
            
            vertexFace = cell2mat(meshData(:,2));
            for i = 1:length(vertexData)
                [faceID] = geoMesh.FaceVertex(i,vertexFace);
                vertexData{i,3} = faceID'; %edge face relation
                vvID = [];
                for j = 1:length(faceID)
                    vvID = [vvID,meshData{faceID(j),2}];
                end
                vertexData{i,4} = vvID;
            end
                
            for vID  = 1:M*N
                vertexData{vID,1}= vID;
                yIndex = fix((vID-1)/M) ;
                xIndex = mod((vID-1),M) ;
                xCoord = xIndex*dx+xL;
                yCoord = yIndex*dy+yL;
                vertexData{vID,2} = [xCoord,yCoord];
            end
            dx = (xR-xL)/(2*(M-1));
            dy = (yU-yL)/(2*(N-1));
            for vID  = M*N+1:M*N+(M-1)*(N-1)
                vertexData{vID,1}= vID;
                vv = vID-M*N;
                a = fix((vv-1)/(M-1))+1;
                b = vv-(a-1)*(M-1);
                xIndex = 2*b-1;
                yIndex = 2*a-1;
                xCoord = xIndex*dx+xL;
                yCoord = yIndex*dy+yL;
                
                vertexData{vID,2} = [xCoord,yCoord];
            end
            %vertexData
        end
        
        function [edgeData] = edgeMesh(geoMesh,meshData)
            if isnan(geoMesh.Order) == 0
                edgeData = edgeAllMesh(geoMesh,meshData);
            else
                edgeData = edgeVHTMesh(geoMesh,meshData);
            end
            
        end
        function [edgeData] = edgeAllMesh(geoMesh,meshData)
            %[edgeID, VID for that edge, faceID, edgeOrder, global_edge_dof] 
            
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            p = geoMesh.Order;
            
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            Grid_size = (M-1)*(N-1)*4;
            %Patch_size*4
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);
            nBlock = length(p);
            
            %total number of edges
            numberEdge = max([meshData{:,3}]);
            edgeData = cell(numberEdge,5);
            meshDataEdge = cell2mat(meshData(:,3));
            for i = 1:Grid_size
                %edge2(meshDataEdge(i,1),:) = [meshData{2}(i,1), meshData{2}(i,2)];
                %edge2(meshDataEdge(i,2),:) = [meshData{2}(i,2), meshData{2}(i,3)];
                %edge2(meshDataEdge(i,3),:) = [meshData{2}(i,3), meshData{2}(i,1)];
                
            edgeData{meshDataEdge(i,1),2} = [meshData{i,2}(2),meshData{i,2}(3)];
            edgeData{meshDataEdge(i,2),2} = [meshData{i,2}(3),meshData{i,2}(1)];
            edgeData{meshDataEdge(i,3),2} = [meshData{i,2}(1),meshData{i,2}(2)];
            end
            
            edgeFace = cell2mat(meshData(:,3));
            for i = 1:numberEdge
                edgeData{i,1} = i; %edgeID
                [faceID] = geoMesh.FaceEdge(i,edgeFace);
                edgeData{i,3} = faceID'; %edge face relation
                edgeData{i,4} = p(1); %edge Order initialize to p(1)
                %edgeData{i,4} = -1;%edge DOF initialize to -1
            end
            
            edgeVertex = cell2mat(edgeData(:,2));
            %process all high order blocks
            if nBlock>1
                for i = 2:nBlock
                    pOrder = p(i);
                    %from block information, calculate affected vertex
                    vertexList =  vertex_affected(geoMesh,dx,dy,i);
                    %from the vertexList, get the high order edges, correct edge orders
                    edgeListA = geoMesh.edgeList(vertexList, edgeVertex);
                    
                    
                    for edge = 1:length(edgeListA)
                        edgeData{edgeListA(edge),4} = max(pOrder,edgeData{edgeListA(edge),4});
                    end
                end
            end
            
            %global edge DOF
            edgeStartCount = max(max(meshData{Grid_size,2}));
            edgeCount = 1;
            for i = 1:numberEdge
                pEdge = edgeData{i,4};
                
                if pEdge>1
                    for nEdge = 1:pEdge-1
                        edgeData{i,5} = (edgeCount-(pEdge-2):edgeCount)+ edgeStartCount;
                        edgeCount = edgeCount+1;
                    end
                end
            end
            
        end
        
        function [faceData] = faceMesh(geoMesh,meshData,edgeData)
            %[fID, fOrder, global_face_dof,faceArea]
            
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            TotalArea = (xR-xL)*(yU-yL);
            Grid_size = (M-1)*(N-1)*4;
            faceData = cell(Grid_size,4);
            
            for i  = 1:Grid_size
                faceData{i,1} = i;
                edgeList = meshData{i,3};
                %faceOrder is the minimum of edge orders
                faceData{i,2} = max([edgeData{edgeList(1),4},edgeData{edgeList(2),4},edgeData{edgeList(3),4}]);
                %initialize face global_DOF to -1
                faceData{i,4} = TotalArea/Grid_size;
            end
            
            %faceGlobalDOF
            faceStartCount = max([edgeData{:,5}]);
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
        
        function [IBC1,IBC2,IBC3,IBC4,V,BCval1,BCval2,BCval3,BCval4,VBCval] = BoundaryCondition(geoMesh,edgeData)
            %All boundary DOFs list
             global BC;
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            
            vertexID1 = zeros(1,M);
            vertexID2 = zeros(1,M);
            vertexID3 = zeros(1,N);
            vertexID4 = zeros(1,N);
            
            %bottom and top vertex
            for y = 1: M
                %x = 1
                vertexID1(y) =  y;
                %x = N
                vertexID2(y) =   (N-1)*M+y;
            end
            
            %left and right vertex
            for x = 1: N
                vertexID3(x) =   (x-1)*M+1;
                vertexID4(x) =   (x-1)*M+M;
            end
            
            %vertexID = [vertexID1,vertexID2,vertexID3(2:end-1),vertexID4(2:end-1)];
            
            %vertexID = [vertexID3,vertexID4];
            edgeV = cell2mat(edgeData(:,2));
            
            edge1 = geoMesh.edgeList( vertexID1, edgeV);
            edge2 = geoMesh.edgeList( vertexID2, edgeV);
            edge3 = geoMesh.edgeList( vertexID3, edgeV);
            edge4 = geoMesh.edgeList( vertexID4, edgeV);
            
            edgeDOF1 = [];
            edgeDOF2 = [];
            edgeDOF3 = [];
            edgeDOF4 = [];
            
            for i = 1: length(edge1)
                edgeDOF1 = [edgeDOF1,edgeData{edge1(i),5}];
            end
            
            for i = 1: length(edge2)
                edgeDOF2 = [edgeDOF2,edgeData{edge2(i),5}];
            end
            for i = 1: length(edge3)
                edgeDOF3 = [edgeDOF3,edgeData{edge3(i),5}];
            end
            for i = 1: length(edge4)
                edgeDOF4 = [edgeDOF4,edgeData{edge4(i),5}];
            end
            
            IBC1 = [vertexID1(2:end-1),edgeDOF1];
            IBC2 = [vertexID2(2:end-1),edgeDOF2];
            IBC3 = [vertexID3(2:end-1),edgeDOF3];
            IBC4 = [vertexID4(2:end-1),edgeDOF4];
            
            BCval1 = zeros(1,length(IBC1))*NaN; 
            BCval2 = zeros(1,length(IBC2))*NaN; 
            BCval3 = zeros(1,length(IBC3))*NaN; 
            BCval4 = zeros(1,length(IBC4))*NaN; 
            
            % Bottom
            count = 1;
            for i = 2:length(vertexID1)-1
                BCval1(count) = BC(1);
                count = count+1;
            end
            for i = 1:length(edgeDOF1)
                BCval1(count) = BC(1);
                count = count+1;
            end
            
            % Top
            count = 1;
            for i = 2:length(vertexID2)-1
                BCval2(count) = BC(2);
                count = count+1;
            end
            for i = 1:length(edgeDOF2)
                BCval2(count) = BC(2);
                count = count+1;
            end
            % Left
            count = 1;
            for i = 2:(length(vertexID3))-1                
                BCval3(count) = BC(3);
                count = count+1;
            end
            for i = 1:length(edgeDOF3)
                    BCval3(count) = BC(3);
                    count = count+1; 
            end 
            % Right
            count = 1;
            for i = 2:(length(vertexID4))-1
                BCval4(count) = BC(4);
                count = count+1;
            end
            for i = 1:length(edgeDOF4)
                BCval4(count) = BC(4);
                count = count+1;
            end
            
            V = [vertexID1(1),vertexID1(end),vertexID2(1),vertexID2(end)];
            VBCval = [max(BC(1),BC(3)),max(BC(1),BC(4)),max(BC(2),BC(3)),max(BC(2),BC(4))];
            
        end
        
        function [edgeData] = edgeVHTMesh(geoMesh,meshData)
            %[edgeID, VID for that edge, edgeOrder, global_edge_dof] 
            
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            pV = geoMesh.PV;
            pH = geoMesh.PH;
            pT = geoMesh.PT;
           
            
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            Grid_size = (M-1)*(N-1)*2;
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);
            nBlock = length(pV);
            
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
                edgeData{i,3} = min([pV(1),pH(1),pT(1)]); %edge Order initialize to p(1)
                %edgeData{i,4} = -1;%edge DOF initialize to -1
            end
                    
            edgeVertex = cell2mat(edgeData(:,2));
            %process all blocks
            
                for i = 1:nBlock
                    %from block information, calculate affected vertex
                    [vertexList,xDis,yDis] =  vertex_affected(geoMesh,dx,dy,i);
                    %from the vertexList, get the high order edges, correct edge orders
%                     
                     edgeListA = geoMesh.edgeList(vertexList, edgeVertex);
                     edgeHList = edgeListA(1:xDis*(yDis+1));
                     edgeTList = edgeListA(xDis*(yDis+1)+1:xDis*(yDis+1)+xDis*yDis);
                     edgeVList = edgeListA(xDis*(yDis+1)+xDis*yDis+1:end);
                     
                    for edge = 1:length(edgeVList)
                        edgeData{edgeVList(edge),3} = max(pV(i),edgeData{edgeVList(edge),3});
                    end
                    for edge = 1:length(edgeHList)
                        edgeData{edgeHList(edge),3} = max(pH(i),edgeData{edgeHList(edge),3});
                    end
                    for edge = 1:length(edgeTList)
                        edgeData{edgeTList(edge),3} = max(pT(i),edgeData{edgeTList(edge),3});
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
        
        
        
        
    end
    methods(Access = private)
        
        
        function [vertexID] = vertex(geoMesh,patch)
            %for each odd element, calculate the vertexID
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            nRow = fix((patch-1)/(M-1));
            nColumn = mod(patch-1,((M-1)))+1;
            c = fix(nColumn);
            a = nRow*(M)+c;
            vertexID = [a,a+1,a+M+1,a+M,M*N+patch];
        end
        
        
        function [edgeMesh] = edge(geoMesh,patch)
            %for each odd element, calculate the Horizontal edgeID
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            
            nRow = fix((patch-1)/(M-1));
            nColumn = mod(patch-1,((M-1)))+1;
            c = fix(nColumn);
            a = nRow*(M)+c;
            
            %count = patch+2*M*(N-1);
            edgeMesh = [patch,M*(N-1)+a+1, patch+M-1,M*(N-1)+a,...
                (patch-1)*4+2*M*(N-1)+1,(patch-1)*4+2*M*(N-1)+2,(patch-1)*4+2*M*(N-1)+3,(patch-1)*4+2*M*(N-1)+4];
            
            
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
        
        
        
        function [vertexNumber,xDis,yDis] = vertex_affected(geoMesh,dx,dy,i)
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
            
            xDis = x2Index - x1Index;
            yDis = y2Index - y1Index;
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
            
        end
        
        function [faceID] = FaceEdge(edgeID,edgeFace)
            %from edgeID, find the two/one face 
            [faceID,~] = find(edgeID == edgeFace);
            %faceID = row;
        end
        
        function [faceID] = FaceVertex(vID,vertexFace)
            %from edgeID, find the two/one face 
            [faceID,~] = find(vID == vertexFace);
            %faceID = row;
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
            
            
        end
        
    end
end