classdef geo_Mesh_Square < handle
    
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
        
        function geoMesh = geo_Mesh_Square(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,varargin)
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
            global Grid_size
            meshData = cell(Grid_size,3);
            %element number and vertex
            for i  = 1:Grid_size
                meshData{i,1} = i;               
                meshData{i,2} = vertex(geoMesh,i);
                [Hedge,startDOF] = edgeH(geoMesh,i);
                
                vEdge = edgeV(geoMesh,i,startDOF);
                meshData{i,3} = [Hedge(1),vEdge(2),Hedge(2),vEdge(1)];
            end
        end
        
        function [vertexData] = vertexMesh(geoMesh,meshData)
            %[vID, vCoordinate, element around vertex, vertex around vertex]
            %global Grid_size
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);

            vertexData = cell(M*N,2);
            for vID  = 1:M*N
                vertexData{vID,1} = vID;
                yIndex = fix((vID-1)/M) ;
                xIndex = mod((vID-1),M) ;
                xCoord = xIndex*dx+xL;
                yCoord = yIndex*dy+yL;
                vertexData{vID,2} = [xCoord,yCoord];
            end
            
            vertexFace = cell2mat(meshData(:,2));
            for i = 1:length(vertexData)
                [faceID] = geoMesh.FaceVertex(i,vertexFace);
                vertexData{i,3} = unique(faceID'); %edge face relation
                vvID = [];
                for j = 1:length(faceID)
                    vvID = [vvID,meshData{faceID(j),2}];
                end
                vertexData{i,4} = unique(vvID);
            end
        end
        
        
        
        function [edgeData] = edgeMesh(geoMesh,meshData)
            if isnan(geoMesh.Order) == 0
                edgeData = edgeAllMesh(geoMesh,meshData);
            else
                edgeData = edgeVHTMesh(geoMesh,meshData);
            end
            
        end
        
        function [edgeData] = edgeAllMesh(geoMesh,meshData)
            %[edgeID, VID for that edge, edgeOrder, global_edge_dof] 
            global Grid_size
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            p = geoMesh.Order;
            
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            %Grid_size =  (M-1)*(N-1);
            dx = (xR-xL)/(M-1);
            dy = (yU-yL)/(N-1);
            nBlock = length(p);
            
            %total number of edges
            numberEdge = max([meshData{:,3}]);
            edgeData = cell(numberEdge,5);
            meshDataEdge = cell2mat(meshData(:,3));
            
            for i = 1:Grid_size
            edgeData{meshDataEdge(i,1),2} = [meshData{i,2}(1),meshData{i,2}(2)];
            edgeData{meshDataEdge(i,2),2} = [meshData{i,2}(2),meshData{i,2}(3)];
            edgeData{meshDataEdge(i,3),2} = [meshData{i,2}(3),meshData{i,2}(4)];
            edgeData{meshDataEdge(i,4),2} = [meshData{i,2}(4),meshData{i,2}(1)];
            end
            edgeFace = cell2mat(meshData(:,3));
            for i = 1:numberEdge
                edgeData{i,1} = i; %edgeID
                [faceID] = geoMesh.FaceEdge(i,edgeFace);
                edgeData{i,3} = faceID; %edge face relation
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
            global Grid_size
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            xL = geoMesh.LeftBound;
            xR = geoMesh.RightBound;
            yL = geoMesh.LowerBound;
            yU = geoMesh.UpperBound;
            
            TotalArea = (xR-xL)*(yU-yL);
            %Grid_size = (M-1)*(N-1);
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
        
        function [IBC1,IBC2,IBC3,IBC4,V,BCval1,BCval2,BCval3,BCval4,VBCval] = BoundaryCondition(geoMesh,edgeData,vertexData)
            %All boundary DOFs list
            global BC;
            global StepFunc;
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

            [BCval1] = geoMesh.AssignBCval(1,vertexID1,edgeDOF1,BCval1,vertexData);
            [BCval2] = geoMesh.AssignBCval(2,vertexID2,edgeDOF2,BCval2,vertexData);
            [BCval3] = geoMesh.AssignBCval(3,vertexID3,edgeDOF3,BCval3,vertexData);
            [BCval4] = geoMesh.AssignBCval(4,vertexID4,edgeDOF4,BCval4,vertexData);
            
            V = [vertexID1(1),vertexID1(end),vertexID2(1),vertexID2(end)];
            %VBCval = [0,0,0,0];
            VBCval = [max(BC{1},StepFunc(1)),max(BC{1},BC{4}),max(BC{2},StepFunc(2)),max(BC{2},BC{4})];
            %VBCval = [max(BC{1},BC{3}),max(BC{1},StepFunc(1)),max(BC{2},BC{3}),max(BC{2},StepFunc(2))];
            
        end
        
    end
    methods(Access = private)
        
        
        function [vertexID] = vertex(geoMesh,ele)
            %for each odd element, calculate the vertexID
            [row,column] = geoMesh.elementDouble(ele);
            M = geoMesh.XPoints;
            a = (row-1)*M+column;
            vertexID = [a,a+1,a+1+M,a+M];
        end
        
        
        function [edgeMesh,lastDOF] = edgeH(geoMesh,ele)
            %for each odd element, calculate the Horizontal edgeID
            M = geoMesh.XPoints;
            N = geoMesh.YPoints;
            edgeMesh =  [ele,ele+M-1];
            lastDOF = (M-1)*N;            
        end
       
        function [edgeMesh] = edgeV(geoMesh,ele,startDOF)
            %for each odd element, calculate the vertical edgeID
            [vertexID] = geoMesh.vertex(ele);
            a = vertexID(1);
            edgeMesh =  [a+startDOF,a+startDOF+1];
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
        
        function [a,b] = elementDouble(geoMesh,n)
            M = geoMesh.XPoints-1;
            a = fix((n-1)/(M))+1;%row
            b = n-M*(a-1);%column
        end
        
        function [a,b] = vertexDouble(geoMesh,n)
            M = geoMesh.XPoints;
            a = fix((n-1)/(M))+1;%row
            b = n-M*(a-1);%column
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
        
        function [BCval] = AssignBCval(BCI,vertexID,edgeDOF,BCval,vertexData)
            global BC;
            global StepFunc;
            if (strcmp(BC(BCI),'StepFunc')==1)
                disp ('StepFunc!!!');
                count = 1;
                for i = 2:(length(vertexID))-1
                    vcord = vertexData{vertexID(i),2};
                    
                    if vcord(2)<StepFunc(3)%only consider stepfunction in y direction
                        BCval(count) = StepFunc(1);
                    else
                        BCval(count) = StepFunc(2);
                    end
                    count = count+1;
                end
            else            
            count = 1;
            for i = 2:(length(vertexID))-1
                BCval(count) = BC{BCI};
                count = count+1;
            end
            for i = 1:length(edgeDOF)
                BCval(count) = BC{BCI};
                count = count+1;
            end
            end
        end
        
    end
end