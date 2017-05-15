classdef IEN < handle
    
    properties (Access = private)
        %input from geo_Mesh class
        MeshData; 
        EdgeData; 
        FaceData;
        %number of points in x and y 
        XPoints;
        YPoints;
        
    end
    
    methods(Access = public)
        function ien = IEN(xPoints,yPoints,Mesh,Edge,Face)
            ien.XPoints = xPoints;    
            ien.YPoints = yPoints;
            ien.MeshData = Mesh;
            ien.EdgeData = Edge;
            ien.FaceData = Face;
            
        end
        function [vertexIEN,edgeIEN,faceIEN] = Construct_IEN(ien)
            %vertexIEN[face, vertexDOF]
            %edgeIEN [edge1Order, edge1DOF],[edge2Order,
            %edge2DOF],[edge3Order, edge3DOF]
            %faceIEN [faceID] [faceOrder, faceDOF]
            global Grid_size
            meshData = ien.MeshData;
            edgeData = ien.EdgeData;
            faceData = ien.FaceData;
            M = ien.XPoints;
            N = ien.YPoints;
            
            %Grid_size = (M-1)*(N-1)*4;
            
            vertexIEN = cell(Grid_size,3);
            edgeIEN = cell(Grid_size,3);
            faceIEN = cell(Grid_size,3);
            
            %loop over all elements, find the correspond vID, edgeID, faceID
            for i  = 1:Grid_size
                vertexIEN{i,1} = [meshData{i,1}];
                vertexIEN{i,2} = [meshData{i,2}];
                
                eID = meshData{i,3};
                if mod(i,2) == 1
                    edgeIEN{i,1} = [edgeData{eID(1),4}, edgeData{eID(1),5}];
                    edgeIEN{i,2} = [edgeData{eID(2),4},edgeData{eID(2),5}];
                    edgeIEN{i,3} = [edgeData{eID(3),4}, edgeData{eID(3),5}];
                %face uses edge in reverse order, will add usage
                %information later
                else
                    edgeIEN{i,1} = [edgeData{eID(1),4}, fliplr(edgeData{eID(1),5})];
                    edgeIEN{i,2} = [edgeData{eID(2),4}, fliplr(edgeData{eID(2),5})];
                    edgeIEN{i,3} = [edgeData{eID(3),4}, fliplr(edgeData{eID(3),5})]; 
                end
               
                fID = i;
                faceIEN{i,1} = [meshData{i,1}];
                faceIEN{i,2} = [faceData{fID,2},faceData{fID,3}];
                
            end
            
            %IENstruct = struct('vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN)
            
        end

    end
    
    
end