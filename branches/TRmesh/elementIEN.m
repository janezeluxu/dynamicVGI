function [IENall,pAll] = elementIEN(ele,IENstruct)
%get the IEN array and P list for a given element
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

