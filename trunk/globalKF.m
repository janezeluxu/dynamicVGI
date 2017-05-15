function [kGlobal, fGlobal] = globalKF(p,IENstruct,IBC,BCval,...
                                vertexData,M,N,iper)
% Assemble of element level matrix   
global TotalDOF; 

%create uniform p shape function tables, p is a list
ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    n = nIntergerPoints(p(j));
    qPoints = TriGaussPoints(n);
    simplexsf = SimplexShapeFunc(qPoints,p(j));
    [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();    
    ShapeFuncTable{p(j)} =  ShapeFunc;
    divSFtable{p(j)} = divShapeFunc;
end

Grid_size = (M-1)*(N-1)*2;

row = [];
column = [];
value = [];
fGlobal = zeros(TotalDOF,1);
Kindex = ones(TotalDOF,TotalDOF)*-1;
%Kindex = spalloc(TotalDOF,TotalDOF,TotalDOF*200);
keyCount = 0;
for ele = 1:Grid_size
%for ele = 1:1    
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions    
    [IENall,pAll] = elementIEN(ele,IENstruct);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll));
    
    if range(pAll) == 0
        %disp('uniform p element');
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
    else
        %disp('nonuniform p element');
        qPoints = TriGaussPoints(n);
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, divShapeFunc ] = simplexsf.variablePShapeFunctionTable();
    
    end
    
    %find element Shape Functions, insert stiffness matrix and force vectors
    %into global matrixs(spars)
    %simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,kappa,q,ShapeFunc,divShapeFunc);
    simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,ShapeFunc,divShapeFunc);

    [elementK,elementF] = simplex.eleStiffMatrix();
    
    %sparse matrix input row, column, and value list, using row and column
    %combination as key
    for i = 1:nssl
        rowTemp = IENall(i);
        fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
        for j = 1:nssl           
            columnTemp = IENall(j);
            %from a given key, find the value index, -1 means it is a new
            %key, thus insert rowTemp, columnTemp, value at the end of the
            %list
            indexKglobal = Kindex(rowTemp,columnTemp);
            %indexKglobal = IENindex(rowTemp, columnTemp, row, column);
            
            if indexKglobal == -1
                keyCount = keyCount+1; %its a new key
                row = [row, rowTemp]; 
                column = [column, columnTemp]; 
                value = [value,elementK(i,j)]; 
                Kindex(rowTemp,columnTemp) = keyCount; %update Kindex
            else
                value(indexKglobal) = elementK(i,j)+value(indexKglobal);
            end
        end
    end
      
end

% Account for BC in forcing vector. 
for i = 1:length(IBC)
    for j = 1:TotalDOF
        % XXX use ien instead of ismember and TotalDOF, will be faster
        if(~ismember(j,IBC) || BCval(i)==0)    % Don't do this if there's a BC here, it'll be 0 anyways
            tmpIndex = IENindex(j,IBC(i),row,column);
            %tmpIndex = Kindex(j,IBC(i));
            if(tmpIndex>0)
                fGlobal(j) = fGlobal(j) - BCval(i)*value(tmpIndex); 
            end
        end       
    end
end

% Account for periodicity 
for i = 1:length(iper)
    if(iper(i) ~= i)
        rowTemp = iper(i);

        % forcing vector
        fGlobal(i) = fGlobal(i) + fGlobal(rowTemp);
        fGlobal(rowTemp) = fGlobal(i); 
        %fGlobal(i) = 0; 
        
        % stiffness matrix
        for j = 1:TotalDOF
            IndexSlave = Kindex(i,j);
            IndexMaster = Kindex(rowTemp,j);            
            if(IndexMaster == -1 && IndexSlave~=-1) 
                row = [row, rowTemp]; 
                column = [column, j]; 
                value = [value,value(IndexSlave)]; 
            elseif(IndexSlave == -1 && IndexMaster ~= -1)
                row = [row, i]; 
                column = [column, j]; 
                value = [value,value(IndexMaster)];                
            elseif(IndexSlave ~= -1 && IndexMaster ~= -1)
                value(IndexMaster) = value(IndexMaster) + value(IndexSlave); 
                value(IndexSlave) = value(IndexMaster); 
            end
%             if IndexSlave ~= -1
%                 %%delete indexSlave
%                 value(IndexSlave)=[];
%                 row(IndexSlave)=[];
%                 column(IndexSlave)=[];
%             end
%             
%             if i==j
%                 row = [row, IndexSlave]; 
%                 column = [column, IndexSlave]; 
%                 value = [value,1]; 
%             end
        end
        
                
    end    
end

% Account for periodicity 
% iper
% fGlobal
% row
% column
% value

% for i = 1:length(iper)
%     if(iper(i) ~= i)
%         rowTemp = iper(i);
% 
%         % forcing vector
%         fGlobal(rowTemp) = fGlobal(i) + fGlobal(rowTemp);
%         fGlobal(i) = 0; 
%         
%         % stiffness matrix
%         for j = 1:TotalDOF
%             IndexSlave = Kindex(i,j);
%             IndexMaster = Kindex(rowTemp,j);            
%             if(IndexMaster == -1 && IndexSlave~=-1) 
%                 row = [row, rowTemp]; 
%                 column = [column, j]; 
%                 value = [value,value(IndexSlave)]; 
%             elseif(IndexSlave == -1 && IndexMaster ~= -1)
%                 row = [row, i]; 
%                 column = [column, j]; 
%                 value = [value,value(IndexMaster)];                
%             elseif(IndexSlave ~= -1 && IndexMaster ~= -1)
%                 value(IndexMaster) = value(IndexMaster) + value(IndexSlave); 
%                 %value(tmpIndex1) = value(tmpIndex2); 
%             end
% %             if IndexSlave ~= -1
% %                 %%delete indexSlave
% %                 value(IndexSlave)=[];
% %                 row(IndexSlave)=[];
% %                 column(IndexSlave)=[];
% %             end
%             
% %             if i==j
% %                 row = [row, IndexSlave]; 
% %                 column = [column, IndexSlave]; 
% %                 value = [value,1]; 
% %             end
%         end
%                 
%     end    
% end

% %add Boundary conditions from IBC
for i = 1:length(row)
        if (any(row(i) == IBC))
            row(i) = -1; 
            column(i) = -1;
            value(i) = NaN;
        end
        if (any(column(i) == IBC))
            row(i) = -1;
            column(i) = -1;
            value(i) = NaN;
        end
end

%get rid of the boundary condition rows and columns since they are 0
row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

%add back value as 1 at [IBC(i),IBC(i)]
for i = 1: length (IBC)
    rowTemp = IBC(i);
    fGlobal(rowTemp) = BCval(i); 
    row = [row,rowTemp];
    column = [column,rowTemp];
    value = [value,1];
end

KindexR = sparse(row, column, 1:length(row));
kGlobal = sparse(row, column, value);

end

function [nInt] = nIntergerPoints(p)
%define intergral rules
nInt = p*2;
end

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


function [Index] = IENindex(rowTemp, columnTemp, row, column)
%from a given row column combination, find the correspond index, store the
%index as another array
%KIndex(rowTemp,columnTemp) = Index;
if isempty(row) ==1
    Index = -1;
else
    indexRow = row == rowTemp;
    indexColumn = column == columnTemp;
    Index = find(indexRow&indexColumn);
    
    %KIndex(rowTemp,columnTemp) = Index;
    
    if isempty(Index)
        Index = -1;
    end
end
end




