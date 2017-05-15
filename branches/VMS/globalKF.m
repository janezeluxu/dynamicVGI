function [kGlobal, fGlobal] = globalKF(p,IENstruct,IBC,BCval,...
                                vertexData,M,N,iper)

global TotalDOF; 

% Assemble of element level matrix   

%create uniform p shape function tables, p is a list
ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    n = nIntegerPoints(p(j));
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
keyCount = 0;
for ele = 1:Grid_size
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions    
    [IENall,pAll] = elementIEN(ele,IENstruct);
    nssl = length(IENall);
    n = nIntegerPoints(max(pAll));
    
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
    simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,ShapeFunc,divShapeFunc);

    [elementK,elementF] = simplex.eleStiffMatrix();
    
    %sparse matrix input row, column, and value list, using row and column
    %combination as key
    for i = 1:nssl
        rowTemp = IENall(i);
        fGlobal(rowTemp) = fGlobal(rowTemp) + elementF(i);
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
                value(indexKglobal) = value(indexKglobal) + elementK(i,j);
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
            if(tmpIndex>0)
                fGlobal(j) = fGlobal(j) - BCval(i)*value(tmpIndex); 
            end
        end       
    end
end

% eliminate entries with BC's in stiffness matrix 
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
row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

% account for BC in stiffness matrix diagonal and forcing vector
for i = 1: length (IBC)
    rowTemp = IBC(i);
    fGlobal(rowTemp) = BCval(i); 
    
    % Set the BC's contribution to the diagonal to 1
    row = [row,rowTemp];
    column = [column,rowTemp];
    value = [value,1];
end

% Account for periodicity 
for i = 1:length(iper)
    if(iper(i) ~= i)
        rowTemp = iper(i);

        % forcing vector
        fGlobal(i) = fGlobal(i) + fGlobal(rowTemp);
        fGlobal(rowTemp) = fGlobal(i); 
        
        % stiffness matrix
        for j = 1:TotalDOF
            tmpIndex1 = IENindex(i,j,row,column);
            tmpIndex2 = IENindex(rowTemp,j,row,column);
            if(tmpIndex2 == -1 && tmpIndex1~=-1) 
                row = [row, rowTemp]; 
                column = [column, j]; 
                value = [value,value(tmpIndex1)]; 
            elseif(tmpIndex1 == -1 && tmpIndex2 ~= -1)
                row = [row, i]; 
                column = [column, j]; 
                value = [value,value(tmpIndex2)];                
            elseif(tmpIndex1 ~= -1 && tmpIndex2 ~= -1)
                value(tmpIndex2) = value(tmpIndex2) + value(tmpIndex1); 
                value(tmpIndex1) = value(tmpIndex2); 
            end
        end
    end    
end

% XXX for periodicity, still need to remove rows and then copy back to
% solution later. When we do this, we can remove one of the if clauses in
% the section above. XXX

% build sparse stiffness matrix
kGlobal = sparse(row, column, value);

end

function [Index] = IENindex(rowTemp, columnTemp, row, column)
%from a given row column combination, find the correspond index
if isempty(row) ==1
    Index = -1;
else
    indexRow = row == rowTemp;
    indexColumn = column == columnTemp;
    Index = find(indexRow&indexColumn);
    if isempty(Index)
        Index = -1;
    end
end
end




