function [kGlobal, fGlobal] = globalKF(p,IENstruct,IBC,BCval,...
                                vertexData,iper,c1,EdgeData,n1,n2,n3,n4)
% Assemble of element level matrix   
global TotalDOF; 
global Grid_size;
%create uniform p shape function tables, p is a list
[ShapeFuncTable,divSFtable] = ShapeTable(0,p);

row = [];
column = [];
value = [];
fGlobal = zeros(TotalDOF,1);
%Kindex = ones(TotalDOF,TotalDOF)*-1;
Kindex = spalloc(TotalDOF,TotalDOF,TotalDOF*200);
keyCount = 0;
for ele = 1:Grid_size  
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions    
    [IENall,pAll] = elementIEN(ele,IENstruct);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll),0);
    
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

    
    c1Ele = (1/4)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3))+c1(IENall(4)));
    qPoints = 0;
    Square = SquareStiffMatrix(ele,IENstruct,IENall,n,qPoints,vertexData,ShapeFunc,divShapeFunc,c1Ele);
    %ShapeFunc
    %divShapeFunc
    %ele
    [elementK,elementF] = Square.eleStiffMatrix();
    
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
            if indexKglobal == 0
            %if indexKglobal == -1
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

%kGlobal = sparse(row, column, value);
%K = full(kGlobal)
%fGlobal

boundary = BC(IBC,BCval,iper,row,column,value,fGlobal,vertexData,EdgeData,p,IENstruct,c1,n1,n2,n3,n4);
[row,column,value,fGlobal] = boundary.ApplyBC();

%kGlobal = sparse(row, column, value);
%K = full(kGlobal)
%fGlobal
 %[row,column,value,fGlobal] = ApplyBC(row,column,value,fGlobal,Kindex,IBC,BCval,iper);
 kGlobal = sparse(row, column, value);
end


function [row,column,value,fGlobal] = ApplyBC(row,column,value,fGlobal,Kindex,IBC,BCval,iper)
global BCType
IBC1 = IBC{1};
IBC2 = IBC{2};
IBC3 = IBC{3};
IBC4 = IBC{4};

BCval1 = BCval{1};
BCval2 = BCval{2};
BCval3 = BCval{3};
BCval4 = BCval{4};

[row,column,value,fGlobal] = addBC(IBC1,BCval1,BCType{1},row,column,value,fGlobal,Kindex,iper);
 
[row,column,value,fGlobal] = addBC(IBC2,BCval2,BCType{2},row,column,value,fGlobal,Kindex,iper);
[row,column,value,fGlobal] = addBC(IBC3,BCval3,BCType{3},row,column,value,fGlobal,Kindex,iper);

[row,column,value,fGlobal] = addBC(IBC4,BCval4,BCType{4},row,column,value,fGlobal,Kindex,iper);

end

function [row,column,value,fGlobal] = addBC(IBC,BCval,BCType,row,column,value,fGlobal,Kindex,iper)
switch BCType
    case 'Strong' 
        [row,column,value,fGlobal] = StrongBoundaryCondition(row,column,value,fGlobal,IBC,BCval);
    case {'PeriodicX'}
        [row,column,value,fGlobal] = PeriodicBC(iper,row,column,value,fGlobal);
    case {'PeriodicY'}
        [row,column,value,fGlobal] = PeriodicBC(iper,row,column,value,fGlobal);
    case {'Inlet'} 
        [row,column,value,fGlobal] = InletBoundaryCondition(p,row,column,value,fGlobal,Kindex,IBC,BCval);
    case {'Outlet'}
        [row,column,value,fGlobal] = OutletBoundaryCondition(p,row,column,value,fGlobal,Kindex,IBC,BCval);
    otherwise
        disp('NO BC specified')
end

end


function [row,column,value,fGlobal] = StrongBoundaryCondition(row,column,value,fGlobal,IBC,BCval)
% Account for BC in forcing vector. 
global TotalDOF

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

%add Boundary conditions from IBC
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

%add back value as 1 at [IBC(i),IBC(i)]
for i = 1: length (IBC)
    rowTemp = IBC(i);
    fGlobal(rowTemp) = BCval(i);
    %fGlobal(rowTemp) = 0;
    row = [row,rowTemp];
    column = [column,rowTemp];
    value = [value,1];
end


row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

end


function [row,column,value,fGlobal] = PeriodicBC(iper,row,column,value,fGlobal)
global TotalDOF
iperBC = [];
% Account for periodicity 
for i = 1:length(iper)
    if(iper(i) ~= i)
        rowTemp = iper(i);
        iperBC = [iperBC,i];
        %rowTemp is primary/accumulator, i is non-primary
        % forcing vector
        fGlobal(rowTemp) = fGlobal(rowTemp) + fGlobal(i);
        %fGlobal(i) = fGlobal(rowTemp);         
        fGlobal(i) = 0; 
        
        % stiffness matrix
        for j = 1:TotalDOF
            %IndexSlave = Kindex(i,j);
            IndexSlave = IENindex(i,j,row,column);
            %IndexMaster = Kindex(rowTemp,j);
            IndexMaster = IENindex(rowTemp,j,row,column);
            columnPrim = iper(j);
            
            if(IndexMaster == 0 && IndexSlave~=0) 
                row = [row, rowTemp]; 
                column = [column, columnPrim]; 
                value = [value,value(IndexSlave)]; 
            elseif(IndexSlave ~= 0 && IndexMaster ~= 0)   
                value(IndexMaster) = value(IndexMaster) + value(IndexSlave); 
            end
        
        end   
    end    
end


%add Boundary conditions from iperBC
for i = 1:length(row)
        if  (any(row(i) == iperBC))  
            row(i) = -1; 
            column(i) = -1;
            value(i) = NaN;
        end
end
%fGlobal

for i =1:length(iperBC)
    row = [row,iperBC(i)];
    column = [column,iperBC(i)];
    value = [value,1];
    
    row = [row,iperBC(i)];
    column = [column,iper(iperBC(i))];
    value = [value,-1];
end

row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

end

function [row,column,value,fGlobal] = InletBoundaryCondition(p,row,column,value,fGlobal,Kindex,IBC,BCval)

global BCType
global TotalDOF
global kappa;
global a;

ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    %n = nIntergerPoints(p(j),nInt);
    qPoints = [0,0.5];
    simplexsf = SimplexShapeFunc(qPoints,p(j));
    [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();    
    ShapeFuncTable{p(j)} =  ShapeFunc;
    divSFtable{p(j)} = divShapeFunc;
end

%% find inlet element list based on IBC
InletBC = [];

%% calculate added Stiffness matrix for inlet
for ii = 1:length(InletBC)
    ele = InletBC(ii);
    [IENall,pAll] = elementIEN(ele,IENstruct);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll));
    
    %disp('uniform p element');
    ShapeFunc = ShapeFuncTable{pAll(1)};
    divShapeFunc = divSFtable{pAll(1)};
    
    c1Ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
    simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,ShapeFunc,divShapeFunc,c1Ele);
    %ShapeFunc
    %divShapeFunc
    
    [elementK,elementF] = simplex.WeakBCInletStiffMatrix();
     %Assemble the inlet BC into global matrix
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
            if indexKglobal == 0
            %if indexKglobal == -1
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
%% calculate added Stiffness matrix for outlet

%% Account for BC in forcing vector. 
for i = 1:length(IBC)
    for j = 1:TotalDOF
        % XXX use ien instead of ismember and TotalDOF, will be faster
        if(~ismember(j,IBC) || BCval(i)==0)    % Don't do this if there's a BC here, it'll be 0 anyways
            %tmpIndex = IENindex(j,IBC(i),row,column);
            tmpIndex = Kindex(j,IBC(i));
            if(tmpIndex>0)
                fGlobal(j) = fGlobal(j) - BCval(i)*value(tmpIndex); 
            end
        end       
    end
end

%kGlobal = sparse(row, column, value);
%K = full(kGlobal)
%fGlobal

if strcmp(BCType{4},'PeriodicX')==1 ||strcmp(BCType{2},'PeriodicY')==1
 [iperBC,row,column,value,fGlobal] = PeriodicBC(iper,row,column,value,fGlobal,Kindex);
 %disp ('BC is Periodic!!!');
else
    %disp ('BC is not Periodic!!!');
    iperBC = [];
end

%add Boundary conditions from IBC and iperBC
for i = 1:length(row)
        if (any(row(i) == IBC)) || (any(row(i) == iperBC))  
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

for i =1:length(iperBC)
    row = [row,iperBC(i)];
    column = [column,iperBC(i)];
    value = [value,1];
    
    row = [row,iperBC(i)];
    column = [column,iper(iperBC(i))];
    value = [value,-1];
end

%get rid of the boundary condition rows and columns since they are 0
row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

%add back value as 1 at [IBC(i),IBC(i)]
for i = 1: length (IBC)
    rowTemp = IBC(i);
    fGlobal(rowTemp) = BCval(i);
    %fGlobal(rowTemp) = 0;
    row = [row,rowTemp];
    column = [column,rowTemp];
    value = [value,1];
end

end

function [row,column,value,fGlobal] = OutletBoundaryCondition(p,row,column,value,fGlobal,Kindex,IBC,BCval)

global BCType
global TotalDOF
global kappa;
global a;

ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    %n = nIntergerPoints(p(j),nInt);
    qPoints = [0,0.5];
    simplexsf = SimplexShapeFunc(qPoints,p(j));
    [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();    
    ShapeFuncTable{p(j)} =  ShapeFunc;
    divSFtable{p(j)} = divShapeFunc;
end
InletBC = [];

%% calculate added Stiffness matrix for inlet
for ii = 1:length(InletBC)
    ele = InletBC(ii);
    [IENall,pAll] = elementIEN(ele,IENstruct);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll));
    
    %disp('uniform p element');
    ShapeFunc = ShapeFuncTable{pAll(1)};
    divShapeFunc = divSFtable{pAll(1)};
    
    c1Ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
    simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,ShapeFunc,divShapeFunc,c1Ele);
    %ShapeFunc
    %divShapeFunc
    
    [elementK,elementF] = simplex.WeakBCInletStiffMatrix();
     %Assemble the inlet BC into global matrix
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
            if indexKglobal == 0
            %if indexKglobal == -1
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
%% calculate added Stiffness matrix for outlet

%% Account for BC in forcing vector. 
for i = 1:length(IBC)
    for j = 1:TotalDOF
        % XXX use ien instead of ismember and TotalDOF, will be faster
        if(~ismember(j,IBC) || BCval(i)==0)    % Don't do this if there's a BC here, it'll be 0 anyways
            %tmpIndex = IENindex(j,IBC(i),row,column);
            tmpIndex = Kindex(j,IBC(i));
            if(tmpIndex>0)
                fGlobal(j) = fGlobal(j) - BCval(i)*value(tmpIndex); 
            end
        end       
    end
end

%kGlobal = sparse(row, column, value);
%K = full(kGlobal)
%fGlobal

if strcmp(BCType{4},'PeriodicX')==1 ||strcmp(BCType{2},'PeriodicY')==1
 [iperBC,row,column,value,fGlobal] = PeriodicBC(iper,row,column,value,fGlobal,Kindex);
 %disp ('BC is Periodic!!!');
else
    %disp ('BC is not Periodic!!!');
    iperBC = [];
end

%add Boundary conditions from IBC and iperBC
for i = 1:length(row)
        if (any(row(i) == IBC)) || (any(row(i) == iperBC))  
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

for i =1:length(iperBC)
    row = [row,iperBC(i)];
    column = [column,iperBC(i)];
    value = [value,1];
    
    row = [row,iperBC(i)];
    column = [column,iper(iperBC(i))];
    value = [value,-1];
end

%get rid of the boundary condition rows and columns since they are 0
row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

%add back value as 1 at [IBC(i),IBC(i)]
for i = 1: length (IBC)
    rowTemp = IBC(i);
    fGlobal(rowTemp) = BCval(i);
    %fGlobal(rowTemp) = 0;
    row = [row,rowTemp];
    column = [column,rowTemp];
    value = [value,1];
end

end

function [Index] = IENindex(rowTemp, columnTemp, row, column)
%from a given row column combination, find the correspond index, store the
%index as another array
%KIndex(rowTemp,columnTemp) = Index;
if isempty(row) ==1
    Index = 0;
else
    indexRow = row == rowTemp;
    indexColumn = column == columnTemp;
    Index = find(indexRow&indexColumn);
    
    %KIndex(rowTemp,columnTemp) = Index;
    
    if isempty(Index)
        Index = 0;
    else
        Index = Index(1);
    end
end
end

