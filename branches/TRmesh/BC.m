classdef BC < handle
    
    properties (Access = private)
        ibc;
        bcval;
        %bctype;
        IPER;
        row;
        column;
        value;
        fglobal;
        %geomesh;
        edge;
        P;
        ienstruct;
        C1;
        VERTEX;
        N1;
        N2;
        N3;
        N4;
%         inlet;
%         ourlet;
    end
    
    methods(Access = public)

function BoundaryCondition = BC(ibc,bcval,Iper,Row,Column,Value,FGlobal,vertex,EdgeData,Order,IENSTRUCT,C,nn1,nn2,nn3,nn4)
    BoundaryCondition.ibc = ibc;
    BoundaryCondition.bcval = bcval;
    %BoundaryCondition.bctype = BCtype;
    BoundaryCondition.IPER = Iper;
    BoundaryCondition.row = Row;
    BoundaryCondition.column = Column;
    BoundaryCondition.value = Value;
    BoundaryCondition.fglobal = FGlobal;
    %BoundaryCondition.geomesh = geoMesh;
    BoundaryCondition.edge = EdgeData;
    BoundaryCondition.P = Order;
    BoundaryCondition.ienstruct = IENSTRUCT;
    BoundaryCondition.C1 = C;
    BoundaryCondition.VERTEX = vertex;
    BoundaryCondition.N1 = nn1;
    BoundaryCondition.N2 = nn2;
    BoundaryCondition.N3 = nn3;
    BoundaryCondition.N4 = nn4;
end

function [row,column,value,fGlobal] = ApplyBC(BoundaryCondition)
    global BCType
    row = BoundaryCondition.row;
    column = BoundaryCondition.column;
    value = BoundaryCondition.value;
    fGlobal = BoundaryCondition.fglobal;
    IBC = BoundaryCondition.ibc;
    BCval = BoundaryCondition.bcval;
    edgeData = BoundaryCondition.edge;
    
    n1 = BoundaryCondition.N1;
    n2 = BoundaryCondition.N2;
    n3 = BoundaryCondition.N3;
    n4 = BoundaryCondition.N4;
    
    IBC1 = IBC{1};
    IBC2 = IBC{2};
    IBC3 = IBC{3};
    IBC4 = IBC{4};
    
    BCval1 = BCval{1};
    BCval2 = BCval{2};
    BCval3 = BCval{3};
    BCval4 = BCval{4};
    
    edgeVertex = cell2mat(edgeData(:,2));
    [row,column,value,fGlobal] = BoundaryCondition.addBC(IBC1,BCval1,BCType{1},row,column,value,fGlobal,edgeVertex,n1);
    %kGlobal = sparse(row, column, value);
    %K1 = full(kGlobal)
    [row,column,value,fGlobal] = BoundaryCondition.addBC(IBC2,BCval2,BCType{2},row,column,value,fGlobal,edgeVertex,n2);
    %kGlobal = sparse(row, column, value);
    %K2 = full(kGlobal)
    [row,column,value,fGlobal] = BoundaryCondition.addBC(IBC3,BCval3,BCType{3},row,column,value,fGlobal,edgeVertex,n3);
    %kGlobal = sparse(row, column, value);
    %K3 = full(kGlobal)
    [row,column,value,fGlobal] = BoundaryCondition.addBC(IBC4,BCval4,BCType{4},row,column,value,fGlobal,edgeVertex,n4);
    %kGlobal = sparse(row, column, value);
    %K4 = full(kGlobal)
    %fGlobal
end
        
        function [row,column,value,fGlobal] = addBC(BoundaryCondition,IBC,BCval,BCType,row,column,value,fGlobal,edgeVertex,n)
            switch BCType
                case 'Strong'
                    [row,column,value,fGlobal] = BoundaryCondition.StrongBC(IBC,BCval,row,column,value,fGlobal);
                case {'PeriodicX'}
                    [row,column,value,fGlobal] = BoundaryCondition.PeriodicBC(row,column,value,fGlobal);
                case {'PeriodicY'}
                    [row,column,value,fGlobal] = BoundaryCondition.PeriodicBC(row,column,value,fGlobal);
                case {'Inlet'}
                    [InletElement,InletEdge] = BoundaryCondition.getElement(IBC,edgeVertex);
                    [row,column,value,fGlobal] = BoundaryCondition.Inlet(BCval,row,column,value,fGlobal,InletElement,InletEdge,n);
                case {'Outlet'}
                    [OutletElement,OutEdge] = BoundaryCondition.getElement(IBC,edgeVertex);
                    [row,column,value,fGlobal] = BoundaryCondition.Outlet(BCval,row,column,value,fGlobal,OutletElement,OutEdge,n);
                case {'flux'}
                    %[OutletElement,OutEdge] = BoundaryCondition.getElement(IBC,edgeVertex);
                    [row,column,value,fGlobal] = BoundaryCondition.flux(row,column,value,fGlobal);
                otherwise
                    disp('NO BC specified')
            end
            
        end
        
        
        function [row,column,value,fGlobal] = StrongBC(BoundaryCondition,IBC,BCval,row,column,value,fGlobal)
            %iper = BoundaryCondition.IPER;
            
            % Account for BC in forcing vector.
            global TotalDOF
            
            for i = 1:length(IBC)
                for j = 1:TotalDOF
                    % XXX use ien instead of ismember and TotalDOF, will be faster
                    if(~ismember(j,IBC) || BCval(i)==0)    % Don't do this if there's a BC here, it'll be 0 anyways
                        tmpIndex = BoundaryCondition.IENindex(j,IBC(i),row,column);
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
        
        
        function [row,column,value,fGlobal] = PeriodicBC(BoundaryCondition,row,column,value,fGlobal)
            iper = BoundaryCondition.IPER;
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
                        IndexSlave = BoundaryCondition.IENindex(i,j,row,column);
                        %IndexMaster = Kindex(rowTemp,j);
                        IndexMaster = BoundaryCondition.IENindex(rowTemp,j,row,column);
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

        function [row,column,value,fGlobal] = Inlet(BoundaryCondition,BCval,row,column,value,fGlobal,InletBC,InletEdge,nin)
            IENstruct = BoundaryCondition.ienstruct;
            p = BoundaryCondition.P;
            c1 = BoundaryCondition.C1;
            vertexData = BoundaryCondition.VERTEX;
            edgeD = BoundaryCondition.edge;           
            %% calculate added Stiffness matrix for inlet
            for ii = 1:length(InletBC)
                ele = InletBC(ii);
                BCvalEle = BCval(ii);
                [IENall,pAll] = elementIEN(ele,IENstruct);                
                vIDs = edgeD{InletEdge(ii),2};
                [~,edgeNumber] = find(vIDs(1) == IENall);
                
                n = nIntergerPoints(p,0);
                [xi, w] = GaussQuad(n, 1);
                qPoints = BC.getQuadrature(xi,w,edgeNumber);
                
                simplexsf = SimplexShapeFunc(qPoints,p);
                [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();
                
                nssl = length(IENall);
                n = nIntergerPoints(max(pAll),0);
                c1Ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
                simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,qPoints,vertexData,ShapeFunc,divShapeFunc,c1Ele);
                
                [elementK,elementF] = simplex.WeakBCInletStiffMatrix(nin,BCvalEle);
                %Assemble the inlet BC into global matrix
                for i = 1:nssl
                    rowTemp = IENall(i);
                    fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
                    for j = 1:nssl
                        columnTemp = IENall(j);
                        indexKglobal = BC.IENindex(rowTemp, columnTemp, row, column);
                        if indexKglobal == 0
                            %if indexKglobal == -1
                            %keyCount = keyCount+1; %its a new key
                            row = [row, rowTemp];
                            column = [column, columnTemp];
                            value = [value,elementK(i,j)];
                            %Kindex(rowTemp,columnTemp) = keyCount; %update Kindex
                        else
                            value(indexKglobal) = elementK(i,j)+value(indexKglobal);
                        end
                    end
                end
            end
            
        end
        
        
        function [row,column,value,fGlobal] = Outlet(BoundaryCondition,BCval,row,column,value,fGlobal,InletBC,OutEdge,nout)
            IENstruct = BoundaryCondition.ienstruct;
            p = BoundaryCondition.P;
            c1 = BoundaryCondition.C1;
            vertexData = BoundaryCondition.VERTEX;
            edgeD = BoundaryCondition.edge;
            
            %% calculate added Stiffness matrix for outlet
            for ii = 1:length(InletBC)
                ele = InletBC(ii);
                BCvalEle = BCval(ii);
                %BCvalEle = 0;
                [IENall,pAll] = elementIEN(ele,IENstruct);
                vIDs = edgeD{OutEdge(ii),2};
                [~,edgeNumber] = find(vIDs(1) == IENall);
                
                n = nIntergerPoints(p,0);
                [xi, w] = GaussQuad(n, 1);
                qPoints = BC.getQuadrature(xi,w,edgeNumber);
                simplexsf = SimplexShapeFunc(qPoints,p);
                [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();
                
                nssl = length(IENall);
                n = nIntergerPoints(max(pAll),0);
                
                c1Ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
                simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,qPoints,vertexData,ShapeFunc,divShapeFunc,c1Ele);
                
                [elementK,elementF] = simplex.WeakBCOutlettiffMatrix(nout,BCvalEle);
                
                %Assemble the outlet BC into global matrix
                for i = 1:nssl
                    rowTemp = IENall(i);
                    fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
                    for j = 1:nssl
                        columnTemp = IENall(j);
                        indexKglobal = BC.IENindex(rowTemp, columnTemp, row, column);
                        if indexKglobal == 0
                            row = [row, rowTemp];
                            column = [column, columnTemp];
                            value = [value,elementK(i,j)];
                        else
                            value(indexKglobal) = elementK(i,j)+value(indexKglobal);
                        end
                    end
                end
            end
            
            %kGlobal = sparse(row, column, value);
            %K = full(kGlobal)
            %fGlobal
            
        end
        
        function [row,column,value,fGlobal] = flux(BoundaryCondition,row,column,value,fGlobal)
            
            row = row;
            column = column;
            value = value;
            fGlobal = fGlobal;
        end
        
        function [Element,edgeList] = getElement(BoundaryCondition,IBC,edgeVertex)
            %edgeVertex = cell2mat(edgeData(:,2));
            Edge = BoundaryCondition.edge ;
            edgeList = geo_Mesh.edgeList(IBC, edgeVertex);
            Element =[];
            for i = 1:length(edgeList)
                Element = [Element,Edge{edgeList(i),3}];
            end
            
        end
        
    end
    methods(Static)
        function qPoints = getQuadrature(xi,w,edgeNumber)
            if (edgeNumber ==1)
                lxi = length(xi);
                qPoints = zeros(lxi,3);
                for i = 1:lxi
                        qPoints(i,:)=[xi(i),1-xi(i),w(i)];
                end
            elseif (edgeNumber ==2)
                lxi = length(xi);
                qPoints = zeros(lxi,3);
                for i = 1:lxi
                        qPoints(i,:)=[0,1-xi(i),w(i)];
                end
            elseif (edgeNumber ==3)
                lxi = length(xi);
                qPoints = zeros(lxi,3);
                for i = 1:lxi
                        qPoints(i,:)=[1-xi(i),0,w(i)];
                end
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
    end
end
