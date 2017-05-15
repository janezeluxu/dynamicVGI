classdef SquareShapeFunc < handle
%calculate shape functions for Bernstain polynomial in lambda space,
%simplex order can be uniform or variable/transition elements
properties (Access = private)
        %variable p element IEN
        elementIEN;        
        %uniform element Order
        Order;             
        %variable p element order list
        OrderAll;          
        %list of quadrature points
        QuadraturePoints;

    end

methods(Access = public)
    function sqaure = SquareShapeFunc(QPoints,varargin)
        sqaure.QuadraturePoints = QPoints;
        %if the simplex has uniform order, only need p information
        if length(varargin) ==1
            sqaure.Order = varargin{1};
            sqaure.elementIEN = nan;
            sqaure.OrderAll = nan;
        %if element is variable order, will need the element IEN array and order list   
        else
            sqaure.Order = nan;
            sqaure.elementIEN = varargin{1};
            sqaure.OrderAll = varargin{2};
        end
    end
    
    function [ShapeFuncTable, divShapeFuncTable] = uniformPShapeFunctionTable(square)
        %this function will call uniformPShapeFunction and give the shape
        %function table of a quadrature list 
        qPoints = square.QuadraturePoints;
        p = square.Order;
        nShapeFunc = 4;
        nQuadPoints = size(qPoints,1);
        ShapeFuncTable = zeros(nShapeFunc,nQuadPoints);
        divShapeFuncTable = zeros(nShapeFunc,2,nQuadPoints);
        for i = 1: nQuadPoints
            IPoint = qPoints(i,:);
            [ ShapeFunc, divShapeFunc ] = uniformPShapeFunction(square,IPoint);
            ShapeFuncTable(:,i) =  ShapeFunc;
            divShapeFuncTable(:,:,i) = divShapeFunc;
        end
        
    end
    
    function [ShapeFuncTable, divShapeFuncTable] = variablePShapeFunctionTable(square)
        %this function will call variablePShapeFunction and give the shape
        %function table of a quadrature list 
        qPoints = square.QuadraturePoints;
        nQuadPoints = size(qPoints,1);
        IENall = square.elementIEN;
        nSF = length(IENall);
        ShapeFuncTable = zeros(nSF,nQuadPoints);
        divShapeFuncTable = zeros(nSF,2,nQuadPoints);
        for i = 1:nQuadPoints
            IPoint = qPoints(i,:);
            [ ShapeFunc, divShapeFunc ] = variablePShapeFunction(square,IPoint);
            ShapeFuncTable(:,i) =  ShapeFunc;
            divShapeFuncTable(:,:,i) = divShapeFunc;
        end
    end
    
    function [ ShapeFunc, divShapeFunc ] = uniformPShapeFunction(square,IntPoint)
        %calculate the uniform order shape functions at one given intergral point, this
        %function calls ShapeFunction
        p = square.Order;
        
        [ ShapeFunc, divShapeFunc ] = square.ShapeFunction(IntPoint);
    end
    
    function [ ShapeFunc, divShapeFunc ] = variablePShapeFunction(square,IntPoint)
        %calculate the variable order shape functions at one given intergral point, this
        %function calls ShapeFunction
        IENall = square.elementIEN;
        pAll = square.OrderAll;
        pE1 = pAll(1);
        pE2 = pAll(2);
        pE3 = pAll(3);
        pF = pAll(4);
        
        nssl = length(IENall);
        ShapeFunc = zeros(nssl,1);
        divShapeFunc = zeros(nssl,2);
        
        if pE1>1
            edgeAlpha1 = square.getLexico(pE1,'edge',1);
            [ SFedge1, divEdge1 ] = square.ShapeFunction(edgeAlpha1,IntPoint);
        else
            SFedge1 = [];
            divEdge1 = [];
        end
        if pE2>1
            edgeAlpha2 = square.getLexico(pE2,'edge',2);
            [ SFedge2, divEdge2 ] = square.ShapeFunction(edgeAlpha2,IntPoint);
        else
            SFedge2 = [];
            divEdge2 = [];
        end
        if pE3>1
            edgeAlpha3 = square.getLexico(pE3,'edge',3);
            [ SFedge3, divEdge3 ] = square.ShapeFunction(edgeAlpha3,IntPoint);
        else
            SFedge3 = [];
            divEdge3 = [];
        end
        if pF>2
            faceAlpha = square.getLexico(pF,'face',1);
            [ SFFace, divFace ] = square.ShapeFunction(faceAlpha,IntPoint);
        else
            SFFace = [];
            divFace = [];
        end
        
        lam2 = IntPoint(1);
        lam3 = IntPoint(2);
        lam1 = 1-lam2-lam3;
        lam = [lam1,lam2,lam3];
        dlam = [-1,-1;1,0;0,1];
        
        %corrected vertex mode
        for j = 1:3
            sumE1 = 0; sumE2 = 0;sumE3 =0;sumF = 0;
            divSum11 = 0;divSum21=0;divSum31=0;divsumF1=0;
            divSum12 = 0;divSum22=0;divSum32=0;divsumF2=0;
            
            if pE1>1
                for i = 1:pE1-1
                    coeff = edgeAlpha1(i,j)/pE1;
                    sumE1 = sumE1+coeff*SFedge1(i);
                    divSum11 = divSum11+coeff*divEdge1(i,1);
                    divSum12 = divSum12+coeff*divEdge1(i,2);
                end
            end
            if pE2>1
                for i = 1:pE2-1
                    coeff = edgeAlpha2(i,j)/pE2;
                    sumE2 = sumE2+coeff*SFedge2(i);
                    divSum21 = divSum21+coeff*divEdge2(i,1);
                    divSum22 = divSum22+coeff*divEdge2(i,2);
                end
            end
            if pE3>1
                for i = 1:pE3-1
                    coeff = edgeAlpha3(i,j)/pE3;
                    sumE3 = sumE3+coeff*SFedge3(i);
                    divSum31 = divSum31+coeff*divEdge3(i,1);
                    divSum32 = divSum32+coeff*divEdge3(i,2);
                end
            end
            if pF>2
                numberFaceMode = (pF-1)*(pF-2)/2;
                for i = 1:numberFaceMode
                    coeff = faceAlpha(i,j)/pF;
                    sumF = sumF+coeff*SFFace(i);
                    divsumF1 = divsumF1 + coeff*divFace(i,1);
                    divsumF2 = divsumF2 + coeff*divFace(i,2);
                end
            end
            ShapeFunc(j) = lam(j) - sumE1 - sumE2 - sumE3 - sumF;
            divShapeFunc(j,1) = dlam(j,1)-divSum11-divSum21-divSum31-divsumF1;
            divShapeFunc(j,2) = dlam(j,2)-divSum12-divSum22-divSum32-divsumF2;
        end
        
        
        ShapeFunc(4:nssl) = [SFedge1;SFedge2;SFedge3;SFFace];
        divShapeFunc(4:nssl,:) = [divEdge1;divEdge2;divEdge3;divFace];
    end


end

methods(Static)
    function [ ShapeFunc, divShapeFunc ] = ShapeFunction(IntPoint)
        %give one shape function and derivative value at one intergral point
        %and one lexico order combination in lambda space
        xi = IntPoint(1);
        eta = IntPoint(2);
        
        ShapeFunc = [(1-xi)*(1-eta);xi*(1-eta);xi*eta;(1-xi)*eta];
        divShapeFunc = [-(1-eta),-(1-xi);1-eta,-xi;eta,xi;-eta,(1-xi)];
        
        %if (a1+a2+a3 == 1)
        %    divShapeFunc = [-1,-1;1,0;0,1];
        %end
        
    end
    function [lexicoOrder] = getLexico(p,type,edgeNumber)
        %get the lexico order list of a given order and type(vertex, edge, or face, 
        %if it is a edge, also need to give the edge numbering.)
        if strcmp(type, 'vertex') ==1
            lexicoOrder = [p,0,0;0,p,0;0,0,p];
        elseif strcmp(type, 'edge') ==1
            lexicoOrder = zeros(p-1,3);
            if edgeNumber ==1
                for i = 1:p-1
                    lexicoOrder(i,:) = [i,p-i,0];
                end
                
                
            elseif edgeNumber ==2
                for i = 1:p-1
                    lexicoOrder(i,:) = [0,i, p-i];
                end
                
            elseif edgeNumber ==3
                for i = 1:p-1
                    lexicoOrder(i,:) = [p-i,0,i];
                end
            end
        elseif strcmp(type, 'face') ==1
            numberFaceMode = (p-1)*(p-2)/2;
            lexicoOrder = zeros(numberFaceMode,3);
            count = 0;
            for i = 1:p-2
                for j = 1:p-1-i
                    count = count+1;
                    lexicoOrder(count,:) = [i,j, p-i-j];
                end
            end
        end
        
    end
      
end
end