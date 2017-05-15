function [Error,uH] = errorEstimate(u,vertexData,IENstruct,nInt,varargin)
        
global kappa; 
global force;
global a;
global Grid_size;
if length(varargin) ==1
    OrderList = varargin{1};
    p = OrderList;
elseif length(varargin) ==3
    PV = varargin{1};
    PH = varargin{2};
    PT = varargin{3};
    p = [PV,PH,PT];
end

[ShapeFuncTable] = ShapeTable(nInt,p);

%loop through all element
Error = 0;
for i = 1:Grid_size
    
    [IENall,pAll] = elementIEN(i,IENstruct);
    uE = zeros(1,length(IENall));
    n = nInt;       
    qPoints = TriGaussPoints(n);
        
    if range(pAll) == 0
        %It is a uniform element 
        ShapeFunc = ShapeFuncTable{pAll(1)};
        %intergral error over the element 
        for j = 1:length(IENall)
            uE(j) = u(IENall(j));
            
        end
       
        uH = uE*ShapeFunc;
        
        
    else
        string = [num2str(i),' is a transient element'];
        %disp(string)
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, ~ ] = simplexsf.variablePShapeFunctionTable();
        
        for j = 1:length(IENall)
            uE(j) = u(IENall(j));        
        end
        
        uH = uE*ShapeFunc;
    end
    
    [~, vIDs,~] = IENstruct(i,:).vertexIEN;
    x1 = vertexData{vIDs(1),2}(1);
    y1 = vertexData{vIDs(1),2}(2);
    x2 = vertexData{vIDs(2),2}(1);
    y2 = vertexData{vIDs(2),2}(2);
    x3 = vertexData{vIDs(3),2}(1);
    y3 = vertexData{vIDs(3),2}(2);
    xCord = [x1,x2,x3];
    yCord = [y1,y2,y3];
    gradx = [(-x1+x2), (-x1+x3);(-y1+y2), (-y1+y3)];
    detJ = 0.5*(gradx(1,1)*gradx(2,2)-gradx(1,2)*gradx(2,1));
    
    lam = [1-qPoints(:,1)-qPoints(:,2),qPoints(:,1:2)];
    x = xCord*lam';
    y = yCord*lam';
    
    uA = Exact(x,y,force,a,kappa);
    
    Error = (((uH - uA)).^2*qPoints(:,3)*detJ)+Error;
    %if Error>1e-6
        %string = [ num2str(i),' Order is ', num2str(pAll), ' error is ',num2str(Error^0.5) ];
        %disp(string)
    %end
        
end
end





