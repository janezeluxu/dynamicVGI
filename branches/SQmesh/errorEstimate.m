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

[ShapeFuncTable,divSFtable] = ShapeTable(nInt,p);

%loop through all element
Error = 0;
for i = 1:Grid_size
    
    [IENall,pAll] = elementIEN(i,IENstruct);
    uE = zeros(1,length(IENall));
    n = nInt;       
    %qPoints = TriGaussPoints(n);
    
    [xi, w] = GaussQuad(n,1);
    %qPoints = [xi(1),xi(1),w(1)*w(1);xi(1),xi(2),w(1)*w(2);xi(2),xi(1),w(2)*w(1);xi(2),xi(2),w(2)*w(2)];
    lxi = length(xi);
    qPoints = zeros(lxi*lxi,3);
    for l = 1:lxi
        for m = 1:lxi
            n = (l-1)*lxi+m;
            qPoints(n,:)=[xi(l),xi(m),w(l)*w(m)];
        end
    end
    
    if range(pAll) == 0
        %It is a uniform element 
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
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

    simplex = SquareStiffMatrix(i,IENstruct,IENall,nInt,0,vertexData,ShapeFunc,divShapeFunc,0);
    [ JInverse, detJ,gij,xCord,yCord] = simplex.jacobian();

    x = zeros(lxi,1);
    y = zeros(lxi,1);
    xy =  zeros(lxi*lxi,2);
    for nQuad = 1:lxi
        x(nQuad) = (xCord(2)-xCord(1))*xi(nQuad);
        y(nQuad) = (yCord(4)-yCord(1))*xi(nQuad);
    end
    for l = 1:lxi
        for m = 1:lxi
            n = (l-1)*lxi+m;
            xy(n,:)=[x(l),y(m)];
        end
    end
   
    uA = Exact(xy(:,1),xy(:,2),force,a,kappa);
    Error = (((uH - uA')).^2*qPoints(:,3)*detJ)+Error;
    %if Error>1e-6
        %string = [ num2str(i),' Order is ', num2str(pAll), ' error is ',num2str(Error^0.5) ];
        %disp(string)
    %end
        
end
end





