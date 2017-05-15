function [c1] = solvec1(p,vertexData,faceData,IBC,nInt,IENstruct,u)

global kappa; 
global a;
global h;
global node_Size;
[ShapeFuncTable,divSFtable] = ShapeTable(nInt,p);

c1 = zeros(node_Size,1);

for node = 1:node_Size  
%for node = 1:2
    %node
    SurroundNode = vertexData{node,4}; %sorrounding node
    nNode = length(SurroundNode);
    CoordMatrix = zeros(nNode,3);
    uE = zeros(nNode,1);
    CoordMatrix(:,1) = 1;
    for micronode = 1:nNode        
        CoordMatrix(micronode,2:3)= vertexData{SurroundNode(micronode),2};
        uE(micronode) = u(SurroundNode(micronode));
    end
      A = CoordMatrix\uE;
      %a0 = A(1);
      a1 = A(2);
      a2 = A(3);
      
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);
    Area = zeros(1,nElement);
    xyc = zeros(nElement,3);
    xyc(:,1) = 1;
       
    AA = [0;0];
    EE = [0;0];
    FF = [0;0];
    for microEle = 1:nElement
        [IENall,pAll] = elementIEN(SurroundEle(microEle),IENstruct);
        c1Ele = 0;
        simplex = SimplexStiffMatrix(SurroundEle(microEle),IENstruct,IENall,nInt,vertexData,ShapeFuncTable,divSFtable,c1Ele);
        [JInverse, detJ,~,xCord,yCord] =  simplex.jacobian();
            
        quadraturePoints = TriGaussPoints(nInt);  
        lam = [1-quadraturePoints(:,1)-quadraturePoints(:,2),quadraturePoints(:,1:2)];
        x = xCord*lam';
        y = yCord*lam';
        xyc(microEle,2) = x;
        xyc(microEle,3) = y;
        
        Area(microEle) = faceData{SurroundEle(microEle),4};
        
        %----------------------------- uh ----------------%
        
        uA = zeros(1,length(IENall));
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
        end
        uh = uA*ShapeFunc;
        duh = uA*divShapeFunc*JInverse;
        
        he = (Area(microEle))^0.5;
        
        AA = AA+a*a'*h*duh'*detJ;
        EE = EE+a'*uh*detJ;
        FF = FF + kappa*duh'*detJ;
        
    end
    uc = xyc*A;
    sumA = sum(Area);
    
    He = 2*h;
    %He = sumA^0.5;
    BB = He*a*a'*[a1;a2]*sumA;
    CC = -a'*Area*uc;
    DD = kappa*[a1;a2]*sumA;
    
    L = AA-BB;
    M = CC+DD+EE-FF;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c1(node) = L\M;
    
     if c1(node) < 0
         c1(node)=0;
     end
end

%BC= [IBC{1},IBC{2},IBC{3},IBC{4}];
%for i = 1:length(BC)
%    c1(BC(i)) = 0;
%end
end