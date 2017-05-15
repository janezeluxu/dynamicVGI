function [c1,L1,L2,M1,M2] = solvec1(p,M,N,vertexData,faceData,IBC,uHBCE,uHBCV,nInt,IENstruct,u)

global kappa; 
global a;
global h;
global node_Size;
global c2;
[ShapeFuncTable,divSFtable] = ShapeTable(nInt,p);

c1 = zeros(M*N,1);
L1 = zeros(M*N,1);
L2 = zeros(M*N,1);
M1 = zeros(M*N,1);
M2 = zeros(M*N,1);

for node = 2:node_Size-1 
    LS = LeastSquare(node,vertexData,u,uHBCE,uHBCV);
    
    a0 = LS(1);
    a1 = LS(2);
    a2 = LS(3);
      
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
        simplex = SimplexStiffMatrix(SurroundEle(microEle),IENstruct,IENall,nInt,0,vertexData,ShapeFuncTable,divSFtable,c1Ele);
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
        
        AA = AA+a*a'*h^(c2)*duh'*detJ;
        EE = EE+a'*uh*detJ;
        FF = FF + kappa*duh'*detJ;
        
    end
    uc = xyc*LS;
    sumA = sum(Area);
    uH = Area*uc;
    
    coordNode = vertexData{node,2};
    uH = sumA*[1,coordNode]*LS;
    
    He = 2*h;
    %He = sumA^0.5;
    BB = He^(c2)*a*a'*[a1;a2]*sumA;
    CC = -a'*uH;
    DD = kappa*[a1;a2]*sumA;
    

    L = AA-BB;
    M = CC+DD+EE-FF;
    
    L1(node) = L(1);
    L2(node) = L(2);
    M1(node) = M(1);
    M2(node) = M(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c1(node) = L\M;
     
     if c1(node) < 0
         c1(node)=0;
         %c1(node)=abs(c1(node));
     end
     
end

% BC= [IBC{1},IBC{2},IBC{3},IBC{4}];
% for i = 1:length(BC)
%     c1(BC(i)) = 0;
% end
end

function LS = LeastSquare(node,vertexData,u,uHBCE,uHBCV)
    global BCType;
    SurroundNode = vertexData{node,4}; %sorrounding node
    nNode = length(SurroundNode);
    CoordMatrix = zeros(nNode,3);
    uE = zeros(nNode,1);
    CoordMatrix(:,1) = 1;
    for micronode = 1:nNode        
        CoordMatrix(micronode,2:3)= vertexData{SurroundNode(micronode),2};
        uE(micronode) = u(SurroundNode(micronode));
    end
    LS = CoordMatrix\uE;
    %%edge bc
    if (strcmp(BCType{1},'Strong')==1||strcmp(BCType{1},'flux')==1)
        %disp ('Strong!!!');
        if any(node ==uHBCE(1,:)) ==1
            a1 = LS(2);
            a2 = LS(3);
            
            index = find(node ==uHBCE(1,:));
            Cord = vertexData{uHBCE(2,index),2};
            index = index(1);
            Value = u(uHBCE(2,index));
            a0 = Value-a1*Cord(1)-a2*Cord(2);
            LS  = [a0;a1;a2];
            %%vertexbc
%         elseif any(node ==uHBCV(1,:)) ==1
%             index = node ==uHBCV(1,:);
%             Cord = vertexData{uHBCV(2,index),2};
%             Value = u(uHBCV(2,index));
%             a1 = LS(2);
%             a2 = LS(3);
%             a0 = Value-a1*Cord(1)-a2*Cord(2);
%             LS  = [a0;a1;a2];
        end
    end
end