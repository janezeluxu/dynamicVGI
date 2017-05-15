function [c1,L1,L2,M1,M2] = solvec1_new(p,M,N,vertexData,faceData,IBC,uHBCE,uHBCV,nInt,IENstruct,u)

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

for node = 1:node_Size 
    LS = LeastSquare(node,vertexData,u,uHBCE,uHBCV);
    
    a0 = LS(1);
    a1 = LS(2);
    a2 = LS(3);
      
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);
    Area = zeros(1,nElement);
       
    AA = [0;0];
    EE = [0;0];
    FF = [0;0];
    uH = 0;
    node
    xyc = zeros(nElement,3);
    for microEle = 1:nElement
        [IENall,pAll] = elementIEN(SurroundEle(microEle),IENstruct);
        c1Ele = 0;
        simplex = SquareStiffMatrix(SurroundEle(microEle),IENstruct,IENall,nInt,0,vertexData,ShapeFuncTable,divSFtable,c1Ele);
        [JInverse, detJ,~,xCord,yCord] =  simplex.jacobian();
        
        [xi, w] = GaussQuad(nInt,1); 
        lxi = length(xi);
            
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
                wxi(n) = w(l)*w(m);
            end
        end
        surrond = SurroundEle(microEle)
        xy
        microEle
        nP = size(xy,1);
        xyc(microEle,1) = 1;
        xyc(microEle,2) = xy(:,1);
        xyc(microEle,3) = xy(:,2);
        Area(microEle) = faceData{SurroundEle(microEle),4};
        %uH = uH+((xyc*LS)'*wxi')*Area(microEle);
        
        %----------------------------- uh ----------------%
        
        uA = zeros(1,length(IENall));
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
        end

        
%         [qPoints,nP] = getQuadrature(nInt);
%         uh = 0;
%         duh = 0;
%         for k = 1:nP
%             duh = duh+(divShapeFunc(:,:,k)'*uA')'*JInverse*detJ*qPoints(k,3);
%             uh = uh+ShapeFunc(:,k)'*uA'*detJ*qPoints(k,3);
%         end
        duh = uA*divShapeFunc*JInverse;
        uh = uA*ShapeFunc;
        
        he = (Area(microEle))^0.5;
        
        AA = AA+a*a'*h^(c2)*duh'*detJ;
        EE = EE+a'*uh*detJ;
        FF = FF + kappa*duh'*detJ;
        
    end
    %----------------------------- uH ----------------%
    uc = xyc*LS;
    sumA = sum(Area);
    
    He = 2*h;
    BB = He^(c2)*(a*a')*[a1;a2]*sumA;
    %CC = -a'*uH;
    CC = -a'*Area*uc
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
         c1(node) = 0;
         %c1(node)=abs(c1(node));
         
     end
end

% for i = 1:length(BC)
%     c1(BC(i)) = 0;
% end
end

function [qPoints,nP] = getQuadrature(n)
    [xi, w] = GaussQuad(n,1);
    lxi = length(xi);
    qPoints = zeros(lxi*lxi,3);
    for i = 1:lxi
        for m = 1:lxi
            n = (i-1)*lxi+m;
            qPoints(n,:)=[xi(i),xi(m),w(i)*w(m)];
        end
    end
    nP = size(qPoints,1);
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
    if (strcmp(BCType{1},'Strong')==1 ||strcmp(BCType{1},'flux')==1)
        %disp ('Strong!!!');
        if any(node ==uHBCE(1,:)) ==1
            a1 = LS(2);
            a2 = LS(3);
            index = node ==uHBCE(1,:);
            Cord = vertexData{uHBCE(2,index),2};
            Value = u(uHBCE(2,index));
            a0 = Value-a1*Cord(1)-a2*Cord(2);
            a0 = a0(1);
            LS  = [a0;a1;a2];
%             %%vertexbc
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
