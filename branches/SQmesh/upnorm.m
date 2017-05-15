function up = upnorm(u,c1,vertexData,IENstruct,nInt,p)

%global kappa; 
%global force;
global a;
global Grid_size;
global h

[ShapeFuncTable,divSFtable] = ShapeTable(nInt,p);
up = 0;

for i = 1:Grid_size
    
    
    [IENall,pAll] = elementIEN(i,IENstruct);
    uE = zeros(1,length(IENall));      
    [xi, w] = GaussQuad(nInt,1);
    lxi = length(xi);
    qPoints = zeros(lxi*lxi,3);
    for l = 1:lxi
        for m = 1:lxi
            n = (l-1)*lxi+m;
            qPoints(n,:)=[xi(l),xi(m),w(l)*w(m)];
        end
    end
    nP = size(qPoints,1);
    %ShapeFunc = ShapeFuncTable{pAll(1)};
    divShapeFunc = divSFtable{pAll(1)};
    %u1 u2 u3 u4
    for j = 1:length(IENall)
        uE(j) = u(IENall(j));
        
    end
        
    [~, vIDs,~] = IENstruct(i,:).vertexIEN;
    x1 = vertexData{vIDs(1),2}(1);
    y1 = vertexData{vIDs(1),2}(2);
    x2 = vertexData{vIDs(2),2}(1);
    y4 = vertexData{vIDs(4),2}(2);
    
    J = [(x2-x1),0;0,(y4-y1)];
    JInverse = J^-1;
    detJ = det(J);
    
    %[JInverse, detJ,gij,xCord,yCord] =  jacobian(simplex);
    c1ele = (1/4)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3))+c1(IENall(4)));
    tau = c1ele*h;

    for k = 1:nP 
        uH = divShapeFunc(:,:,k)'*uE';
        du = uH'*JInverse;     
        %(tau*a*du')^2*detJ*qPoints(:,3)
        up = up+(tau*a*du')^2*detJ*qPoints(k,3);
    end
end

up = up^0.5;

end