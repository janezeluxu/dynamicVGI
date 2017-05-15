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
    qPoints = TriGaussPoints(nInt);
    
    %ShapeFunc = ShapeFuncTable{pAll(1)};
    divShapeFunc = divSFtable{pAll(1)};
    %u1 u2 u3
    for j = 1:length(IENall)
        uE(j) = u(IENall(j));
        
    end
        
    [~, vIDs,~] = IENstruct(i,:).vertexIEN;
    x1 = vertexData{vIDs(1),2}(1);
    y1 = vertexData{vIDs(1),2}(2);
    x2 = vertexData{vIDs(2),2}(1);
    y2 = vertexData{vIDs(2),2}(2);
    x3 = vertexData{vIDs(3),2}(1);
    y3 = vertexData{vIDs(3),2}(2);
    gradx = [(-x1+x2), (-x1+x3);(-y1+y2), (-y1+y3)];
    A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
    JInverse = (0.5/A)*[(y3-y1),(x1-x3);(y1-y2),(x2-x1);];
    detJ = 0.5*(gradx(1,1)*gradx(2,2)-gradx(1,2)*gradx(2,1));
    
    %[JInverse, detJ,gij,xCord,yCord] =  jacobian(simplex);
    c1ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
    tau = c1ele*h;

    for k = 1:nInt 
        uH = divShapeFunc(:,:,k)'*uE';
        du = uH'*JInverse;        
        up = up+(tau*a*du')^2*detJ*qPoints(:,3);
    end
end

up = up^0.5;

end