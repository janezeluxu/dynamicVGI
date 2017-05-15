
function [] = testGalerkinApprox()
syms x y
xCord = [0,1,0]; yCord = [0,0,1];
x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);
a = x2+y2;

A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
lam1 = 1-lam2-lam3;

pAll = [2,3,4,4];
pF = pAll(4);
nssl = pAll(1)+pAll(2)+pAll(3)+(pF-1)*(pF-2)/2;
[ShapeFunc,~] = variableSFAnalytic(nssl,pAll, xCord,yCord);
%[ShapeFunc] = variableSFAnalytic(xCord,yCord)
%[ShapeFunc] = {lam1^2,lam2^2,lam3^2,2*lam1*lam2,2*lam2*lam3,2*lam1*lam3};
%[ShapeFunc] = {lam1^3,lam2^3,lam3^3,3*lam1^2*lam2,3*lam1*lam2^2,3*lam2^2*lam3,...
%    3*lam2*lam3^2, 3*lam1^2*lam3,3*lam1*lam3^2,6*lam1*lam2*lam3};

%  ShapeFunc{1}
%  ShapeFunc{2}
%  ShapeFunc{3}
%  ShapeFunc{4}
%  ShapeFunc{5}
%  ShapeFunc{6}
%  ShapeFunc{7}
%  %ShapeFunc{8}


for i = 1:length(ShapeFunc)
    fe = ShapeFunc{i}*(y^3+x^3+x^2*y+x*y^2);
    f(i) = int(int(fe,x,x1,a-y),y,y1,y3);
    for j = 1:length(ShapeFunc)
        Me = ShapeFunc{i}*ShapeFunc{j};
        M(i,j) =  int(int(Me,x,x1,a-y),y,y1,y3);
    end
end

d = M\f'
ans = d'*ShapeFunc;
simplify(ans)

%xA = [0,1,0,1/2,1/2,0];
%yA = [0,0,1,0,1/2,1/2];
xA = [0,1,0,1/3,2/3,2/3,1/3,0,0,1/3];
%yA = 
%dA = (xA.^2-1).*(yA.^2-1)
%dA = xA+yA.^2

end


function [ShapeFunc,divShapeFunc] = variableSFAnalytic(nssl,pAll, xCord,yCord)
%[lexicoOrder] = testLexiOrder(a,p);

x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);

syms x y 
A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
lam1 = 1-lam2-lam3;

pE1 = pAll(1);
pE2 = pAll(2);
pE3 = pAll(3);
pF = pAll(4);

%nssl = length(IENall);
ShapeFunc = cell(nssl,1);
divShapeFunc = cell(nssl,2);

if pE1>1
    edgeAlpha1 = getLexico(pE1,'edge',1);
    [ SFedge1, divEdge1 ] = SFAnalytical(edgeAlpha1,xCord,yCord);
   
else
    SFedge1 = [];
    divEdge1 = [];
end
if pE2>1
    edgeAlpha2 = getLexico(pE2,'edge',2);
    [ SFedge2, divEdge2 ] = SFAnalytical(edgeAlpha2,xCord,yCord);
   
else
    SFedge2 = [];
    divEdge2 = [];
end
if pE3>1
    edgeAlpha3 = getLexico(pE3,'edge',3);
    [ SFedge3, divEdge3 ] = SFAnalytical(edgeAlpha3,xCord,yCord);
else
    SFedge3 = [];
    divEdge3 = [];
end
if pF>2
    faceAlpha = getLexico(pF,'face',1);
    [ SFFace, divFace ] = SFAnalytical(faceAlpha,xCord,yCord);
else
    SFFace = [];
    divFace = [];
end


lam = [lam1,lam2,lam3];
dlam = [diff(lam1,x),diff(lam1,y);diff(lam2,x),diff(lam2,y);diff(lam3,x),diff(lam3,y)];

% %corrected vertex mode
for j = 1:3
    sumE1 = 0; sumE2 = 0;sumE3 =0;sumF = 0;
    divSum11 = 0;divSum21=0;divSum31=0;divsumF1=0;
    divSum12 = 0;divSum22=0;divSum32=0;divsumF2=0;
    if pE1>1
        for i = 1:pE1-1
            coeff = edgeAlpha1(i,j)/pE1;
            %SFedge1{1}
            sumE1 = sumE1+coeff*SFedge1{i};
            divSum11 = divSum11+coeff*divEdge1{i,1};
            divSum12 = divSum12+coeff*divEdge1{i,2};
        end
    end
    if pE2>1
        for i = 1:pE2-1
            coeff = edgeAlpha2(i,j)/pE2;
            sumE2 = sumE2+coeff*SFedge2{i};
            divSum21 = divSum21+coeff*divEdge2{i,1};
            divSum22 = divSum22+coeff*divEdge2{i,2};
        end
    end
    if pE3>1
        for i = 1:pE3-1
            coeff = edgeAlpha3(i,j)/pE3;
            sumE3 = sumE3+coeff*SFedge3{i};
            divSum31 = divSum31+coeff*divEdge3{i,1};
            divSum32 = divSum32+coeff*divEdge3{i,2};
        end
    end
    if pF>2
        numberFaceMode = (pF-1)*(pF-2)/2;
        for i = 1:numberFaceMode
            coeff = faceAlpha(i,j)/pF;
            sumF = sumF+coeff*SFFace{i};
            divsumF1 = divsumF1 + coeff*divFace{i,1};
            divsumF2 = divsumF2 + coeff*divFace{i,2};
        end
    end
    ShapeFunc{j} = lam(j) - sumE1 - sumE2 - sumE3 - sumF;
    divShapeFunc{j,1} = dlam(j,1)-divSum11-divSum21-divSum31-divsumF1;
    divShapeFunc{j,2} = dlam(j,2)-divSum12-divSum22-divSum32-divsumF2;
end

if nssl>3
ShapeFunc(4:nssl) = [SFedge1;SFedge2;SFedge3;SFFace];
divShapeFunc(4:nssl,:) = [divEdge1;divEdge2;divEdge3;divFace];
end
end

function [ShapeFunc,divShapeFunc] = SFAnalytical(lexicoOrder,xCord, yCord)

x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);

syms x y 
A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
lam1 = 1-lam2-lam3;


nsSF = size(lexicoOrder,1);
ShapeFunc = cell(nsSF,1);
divShapeFunc = cell(nsSF,2);

for i = 1: nsSF
    a1 = lexicoOrder(i,1);
    a2 = lexicoOrder(i,2);
    a3 = lexicoOrder(i,3);
    
    b1 = lam1.^a1;
    b2 = lam2.^a2;
    b3 = lam3.^a3;
    
    p = a1+a2+a3;
    coeff = factorial(p)/(factorial(a1)*factorial(a2)*factorial(a3));
    ShapeFunc{i} = coeff*b1*b2*b3;
    divShapeFunc{i,1} = diff(ShapeFunc(i),x);
    divShapeFunc{i,2} = diff(ShapeFunc(i),y);
end


end

function [lexicoOrder] = getLexico(p,type,edgeNumber)
        %get the lexico order list of a given order and type(vertex, edge, or face, 
        %if it is a edge, also need to give the edge numbering.)
        if strcmp(type, 'vertex') ==1;
            lexicoOrder = [p,0,0;0,p,0;0,0,p];
        elseif strcmp(type, 'edge') ==1;
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
      


function [ShapeFunc223] = variableAnalytic3D(xCord,yCord)
%[lexicoOrder] = testLexiOrder(a,p);

x1 = xCord(1); x2 = xCord(2); x3 = xCord(3); 
y1 = yCord(1); y2 = yCord(2); y3 = yCord(3);

syms x y 
A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
lam2 = (0.5/A)*(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y);
lam3 = (0.5/A)*(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y);
lam1 = 1-lam2-lam3;
%ShapeFunc
ShapeFunc223 = {lam1-lam1*lam2-2*lam1^2*lam3 - lam3^2*lam1-2*lam1*lam2*lam3,...
    lam2-lam2*lam3-lam1*lam2 - 2*lam1*lam2*lam3,...
    lam3-lam2*lam3-2*lam3^2*lam1 - lam1^2*lam3- 2*lam1*lam2*lam3,...
    2*lam1*lam2,2*lam2*lam3,3*lam3^2*lam1,3*lam3*lam1^2,6*lam1*lam2*lam3};

ShapeFunc223NoFace = {lam1-lam1*lam2-2*lam1^2*lam3 - lam3^2*lam1,...
    lam2-lam2*lam3-lam1*lam2 ,...
    lam3-lam2*lam3-2*lam3^2*lam1 - lam1^2*lam3,...
    2*lam1*lam2,2*lam2*lam3,3*lam3^2*lam1,3*lam3*lam1^2};

ShapeFunc224 = {lam1-lam1*lam2 - 3*lam3*lam1^3-3*lam3^2*lam1^2-lam3^3*lam1,...
    lam2-lam2*lam3-lam1*lam2,...
    lam3-lam2*lam3-lam3*lam1^3-3*lam3^2*lam1^2-3*lam3^3*lam1,...     
    2*lam1*lam2,2*lam2*lam3,4*lam3*lam1^3,6*lam3^2*lam1^2,4*lam3^3*lam1};
end


