
function [ShapeFuncTable,divSFtable] = ShapeTable(nInt,p)
ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    n = nIntergerPoints(p(j),nInt);
    
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
            
    Squaresf = SquareShapeFunc(qPoints,p(j));
    [ShapeFunc, divShapeFunc ] = Squaresf.uniformPShapeFunctionTable();    
    ShapeFuncTable{p(j)} =  ShapeFunc;
    divSFtable{p(j)} = divShapeFunc;
end
end