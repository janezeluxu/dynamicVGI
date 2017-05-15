function [gradu] = reconstrDerivs(Grid_size,ShapeFuncTable,divSFtable,vertexData,IENstruct,u0)

global TotalDOF; 

gradu = zeros(TotalDOF,2); 
rmass = zeros(TotalDOF,1);

for ele = 1:Grid_size
    % Initialize local variables
    rmassl = 0;
    gradul = zeros(1,2); 
    
    % Get the element IEN and orders
    [IENall,pAll] = elementIEN(ele,IENstruct);
    
    % Get shape functions
    n = nIntegerPoints(max(pAll));
    if range(pAll) == 0     % uniform p 
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
    else
        qPoints = TriGaussPoints(n);    % variable p
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, divShapeFunc ] = simplexsf.variablePShapeFunctionTable();    
    end
      
    %find element Shape Functions
    simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,vertexData,ShapeFunc,divShapeFunc,u0,gradu);
    
    % Get the quadrature points
    quadraturePoints = TriGaussPoints(n);            
    nP = size(quadraturePoints,1);
    sizeN = length(IENall);    
    Jw = zeros(nP,1);
    [JInverse, detJ,~,~] =  jacobian(simplex);
    
    % Build local rmass and gradu
    for k = 1:nP     % Loop over quad pts
        Jw(k) = detJ*quadraturePoints(k,3);   
        
        % Get ugrad at the quad pt
        ugradtmp = zeros(1,2); 
        for i = 1:sizeN
            gradNaGlobal = divShapeFunc(i,:,k)*JInverse;

            ugradtmp = ugradtmp + gradNaGlobal*u0(IENall(i));
        end
        for i = 1:sizeN       
            NaGlobal = ShapeFunc(i,k);             
            gradul = gradul + NaGlobal*ugradtmp*Jw(k);
            rmassl = rmassl + NaGlobal*Jw(k);
        end
    end
    
    % Assemble to global DOF    
    for i = 1:length(IENall)
        rmass(IENall(i)) = rmass(IENall(i)) + rmassl;
        gradu(IENall(i),:) = gradu(IENall(i),:) + gradul;
    end
end    

% Divide gradu by rmass to get grad(u) at nodes
gradu(:,1) = gradu(:,1)./(rmass+1e-16);
gradu(:,2) = gradu(:,2)./(rmass+1e-16);    

% Check for NaN's
for i = 1:1:TotalDOF
    for j = 1:2
        if(isnan(gradu(i,j)))
            gradu(i,j) = 0.0;
        end
    end
end

end