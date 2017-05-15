classdef SimplexStiffMatrix < handle
    %calculate element level matrix, stiffness matrix and force matrix
    
    properties
        %element ID number
        EleID;
        %IEN structure including vertexIEN, edgeIEN, faceIEN
        IENstructure;
        %element level IEN arra
        IENs;
        %number of intergral points;
        nInt;
        %vertex coordinates
        VertexData;
        %element shape functions
        EleShapeFunc;
        %element shape function derivitives
        EledivShapeFunc;
        %predicted solution
        u;
        %gradient of solution at nodes
        ugrad;
        
    end
    
    methods (Access = public)
        
        function simplex = SimplexStiffMatrix(elementID,ienstruct,IENall,nInterger, vertexCord,ShapeFunction,divSF)
            %contruction
            simplex.EleID = elementID;
            simplex.IENs = IENall;
            simplex.IENstructure = ienstruct;
            simplex.nInt = nInterger;
            simplex.VertexData = vertexCord;
            simplex.EleShapeFunc = ShapeFunction;
            simplex.EledivShapeFunc = divSF;
        end
        
        function [elementK,elementF] = eleStiffMatrix(simplex)
            global stabflag; 
            global itau; 
            global kappa;
            global a; 
            global force; 
            global nsd; 
            global h; 
            
            %element stiffness matrix and force
            IENall = simplex.IENs;
            ShapeFunc = simplex.EleShapeFunc;
            divShapeFunc = simplex.EledivShapeFunc;            
            n = simplex.nInt;
            
            quadraturePoints = TriGaussPoints(n);            
            nP = size(quadraturePoints,1);
            sizeN = length(IENall);
            elementK = zeros(sizeN, sizeN);
            elementF = zeros(sizeN,1);
            Jw = zeros(nP,1);

            % Get element size h = sqrt(A)= sqrt(int(1))
            [JInverse, detJ,gij,~,~] =  jacobian(simplex);
            
            % Loop over quad points            
            for k = 1:nP     
                Jw(k) = detJ*quadraturePoints(k,3);                               
                
                % build element stiffness matrix and force vector
                for i = 1:sizeN
                    NaGlobal = ShapeFunc(i,k); 
                    gradNaGlobal = divShapeFunc(i,:,k)*JInverse;
                    elementF(i) = elementF(i)+Jw(k)*ShapeFunc(i,k)*force;
                    
                    % Stability part of force vector
                    if(stabflag==1)
                        if(itau == 1)
                            alpha = norm(a)*h/(2*norm(kappa));
                            zi = coth(alpha)-1/alpha;           
                            if(nsd==1)
                                tau = h*zi/(2*a(1));   
                            else
                                tau = h*zi*norm(a)/norm(kappa);     
                            end       
                        elseif(itau == 2)
                            if(nsd==1)
                                tau1 = (h/(2*a(1)));
                                tau2 = h^2/(4*kappa(1,1)); 
                                tau = 1/sqrt(tau1^-2 + 9*tau2^-2);
                            else
                                tau1 = a*gij*a';
                                tau2 = 9*sum(dot(kappa*gij,kappa*gij)); 
                                tau = (tau1+tau2)^-.5;
                            end
                        end
                        elementF(i) = elementF(i) + Jw(k)*gradNaGlobal*a'*tau*force; 
                    end
                    
                    for j = 1:sizeN                        
                        NbGlobal = ShapeFunc(j,k); 
                        gradNbGlobal = divShapeFunc(j,:,k)*JInverse;
                        
                        % Diffusion
                        elementK(i,j) = elementK(i,j)+gradNaGlobal*kappa*gradNbGlobal'*Jw(k);

                        % Advection
                        elementK(i,j) = elementK(i,j)-gradNaGlobal*a'*NbGlobal*Jw(k);     

                        % Stab
                        if(stabflag ==1)
                            elementK(i,j) = elementK(i,j) + gradNaGlobal*a'*tau*gradNbGlobal*a'*Jw(k); 
                        end
                    end % end of j loop
                end % end of i loop    
            end %end of k loop
        end % end of fumction        
    end  
end

