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
            %element stiffness matrix and force
            
            global stabflag; 
            global itau; 
            global kappa;
            global a; 
            global force; 
            global nsd; 
            global h;
            
            IENall = simplex.IENs;
            ShapeFunc = simplex.EleShapeFunc;
            divShapeFunc = simplex.EledivShapeFunc;            
            n = simplex.nInt;
            
            quadraturePoints = TriGaussPoints(n);            
            nP = size(quadraturePoints,1);
            sizeN = length(IENall);
            elementK = zeros(sizeN, sizeN);
            elementF = zeros(sizeN,1);
            
            %divShapeFunc;
            [JInverse, detJ,gij,xCord,yCord] =  jacobian(simplex);
            %h = sqrt(sum(detJ*quadraturePoints(:,3)));
            
            Jw = detJ*quadraturePoints(:,3);
            
            lam = [1-quadraturePoints(:,1)-quadraturePoints(:,2),quadraturePoints(:,1:2)];
            x = xCord*lam';
            y = yCord*lam';
            
            sourceTerm = simplex.MMS(x,y,force);
            %sourceTerm = kappa*simplex.MMS_Diffusion(x,y,force)+a*simplex.MMS_Adv(x,y,force);
            
            for k = 1:nP                
              elementF = elementF+ShapeFunc(:,k)*Jw(k)*sourceTerm(k);
                gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
                
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
                                tau1 = h/2*a(1);
                                tau2 = h^2/(4*kappa(1,1)); 
                                tau = 1/sqrt(tau1^-2 + 9*tau2^-2);

                            else
                                tau1 = a*gij*a';
                                tau2 = 9*sum(dot(kappa*gij,kappa*gij)); 
                                tau = (tau1+tau2)^-.5;
                            end
                        end
                        elementF = elementF + Jw(k)*gradNaGlobal*a'*tau*sourceTerm(k); 
                    end
                    
                NbGlobal = ShapeFunc(:,k);
                
                %diffusion
                elementK = elementK + gradNaGlobal*kappa*gradNaGlobal'*Jw(k);
                % Advection
                elementK = elementK-gradNaGlobal*a'*NbGlobal'*Jw(k);
                
                % Stabilizer                
                if(stabflag ==1)
                    elementK = elementK + (gradNaGlobal*a')*tau*(gradNaGlobal*a')'*Jw(k);
                end
            end
        end

        
        function [ JInverse, detJ,gij,x,y] = jacobian(simplex)
            %calculate jacobian inverse and det J in lambda space, can put
            %this function into private access
            global nsd; 
            elementID = simplex.EleID;
            IENstruct = simplex.IENstructure;
            vertexData = simplex.VertexData;
            
            [~, vIDs,~] = IENstruct(elementID,:).vertexIEN;
            x1 = vertexData{vIDs(1),2}(1);
            y1 = vertexData{vIDs(1),2}(2);
            x2 = vertexData{vIDs(2),2}(1);
            y2 = vertexData{vIDs(2),2}(2);
            x3 = vertexData{vIDs(3),2}(1);
            y3 = vertexData{vIDs(3),2}(2);
            
            
            x = [x1,x2,x3];
            y = [y1,y2,y3];
            A = 0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1));
            
            JInverse = (0.5/A)*[(y3-y1),(x1-x3);(y1-y2),(x2-x1);];
            
            gradx = [(-x1+x2), (-x1+x3);(-y1+y2), (-y1+y3)];
            detJ = 0.5*(gradx(1,1)*gradx(2,2)-gradx(1,2)*gradx(2,1));
            
            gij = zeros(nsd,nsd); 
            dxi1(1) = 0.5/A*(y2-y3);
            dxi1(2) = 0.5/A*(x3-x2);
            dxi2(1) = 0.5/A*(y3-y1);
            dxi2(2) = 0.5/A*(x1-x3);
            
            for i = 1:nsd
                for j = 1:nsd
                    gij(i,j) = dxi1(i)*dxi1(j) + dxi2(i)*dxi2(j);
                end
            end

        end
    end
    methods(Static)
        
        function [Source] = MMS(x,y,q)
            global kappa;
            global a;
            if q ==0
                Source = zeros(length(x));
            elseif q == 1
                Source = 1*ones(length(x));                
            elseif q == 2
                Source_diff = kappa(1,1)*(-2*(x.^2-1))+kappa(2,2)*(-2*(y.^2-1));
                Source_adv = a(1)*2*x.*(y.^2-1)+a(2)*2*y.*(x.^2-1);
                Source = Source_diff+Source_adv;
            elseif q ==3
                Source = 12*x.^2;
            elseif q ==4
                Source = 2*sin(x).*sin(y);
            elseif q == 5
                Source = zeros(1,length(x));
                for i = 1:length(x)
                        if x(i)<= -0.6                           
                            Source = 0*ones(length(x),1);
                        elseif x(i) >0.6
                            Source = 0*ones(length(x),1);
                        else
                            Source = (5/3)*ones(length(x),1);
                        end
                end
                
                elseif q == 6
                Source = zeros(1,length(x));
                for i = 1:length(x)
                        if x(i)<= -0.6
                            Source = 0*ones(length(x),1);
                        elseif x(i) >0.6
                            Source = 0*ones(length(x),1);
                        else
                            a = 1/(4*0.6^3);
                            Source(i) = 12*a*x(i)^2;
                        end
                end
            elseif q == 7
                Source = zeros(1,length(x));
                for i = 1:length(x)
                    if x(i)<= -0.6
                        Source(i) = 2*(x(i)+1);
                    elseif x(i) >0.6
                        Source(i) = -2*(x(i)-1);
                    else
                        Source(i) = 2*((-5/6)*x(i)^2+0.7)+(5/3)*(1-y(i)^2);
                    end
                end
            
            elseif q == 8
                Source = zeros(1,length(x));
                for i = 1:length(x)
                    if (x(i)> -0.6 && x(i)< 0.6 && y(i)> -0.6 && y(i)< 0.6)
                        Source(i) = -2*(x(i)^2+y(i)^2-0.72);
                    else
                        Source(i) = 0;
                    end
                end
            end
        end

        
    end
       
end

