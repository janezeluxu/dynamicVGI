function [c1,L1,L2,M1,M2] = update_c1(c1,L1,L2, M1,M2, L1new,L2new,M1new,M2new, dt,M,N,vertexData,IBC,updateType)
global node_Size
c1new = zeros(node_Size,1);

if strcmp(updateType, 'pathLinec1') == 1
    for node = 1:M*N 
            Lnew = [L1new(node);L2new(node)];
            Mnew = [M1new(node);M2new(node)];
            c1new(node) = Lnew\Mnew;
    
     if c1new(node) < 0
         c1new(node)=0;
         %cc = c1new(node)
         %c1new(node)=abs(c1new(node));
     end
    end
    
    %BC= [IBC{1},IBC{2},IBC{3},IBC{4}];
    BC= [IBC{1},IBC{3}];
    for i = 1:length(BC)
        c1new(BC(i)) = 0;
    end
    
    %upwind c1 outlet BCs
%     BC1 = IBC{1};
%     BC1 = BC1(1:end-1);
%     for i = 1:length(BC1)
%         c1new(BC1(i)) = c1new(BC1(i)+M+1);
%     end
%     
%     BC3 = IBC{3};
%     BC3 = BC3(1:end-1);
%     for i = 1:length(BC3)
%         c1new(BC3(i)) = c1new(BC3(i)+M+1);
%     end
    c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
    
elseif strcmp(updateType, 'pathTubec1') == 1
    for node = 1:M*N
        
        Lnew = [L1new(node);L2new(node)];
        Mnew = [M1new(node);M2new(node)];
        c1new(node) = Lnew\Mnew;
            
        SurroundNode = vertexData{node,4};
        nNode = length(SurroundNode);
        
        for micronode = 1:nNode 
            c1new = c1new+c1new(SurroundNode(micronode));

        end
        c1new = c1new/(nNode+1);
        
        if c1new(node) < 0
            c1new(node)=0;
        end
    end
    
    BC= [IBC{1},IBC{2},IBC{3},IBC{4}];
    for i = 1:length(BC)
        c1(BC(i)) = 0;
    end
    
    c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;


end