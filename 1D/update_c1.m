function [c1,l1,m1] = update_c1(c1,l1, m1, l1new,m1new, dt,a,nPoints,updateType)
c1new = zeros(nPoints,1);

if strcmp(updateType, 'pathLinec1') == 1
    for node = 2:nPoints-1
        c1new(node) = (m1new(node)/l1new(node))/(a^2);
        if c1new(node)<0|| l1new(node) ==0
            c1new(node) = 0;
        end
        %c1new(1) = c1new(2);
    end
    %c1(1) = c1(2);
    %c1(nPoints) = c1(nPoints-1);
    c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
    
elseif strcmp(updateType, 'pathTubec1') == 1
    for node = 2:nPoints-1
            %l1new(node) = (l1new(node-1)+2*l1new(node)+l1new(node+1));
            %m1new(node) = (m1new(node-1)+2*m1new(node)+m1new(node+1));
            c1new(node) = (m1new(node)/l1new(node))/(a^2);
    end
    for node = 2:nPoints-1
        c1new(node) = (c1new(node-1)+2*c1new(node)+c1new(node+1))/4;
        %c1new(node) = (m1new(node)/l1new(node))/(a^2);
        if c1new(node)<0|| l1new(node) ==0
            c1new(node) = 0;
        end
    end
    c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
  
elseif strcmp(updateType, 'pathLineL1M1') == 1
        l1 = (1/(1+dt))*l1+(dt/(1+dt))*l1new;
        m1 = (1/(1+dt))*m1+(dt/(1+dt))*m1new;          
        for node = 2:nPoints-1
            c1(node) = ((m1(node))/(l1(node)))/(a^2);
            if c1(node)<0|| l1(node) ==0
                c1(node) = 0;
            end
        end
        
elseif strcmp(updateType, 'pathTubeL1M1') == 1
        l1 = (1/(1+dt))*l1+(dt/(1+dt))*l1new;
        m1 = (1/(1+dt))*m1+(dt/(1+dt))*m1new;          
        for node = 2:nPoints-1
            c1(node) = ((m1(node-1)+2*m1(node)+m1(node+1))/(l1(node-1)+2*l1(node)+l1(node+1)))/(a^2);
            if c1(node)<0|| l1(node) ==0
                c1(node) = 0;
            end
        end
end

end