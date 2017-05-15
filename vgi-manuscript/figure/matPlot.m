function [] = matPlot()

[x, y,u3Pe125, u2Pe125, u0Pe125,u3Pe12p5, u2Pe12p5, u0Pe12p5,...
    u3Pe1p25, u2Pe1p25, u0Pe1p25,u3Pep125, u2Pep125, u0Pep125,...
    TRu3Pe125, TRu2Pe125, TRu0Pe125,TRu3Pe12p5, TRu2Pe12p5, TRu0Pe12p5,...
    TRu3Pe1p25, TRu2Pe1p25, TRu0Pe1p25,TRu3Pep125, TRu2Pep125, TRu0Pep125,Trx,Try] = uValue();

[xDiag, yDiag,Pe125value3, Pe125value2, Pe125value0,Pe12p5value3, Pe12p5value2, Pe12p5value0,...
    Pe1p25value3, Pe1p25value2, Pe1p25value0,Pep125value3, Pep125value2, Pep125value0,...
    TRPe125value3, TRPe125value2, TRPe125value0,TRPe12p5value3, TRPe12p5value2, TRPe12p5value0,...
    TRPe1p25value3, TRPe1p25value2, TRPe1p25value0,TRPep125value3, TRPep125value2, TRPep125value0,...
    xTRdiag, yTRdiag] = DiagValue();

xDomain = [0,1];   % limit of domain in x
yDomain = [0,1]; 
count = 1; 
xu = unique(x); 
yu = unique(y);  
zz = zeros(length(xu),length(yu)); 

minX = min(xDomain); 
maxX = max(xDomain); 
minY = min(yDomain); 
maxY = max(yDomain); 
%colormap(jet)
figure(1)
 %% first 4 plots Pe125
subplot(4,4,1) 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u3Pe125(count); 
        count = count+1; 
    end
end

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,2) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u2Pe125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,3) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u0Pe125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])
subplot(4,4,4) 
 for i = 1:length(xDiag)
     s(i) = ((xDiag(i) - xDiag(1))^2 + (yDiag(i) - yDiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),Pe125value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),Pe125value2(1:2:l),'b*')
 hold on
 plot(s,Pe125value0,'k-') 
 ylim([-0.2 1.2])
 xlim([0 1.5])
 %% next 4 plots Pe12.5
 subplot(4,4,5) 
 count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u3Pe12p5(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,6) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u2Pe12p5(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,7) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u0Pe12p5(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,8) 
 for i = 1:length(xDiag)
     s(i) = ((xDiag(i) - xDiag(1))^2 + (yDiag(i) - yDiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),Pe12p5value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),Pe12p5value2(1:2:l),'b*')
 hold on
 plot(s,Pe12p5value0,'k-') 
 ylim([-0.2 1.2])
 xlim([0 1.5])
  %% 4 plots Pe1.25
 subplot(4,4,9) 
 count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u3Pe1p25(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,10) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u2Pep125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,11) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u0Pe1p25(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,12) 
 for i = 1:length(xDiag)
     s(i) = ((xDiag(i) - xDiag(1))^2 + (yDiag(i) - yDiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),Pe1p25value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),Pe1p25value2(1:2:l),'b*')
 hold on
 plot(s(1:2:l),Pe1p25value0(1:2:l),'k-') 
 ylim([-0.2 1.2])
 xlim([0 1.5]) 
   %% 4 plots Pe.125
 subplot(4,4,13) 
 count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u3Pep125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,14) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u2Pe1p25(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,15) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u0Pep125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,16) 
 for i = 1:length(xDiag)
     s(i) = ((xDiag(i) - xDiag(1))^2 + (yDiag(i) - yDiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),Pep125value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),Pep125value2(1:2:l),'b*')
 hold on
 plot(s(1:2:l),Pep125value0(1:2:l),'k-') 
 ylim([-0.2 1.2])
 xlim([0 1.5]) 
 %% set domain
xDomain = [0,1];   % limit of domain in x
yDomain = [0,1]; 
count = 1; 
x = Trx;
y = Try;
xu = unique(x); 
yu = unique(y);  
zz = zeros(length(xu),length(yu)); 

minX = min(xDomain); 
maxX = max(xDomain); 
minY = min(yDomain); 
maxY = max(yDomain); 
figure(2)
 %% first 4 plots Pe125
subplot(4,4,1) 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu3Pe125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,2) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu2Pe125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,3) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu0Pe125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])
subplot(4,4,4) 
 for i = 1:length(xTRdiag)
     s(i) = ((xTRdiag(i) - xTRdiag(1))^2 + (yTRdiag(i) - yTRdiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),TRPe125value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),TRPe125value2(1:2:l),'b*')
 hold on
 plot(s,TRPe125value0,'k-') 
 ylim([-0.2 1.2])
 xlim([0 1.5]) 
 %% next 4 plots Pe12.5
 subplot(4,4,5) 
 count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu3Pe12p5(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,6) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu2Pe12p5(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,7) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu0Pe12p5(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,8) 
 for i = 1:length(xTRdiag)
     s(i) = ((xTRdiag(i) - xTRdiag(1))^2 + (yTRdiag(i) - yTRdiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),TRPe12p5value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),TRPe12p5value2(1:2:l),'b*')
 hold on
 plot(s,TRPe12p5value0,'k-') 
 ylim([-0.2 1.2])
  xlim([0 1.5])
  %% 4 plots Pe1.25
 subplot(4,4,9) 
 count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu3Pe1p25(count); 
        count = count+1; 
    end
end
grid on
sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,10) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu2Pe1p25(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,11) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu0Pe1p25(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,12) 
 for i = 1:length(xTRdiag)
     s(i) = ((xTRdiag(i) - xTRdiag(1))^2 + (yTRdiag(i) - yTRdiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),TRPe1p25value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),TRPe1p25value2(1:2:l),'b*')
 hold on
 plot(s(1:2:l),TRPe1p25value0(1:2:l),'k-') 
 ylim([-0.2 1.2])
 xlim([0 1.5]) 
   %% 4 plots Pe.125
 subplot(4,4,13) 
 count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu3Pep125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,14) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu2Pep125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,15) 
count = 1; 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = TRu0Pep125(count); 
        count = count+1; 
    end
end
grid on

sp = surf(xu,yu,zz);
sp.EdgeColor = 'none';
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])

subplot(4,4,16) 
 for i = 1:length(xTRdiag)
     s(i) = ((xTRdiag(i) - xTRdiag(1))^2 + (yTRdiag(i) - yTRdiag(1))^2)^0.5;
 end
 l = length(s);
 plot(s(1:2:l),TRPep125value3(1:2:l),'ro')
 hold on
 plot(s(1:2:l),TRPep125value2(1:2:l),'b*')
 hold on
 plot(s(1:2:l),TRPep125value0(1:2:l),'k-') 
 ylim([-0.2 1.2])
  xlim([0 1.5])
end