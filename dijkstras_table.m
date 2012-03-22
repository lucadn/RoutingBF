function [d_table,index,xVector,yVector,source,destination] = dijkstras_table(nn,DIM_MAX,choice)
%samples = 1;
% nn = [100];
% DIM_MAX = 18;
for c_mm = 1:1:length(nn)
%for c_m = 1:1:samples
    [x,y] = topologySelection(nn, DIM_MAX, 1, 2, choice);
    location = x + 1i*y;
    location(1) = 0;
     location(nn(c_mm)) = DIM_MAX+1i*DIM_MAX;
     destination = DIM_MAX+1i*DIM_MAX;
for c1 = 1:1:length(location)
    for c2 = 1:1:length(location)
        if c2~=c1
        table(c1,c2) = abs(location(c1)-location(c2))^2;
        else
            table(c1,c2) = 1e1000;
        end
        
    end
end
[table1,d_table,index] = dijkstras(table,location,DIM_MAX);
%end
end
source = location(1);
xVector = real(location);
yVector = imag(location);
end
% scatter(real(location),imag(location))
% hold on
% for c1 = 2:1:length(rt)
%  plot([real(rt(c1)) real(rt(c1-1))],[imag(rt(c1)) imag(rt(c1-1))],'r','LineWidth',2)
% end