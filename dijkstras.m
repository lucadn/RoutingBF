function [table,rt,index] = dijkstras(table,location,DIM_MAX)
P(1) = location(1);
n = 2;
nh = 2;
rt(1) = location(1);
table2 = table;
for c1 = 1:1:length(location)
    for c2 = 1:1:length(location)
        if c2~=c1
        rtt(c1,c2) = location(c2);
        else
            rtt(c1,c2) = 1e1000;
        end       
    end
end
while length(P) ~= length(location)
    temp = table(1,:);
    test = 0;
    while test == 0
        test = 1;
    [rows,cols] = find(temp==min(temp));
    next_hop = location(cols(1));
    for c4 = 1:1:length(P)
        if P(c4) == next_hop
           test = 0;
           temp(cols(1)) = 1e1000;
        end
    end
    if test == 1
        P(n) = next_hop;
    end
    end
    for c10 = 1:1:length(location)
    for c3 = 1:1:length(location)
        
        table2 = table;
        if length(find(P==location(c3)))==0 && c10 ~= c3   
            if table2(c10,c3) > table(cols(1),c3)+table(c10,cols(1))
            rtt(c10,c3) = P(n);
        end
        table(c10,c3) = min(table(c10,c3),table(cols(1),c3)+table(c10,cols(1)));

        end
        end
    end
    
    n = n+1;
end
rt(1) = DIM_MAX+1i*DIM_MAX;
n = 2;
y = length(location);
index(1) = y;
x = 1e1000;
source = 0+1i*0;
while x ~= source
y = find(location == rtt(1,y));
if x ~= location(y)
rt(n) = location(y);
index(n) = y;
x = location(y);
n = n+1;
else
x = source;
end
end
ph = length(rt);
rt(ph+1) = 0+1i*0;
index(ph+1) = 1;
rt = fliplr(rt);
index = fliplr(index);
end