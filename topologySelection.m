function [x,y]= topologySelection(totalNodes, DIM_MAX, R, Spacing, choice)
%int i,j;
% Circle = 1
% Random placement = 2
% Grid placement = 3
% Load from file = 4
x=zeros(1,totalNodes);
y=zeros(1,totalNodes);
MIN_DISTANCE=10;
PIGR=3.1415;
%FILE *fdispo;
%float dist_x, dist_y, dist_quad;
switch(choice)
    
    case 1 %Circle
        for i=1:totalNodes
            x(i)=R*cos(2*PIGR*i/totalNodes)+R;
            y(i)=R*sin(2*PIGR*i/totalNodes)+R;
        end
    case 2 %Random
        for i=1:totalNodes
            x(i)=DIM_MAX*rand(1,1);
            y(i)=DIM_MAX*rand(1,1);
            for j=1:i
                if(sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2) < MIN_DISTANCE)
                    x(i)=DIM_MAX*rand(1,1);
                    y(i)=DIM_MAX*rand(1,1);
                    j=0;
                end
            end
        end
        %char fName[50];
        fName=sprintf('Topology.txt');
        fdispo=fopen(fName, 'w+');
        for i=1:totalNodes
            fprintf(fdispo, '%f %f\n', [x(i),y(i)]);
        end
        fclose(fdispo);
    case 3 %Grid
        N = round(sqrt(totalNodes));
        for i=1:totalNodes
            x(i)=Spacing*mod(i-1,N);
            y(i)=Spacing*floor((i-1)/N);
        end
        
    case 4 %From file
        
        %char fName[50];
        fName=sprintf('Topology.txt');
        fdispo=fopen(fName,'r');
        %for i=1:totalNodes
            matrix_temp=fscanf(fdispo, '%f %f',[2,inf]);
            x=matrix_temp(1,:);
            y=matrix_temp(2,:);
            
        %end
        fclose(fdispo);
end