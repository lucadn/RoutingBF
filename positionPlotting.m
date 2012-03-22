function positionPlotting(X,Y,ref_index,funflag,xmin,xmax,ymin,ymax,maxDistance);
S=40;
C='b';
zoom=1;
% xmin=min(0,zoom*min(X));
% xmax=zoom*max(X);
% ymin=min(0,zoom*min(Y));
% ymax=zoom*max(Y);
xrange=xmax-xmin;
yrange=ymax-ymin;
scatter(X,Y,S,C,'filled');
Xord=sort(X);
Yord=sort(Y);
hold on
for i=1:length(X)
    switch(funflag)
        case 0
            if i~=ref_index
                rang=sprintf('RANG_{N%d}(N%d)=%g',ref_index,i,sqrt((X(i)-X(ref_index))^2+(Y(i)-Y(ref_index))^2));
                delta=(Y(i)-Y(ref_index))/Y(ref_index);
                text((X(i)+X(ref_index))/2,(Y(i)+Y(ref_index))/2+delta,rang);
                line([X(i) X(ref_index)], [Y(i) Y(ref_index)]);
            end
        case 1
            if i~=ref_index
                pos=sprintf('POS_{N%d}(N%d)=(%g,%g)',ref_index,i,X(i)-X(ref_index),Y(i)-Y(ref_index));
                delta=(Y(i)-Y(ref_index))/Y(ref_index);
                text((X(i)+X(ref_index))/2,(Y(i)+Y(ref_index))/2+delta,pos);
                line([X(i) X(ref_index)], [Y(i) Y(ref_index)]);
            end
        case 2
            pos=sprintf('POS(N%d)=(%g,%g)', i,X(i)-X(ref_index),Y(i)-Y(ref_index));
            delta=(Y(i)-Y(ref_index))/Y(ref_index);
            text((X(i)+X(ref_index))/2,(Y(i)+Y(ref_index))/2+delta,pos);
            if i~=ref_index
                line([X(i) X(ref_index)], [Y(i) Y(ref_index)]);
            end
        case 3
            pos=sprintf('(%g, %g)', X(i)+1,Y(i));
            delta=(Y(i)-Y(ref_index))/Y(ref_index);
            circle([X(i),Y(i)],maxDistance,200,'-r');
            text((X(i)-10),(Y(i)-3),pos);

    end
    nodeid=sprintf('N%d',i-1);
    text(X(i)+0.5,Y(i),nodeid);
end
axis([xmin xmax ymin ymax]);
%set(gca,'XTick', round(Xord));
%set(gca,'YTick', round(Yord));
xlabel('x','Position', [xmax (ymin-ymax/yrange)], 'FontSize', 14 );
ylabel('y','Position', [(xmin-xmax/xrange) ymax], 'FontSize', 14 , 'Rotation', 0.0);