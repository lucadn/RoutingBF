function [DOAVect, numNodes,affectedNodesIndexes]=evaluate_DOA_Vector(xVectorReal,yVectorReal,transmitterID,receiverID,maxDistance)
xVector=xVectorReal-xVectorReal(transmitterID);
yVector=yVectorReal-yVectorReal(transmitterID);
numNodes=0;
fprintf(1,'Transmitter ID: %d \n', transmitterID-1);
fprintf(1,'Receiver ID: %d \n', receiverID-1);
for i=1:length(xVector)
    d=distance(xVector,yVector,transmitterID,i);
    if d<maxDistance && i~=receiverID && i~=transmitterID
        d
        numNodes=numNodes+1;
        if(xVector(transmitterID)>xVector(i))
                        theta=180+atand((yVector(transmitterID)-yVector(i))/(xVector(transmitterID)-xVector(i)));
                    elseif (yVector(transmitterID)>yVector(i))
                        theta = 360+atand((yVector(transmitterID)-yVector(i))/(xVector(transmitterID)-xVector(i)));
                    else
                        theta=atand((yVector(i)-yVector(transmitterID))/(xVector(i)-xVector(transmitterID)));
                    end
        DOAVect(numNodes)=theta;
%         if(xVector(i)<xVector(j))
%                         DOAVect(numNodes) = atand((yVector(i)-yVector(j))/d);
%                     else
%                         theta=180+atand((yVector(i)-yVector(j))/(d));
%                     end
        affectedNodesIndexes(numNodes)=i;
        fprintf(1,'Node %d normxCoord = %f at distance %f with angle %f\n', i-1, xVector(i),d, DOAVect(numNodes));
    end
end
fprintf(1,'out\n');
fprintf(1,'out\n');

