function d = distance(X,Y,i,j)

diffx = (X(i)-X(j));

diffy = (Y(i)-Y(j));

d = sqrt(diffx^2 + diffy^2);
