function drg_pcolor(X,Y,C)
%pcolor does not plot the last row/column

szC=size(C);

X_extended=zeros(szC(1)+1,szC(2)+1);
X_extended(1:end-1,1:end-1)=X;
X_delta=mean(X(1,2:end)-X(1,1:end-1));
X_extended(:,end)=X(1,end)+X_delta;

Y_extended=zeros(szC(1)+1,szC(2)+1);
Y_extended(1:end-1,1:end-1)=Y;
Y_delta=mean(Y(2:end,1)-Y(1:end-1,1));
Y_extended(end,:)=Y(end,1)+Y_delta;

C_extended=zeros(szC(1)+1,szC(2)+1);
C_extended(1:end-1,1:end-1)=C;
pcolor(X,Y,C)