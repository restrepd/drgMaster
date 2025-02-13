function drg_pcolor_new(X,Y,C)
%pcolor does not plot the last row/column

szC=size(C);

X_extended=zeros(szC(1)+1,szC(2)+1);
X_extended(1:end-1,1:end-1)=X;
X_delta=X(2,1)-X(1,1);
X_extended(end,:)=X(end,1)+X_delta;
X_extended(:,end)=X_extended(:,end-1);


Y_extended=zeros(szC(1)+1,szC(2)+1);
Y_extended(1:end-1,1:end-1)=Y;
Y_delta=Y(1,2)-Y(1,1);
Y_extended(:,end)=Y(1,end)+Y_delta;
Y_extended(end,:)=Y_extended(end-1,:);


C_extended=zeros(szC(1)+1,szC(2)+1);
C_extended(1:end-1,1:end-1)=C;
C_extended(end,1:end-1)=C(end,:);
C_extended(1:end-1,end)=C(:,end);
C_extended(end,end)=C(end,end);
 
pcolor(X_extended,Y_extended,C_extended)

pffft=1;
