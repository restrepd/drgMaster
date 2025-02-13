function drg_pcolor(X,Y,C)
%pcolor does not plot the last row/column

szC=size(C);

X_extended=zeros(szC(1)+1,szC(2)+1);
X_extended(1:end-1,1:end-1)=X;
X_delta=X(1,2)-X(1,1);
X_extended(1:end-1,end)=X(:,end)+X_delta;
X_extended(end,:)=X_extended(end-1,:);


Y_extended=zeros(szC(1)+1,szC(2)+1);
Y_extended(1:end-1,1:end-1)=Y;
Y_delta=Y(2,1)-Y(1,1);
Y_extended(end,:)=Y_extended(end-1,:)+Y_delta;
Y_extended(:,end)=Y_extended(:,end-1);

 
C_extended=zeros(szC(1)+1,szC(2)+1);
C_extended(1:end-1,1:end-1)=C;
C_delta_column=C(1,2)-C(1,1);
C_delta_row=C(2,1)-C(1,1);
C_extended(:,end)=[C(:,end); C(end,end)+C_delta_row];
C_extended(end,1:end-1)=C(end,:)+C_delta_column;

pcolor(X_extended,Y_extended,C_extended)
