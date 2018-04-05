function drg_pcolor(X,Y,C)
%pcolor does not plot the last row/column

szC=size(C);

X_extended=zeros(szC(1)+1,szC(2)+1);
X_extended(1:end-1,1:end-1)=X;
X_delta=X(1,2)-X(1,1);
x_bottom=[1:szC(2)+1]*X_delta;
X_extended(end,:)=x_bottom;
x_last_column=(X_delta+X(end,end))*ones(szC(1)+1,1);
X_extended(:,end)=x_last_column;

Y_extended=zeros(szC(1)+1,szC(2)+1);
Y_extended(1:end-1,1:end-1)=Y;
y_last_row=x_last_column';
Y_extended(end,:)=y_last_row;
y_last_column=x_bottom';
Y_extended(:,end)=y_last_column;

C_extended=zeros(szC(1)+1,szC(2)+1);
C_extended(1:end-1,1:end-1)=C;
pcolor(X_extended,Y_extended,C_extended)