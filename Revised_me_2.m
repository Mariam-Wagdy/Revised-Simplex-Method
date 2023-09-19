clear
clc
A=[-3 -2;1 -1;9 7 ;3 7;-2 5];
C=[-1 8];
b= [-6 ;6; 108; 70; 35];
s=[-1;-1;-1;-1;-1];
p=1;
key=1;
D=zeros(1,2*length(A));
for x=1:length(s)
    switch s(x)
        case -1
            point=size(A);
            A(x,point(2)+1)=1;
            C(length(C)+1)=0;
        case 0
            point=size(A);
            A(x, point(2)+1)=1;
            D(point(2)+1)=1;
            p=p+1;
            C(length(C)+1)=0;
            key=2;
        case 1
            point=size(A);
            A(x, point(2)+1)=-1;
            A(x, point(2)+2)=1;
            D(point(2)+2)=1;
            p=p+1;
            C(length(C)+2)=0;
            key=2;
        otherwise
            disp(" wrong input");
    end
end
[m,n]=size(A); %% m: number of constraints (basic variables)
               %%n: number of vartiables
k=n-m;         %% number of nonbasic variables
o=1;
o2=1;
%%%%% getting basis%%%%%
for x=1:n
    if (nnz(A(:,x)) == 1) && (max( A(:,x)) == 1)
        B(:,o)=A(:,x);
        CB(o)=C(x);
        DB(o)=D(x);
        L(o)=x; %% vector of basic elements indecies
        o=o+1;
    else
        J(o2)=x; %% vector of nonbasic elements indecies
        o2=o2+1;
    end
end

while (key)
    switch key
        case 1
            B_inv=inv(B);
            XB=B_inv*b;
            z=-CB*XB;
            
            for j=1:n-m
                zj_cj(j)=CB*B_inv*A(:,J(j))-C(J(j));
            end
            if all(zj_cj<=0)
                disp("solution is reached");
                disp("indeces of basic elements are");
                disp(L);
                disp("their value is");
                disp(XB);
                disp("optimum value is:");
                disp(z);
                break;
            end
            [entering,I]= max(zj_cj);%% get index of maximum positive ratio
            B_invPj=B_inv*A(:,I);
            
            if any(B_invPj > 0)
                M=(B_inv*b)./B_invPj;
                [leaving, E]=min(M);
                temp=J(I);
                J(I)=L(E); %% switching indeces between leaving and entering
                L(E)=temp; %%vectors in basic and nonbasic indeces vectors
                for x=1:length(L)
                    B(:,x)=A(:,L(x));
                    CB(x)=C(L(x));
                end
            else
                disp(" solution is unbounded");
                break;
            end
            
        case 2
            B_inv=inv(B);
            XB=B_inv*b;
            w=-DB*XB;
            
            for j=1:n-m
                wj_dj(j)=DB*B_inv*A(:,J(j))-D(J(j));
            end
            if all(wj_dj<=0)
                disp("solution is reached");
                disp("indeces of basic elements are");
                disp(L);
                disp("their value is");
                disp(XB);
                disp("optimum value is:");
                disp(w);
            end
            [entering,I]= max(wj_dj);%% get index of maximum positive ratio
            B_invPj=B_inv*A(:,I);
            
            if any(B_invPj > 0)
                M=(B_inv*b)./B_invPj;
                [leaving, E]=min(M);
                temp=J(I);
                J(I)=L(E); %% switching indeces between leaving and entering
                L(E)=temp; %%vectors in basic and nonbasic indeces vectors
                for x=1:length(L)
                    B(:,x)=A(:,L(x));
                    DB(x)=D(L(x));
                    CB(x)=C(L(x));
                end
            else
                disp(" solution is unbounded");
                break;
            end
            if w==0
                key=1;
            end
    end
end
tic;toc
