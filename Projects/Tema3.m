 %% Exercitiul 1

A=[0 1 1;2 1 5;4 2 1];
b=[3;5;1];
GaussFaraPiv(A,b)
GaussPivPart(A,b)
GaussPivTot(A,b)

A=[0 1 -2;1 -1 1;1 0 -1];
b=[4;6;2];
GaussFaraPiv(A,b)
GaussPivPart(A,b)
GaussPivTot(A,b)

%% Exercitiul 3

%c)
eps=10^(-20);
A=[eps 1;1 1];
b=[1;2];
GaussFaraPiv(A,b)
GaussPivPart(A,b)

C=10^(20);
A=[1 C;1 1];
b=[C;2];
GaussPivPart(A,b)
GaussPivTot(A,b)

%Se poate observa ca Metoda GFP nu functioneaza pentru valori foarte mici
%ale coeficientilor, iar GPP nu functioneaza pentru valori foarte mari ale
%coeficientilor

%% Exercitiul 5

A=[1 2 -1;2 4 7;-1 2 5];
b=[2;13;10];
[L,U,w]=LU_GPP(A);
n=length(w);
for i=1:n-1
    b([i,w(i)])=b([w(i),i]);
end
y=SubsAsc(L,b)
x=SubsDesc(U,b)

%% Exercitiul 8

A=[1 2 3; 2 5 8;3 8 14];
b=[-5;-14;-25];
L=FactCholesky(A)
y=SubsAsc(L,b)
x=SubsDesc(L',y)

%% Metoda Gauss fara pivotare
function [x] = GaussFaraPiv(A,b)

p=-1;
n=size(A,1);
A=[A,b];
for k=1:n-1
    for j=k:n
        if A(j,k)~=0
            p=j;
            break;
        end
    end

    if p==-1
        fprintf('Sist. incomp. sau comp. nedet');
        x='error';
        return;
    end
    if p~=k
        A([p,k],:)=A([k,p],:);
    end
    for l=k+1:n
        mlk=A(l,k)/A(k,k);
        A(l,:)=A(l,:)-mlk*A(k,:);
    end
end

if A(n,n) == 0
    fprintf('Sist. incomp. sau comp. nedet.');
    x='error';
    return;
end

x=SubsDesc(A(1:n,1:n),A(:,n+1));
end

%% Metoda Gauss cu pivotare partiala

function [x] = GaussPivPart(A,b)

n=size(A,1);
A=[A,b];

for k=1:n-1 
    max = abs(A(k,k));
    p=k;
   for i=k:n
       
           if abs(A(i,k)) > max
               max = abs(A(i,k));
               p=i;
           
       end
   end
if(max==0)
    fprintf('Sistem incompatibil sau nedeterminat');
    x='error';
    return;
end
if p~=k
    A([p,k], :) = A([k,p], :); 
end

for l=k+1:n
    M(l,k) = A(l,k)/A(k,k);
    A(l,:) = A(l,:) - M(l,k)*A(k,:);
end
end

if A(n,n) == 0
    fprintf('Sistem incompatibil sau nedeterminat');
    x='error';
    return;
end

x = SubsDesc(A(1:n, 1:n), A(:, n+1));

end


%% Metoda Gauss cu pivotare totala

function [x] = GaussPivTot(A,b)
n=length(b); %n=size(A,1)
index=1:n;
A=[A,b];
for k=1:n-1
   max = abs(A(k,k));
   for i=k:n
       for j=k:n
           if abs(A(i,j)) > max
               max = abs(A(i,j));
               p=i;
               m=j;
           end
       end
   end
if(max==0)
    fprintf('Sistem incompatibil sau nedeterminat');
    x='error';
    return;
end
if p~=k
    A([p,k], :) = A([k,p], :); 
end

if m~=k
    A(:, [m,k]) = A (:, [k,m]);
    index([m,k])=index([k,m]);
end

for l=k+1:n
    mlk = A(l,k)/A(k,k);
    A(l,:) = A(l,:) - mlk*A(k,:);
end
end

if A(n,n) == 0
    fprintf('Sistem incompatibil sau nedeterminat');
    x='error';
    return;
end

y = SubsDesc(A(1:n, 1:n), A(:, n+1));

for i=1:n
    x(index(i)) = y(i);
end
end

%% Metoda Substitutiei Descendente

function [x] = SubsDesc(A,b)
n = length(b);
x(n) = 1/A(n,n) * b(n);
k = n - 1;

while k>0
    sum=0;
    for j=k+1:n
        sum = sum + A(k,j)*x(j);
    end
    x(k) = 1/A(k,k) * (b(k) - sum);
    k=k-1;
end
end

%% Metoda Substitutiei Ascendente

function [x] = SubsAsc(A,b)
n = length(b);
x(1) = 1/A(1,1) * b(1);
k=1;
for k=2:n-1
    
    sum = 0;
    for j=1:k-1
        sum=sum + A(k,j)*x(j);
    end
    x(k) = 1/A(k,k)*(b(k) - sum);
end
end

%% Factorizarea LU cu Gauss pivotare partiala

function [L,U,w] = LU_GPP(A)
n = size(A,1);
L=eye(n);
for k=1:n-1
    max = abs(A(k,k));
       for j=k:n
           if abs(A(j,k)) > max
               max = abs(A(j,k));
               p=j;
           end
       end
           w(k)=p;
           if p~=k
                A([p,k], :) = A([k,p], :);
           end
           for l=k+1:n 
               llk=A(l,k)/A(k,k);
               A(l,:)=A(l,:)-llk*A(k,:);
           end
           if k>1
               L([p,k],1:k-1)=L([k,p],1:k-1);
           end
       
end
       if A(n,n)==0
           fprintf('Sist. incomp. sau comp. nedet.');
           x='error';
           return;
       end
       U=A;
end

%% Fatorizarea Cholesky
function [L]=FactCholesky (A)

n=size(A,1);
alfa=A(1,1);
if alfa<= 0
    fprintf('A nu admite factorizare Cholesky');
    L='error';
    return;
end
L(1,1)=sqrt(alfa);
for i=2:n
    L(i,1)=A(i,1)/L(1,1);
end
for k=2:n
    sum=0;
    for s=1:k-1
        sum=sum+L(k,s)*L(k,s);
    end
    alfa=A(k,k)-sum;
    if alfa <=0
        fprintf('A nu admite factorizare Cholesky');
        L='error';
        return;
    end
    L(k,k)=sqrt(alfa);
    for i=k+1:n
        sum=0;
        for s=1:k-1
            sum=sum+L(i,s)*L(k,s);
        end
        L(i,k)=(A(i,k)-sum)/L(k,k);
    end
end

end

