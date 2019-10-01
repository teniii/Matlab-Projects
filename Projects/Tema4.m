%% problema 2

f=@(x,y)x^2+y^2-r^2;
g=@(x,y)x^2/8-y;

fimplicit(f,[-3,3,-3,3])
hold on
grid on
fimplicit(g,[-3,3,-3,3])
axis equal

syms x y
F1=x.^2+y.^2-4;
F2=x.^2/8 -y;
J=[diff(F1,x) diff(F1,y); diff(F2,x) diff(F2,y)];
subs(J,{x,y},{1,2});
J=matlabFunction(J,'vars',{x,y});
J(1,2)
F=@(x,y)[x.^2+y.^2-4;x.^2/8-y];
F(1,2)
    
x0=[-2;0];
eps=10^(-5);
tic;
[xaprox,N]=Newton(F,J,x0,0.01)
t=toc
plot(xaprox(1),xaprox(2),'o','MArkerFaceColor','g','MarkerSize',10)
x0=[2,0];
[xaprox,N]=Newton(F,J,x0,eps)
plot(xaprox(1),xaprox(2),'o','MArkerFaceColor','g','MarkerSize',10)

%% Interpolarea Lagrange

x=linspace(-pi/2,pi/2,5);
y=sin(x);
a=MetDirecta(x,y)
plot(x,y,'o','MArkerFaceColor','y','MarkerSize',10)
syms X;
Pn=0;
n=length(a)-1;
for i=1:n+1
    Pn=Pn+a(i)*X^(i-1);
end
hold on
Pn=matlabFunction(Pn,'vars',X);
fplot(Pn,[-pi/2,pi/2]);


%% comanda fplot
r=2;
y1=@(x)sqrt(r^2-x.^2);
y2=@(x)-sqrt(r^2-x^2);

%fplot(f,[xmin,xmax])
figure(1)
fplot(y1,[-2,2])
hold on
fplot(y2,[-r,r])
axis equal
grid on
x=@(t)r*cos(t);
y=@(t)r*sin(t);
%fplot(x,y,[tmin,tmax])
fplot(x,y,[0,2*pi])
grid on 
axis equal
F=@(x,y)x^2+y^2-r^2;
%fimplicit(F,[xmin,xmax,ymin,ymax])
fimplicit (F,[-r,r,-r,r])

f=@(x,y)x^2+y^2;
g=@(x,y)x^2/8-y;
fimplicit(f,[-3,3,-3,3])
fimplicit(g,[-3,3,-3,3])

F=@(x,y)[x.^2+y.^2-r^2;x.^2/8-y];
syms x y
f1=x.^2+y.^2-4;
f2=x.^2/8-y;
J=[diff(f1,x) diff(f1,y);diff(f2,x) diff(f2,y)];
subs(f1,{x,y},{1,2})
J=matlabFunction(J,'vars',{x,y});
J(1,2)
x0=[-3;0];
eps=10*(-5);
[xaprox,N]=Newton(F,J,x0,eps);
plot(xaprox(1),xaprox(2),'o','MarkerFaceColor','g','MarkerSize',10)

%% interp Lagrange
% 1.Metoda directa
% Date de intrare:f, x-vector
%Date de iesire: Pn
f=@(x)sin(x);
x=linspace(-pi/2,pi/2,5);
y=f(x);
plot(x,y,'o','MarkerFaceColor','g','MarkerSize',10)
hold on
tic;
a=MetDirecta(x,y);
syms X
Pn=0;
for i=1:length(a)
    Pn=Pn+a(i)*X^(i-1);
end
Pn=matlabFunction(Pn,'vars',{X});
fplot(Pn,[-pi/2,pi/2]);
t=toc;

function [a]=MetDirecta(x,y)
n=length(x)-1;
for i=1:n+1
    A(i,1)=1;
end
for i=1:n+1
    for j=2:n+1
        A(i,j)=x(i)^(j-1);
    end
end
a=A\transpose(A)
end


function [xaprox,N]=Newton(F,J,x0,eps)
   k=1;
    x(:,k)=x0;
    while true
        k=k+1;
        %z=J(x(1,k-1),x(2,k-1))\(-F(x(1,k-1),x(2,k-1)));
        z=GPT(J(x(1,k-1),x(2,k-1)),-F(x(1,k-1),x(2,k-1)));
        x(:,k)=x(:,k-1)+z;
        if norm(z,2)<eps
            break;
        end
    end
    xaprox=x(:,k);
    N=k;
end

function [x] = GPT(A,b)
  n = length(b);
  A = [A b];
  index = 1:n;
  for k=1:n-1
      apm = 0;
      for i=k:n
          for j=k:n
              if abs(A(i,j))>apm
                  apm = abs(A(i,j));
                  p = i;
                  m = j;
              end
          end
      end
      if apm==0
        disp 'Sistem incompatibil sau sistem compatibil nedeterminat';
        return;
      end
      if p~=k
          A([p k],:) = A([k p],:);
      end
      if m~=k
          A(:,[m k]) = A(:,[k m]);
          index([m k]) = index([k m]);
      end
      for l=k+1:n
          mlk = A(l,k)/A(k,k);
          A(l,:) = A(l,:) - mlk*A(k,:);
      end
  end
  if A(n,n)==0
      disp 'Sistem incompatibil sau sistem compatibil nedeterminat';
      return;
  end
      xtemp=SubsDesc(A(:,1:n),A(:,n+1));
  for i=1:n
      x(index(i)) = xtemp(i);
  end
  x = x';
end

function [x] = SubsDesc(A, b)
   n = length(b);
   x(n) = b(n)/A(n,n);
   for k=n-1:-1:1
      suma = 0;
      for j=k+1:n
          suma = suma + A(k,j)*x(j);
      end
      x(k) = (b(k) - suma)/A(k,k);
   end
   x = x';
end
