%% Exercitiul 1

syms x;
f=8*x^3+4*x-1;
f=matlabFunction(f,'vars',x);
a(1)=0;
b(1)=1;
x(1)=(a(1)+b(1))/2;

for k=2:3
    if f(x(k-1))==0 
        break;
    elseif f(a(k-1))*f(x(k-1))<0
        a(k)=a(k-1);
        b(k)=x(k-1);
        x(k)=(a(k)+b(k))/2;
    elseif f(a(k-1))f(x(k-1))>0 
        a(k)=x(k-1);
        b(k)=b(k-1);
        x(k)=(a(k)+b(k))/2;
    end
end
X_aprox=x(k)
eps=(b(1)-a(1))/(2^2)

%% Exercitiul 2

syms x;
f=x.^3-7*x.^2+14*x-6;
f=matlabFunction(f,'vars',x);
y=linspace(0,4,75);
figure(2)
plot(y,f(y),'--r','Linewidth',1.5);
grid on
hold on
xa=MetBisectiei(f,0,1,10^(-5))
plot(xa,f(xa),'o');
xa=MetBisectiei(f,1,3.2,10^(-5))
plot(xa,f(xa),'o');
xa=MetBisectiei(f,3.2,4,0.00001)
plot(xa,f(xa),'o');


%% Exercitiul 3
syms x;
y=exp(x)-2;
y=matlabFunction(y,'vars',x);
z=cos(exp(x)-2);
z=matlabFunction(y-z,'vars',x);
xa=MetBisectiei(z,0.5,1.5,10^(-5))
a=linspace(0.5,1.5,100);
figure(3)
plot(a,z(a));
grid on
hold on
plot(xa,z(xa),"o")


%% Exercitiul 4
syms x;
y=x^2-3;
y=matlabFunction(y,'vars',x);
xa=MetBisectiei(y,0,10,10^(-5))

%% Exercitiul 5
syms x;
f=x^3-7*x^2+14*x-6;
g=diff(f);
f=matlabFunction(f,'vars',x);
g=matlabFunction(g,'vars',x);
xa=MetNR(f,g,2,0.00001)
y=linspace(0,4,100);
figure(5)
plot(y,f(y))
grid on
hold on
plot(xa,f(xa),"o")
%sirul generat converge spre radacina xa=3 deoarece functia este
%descrescatoare in 2
xa=MetNR(f,g,1,0.00001)
plot(xa,f(xa),"o")
%ca sirul generat sa convearga la radacina din intervalul [0;2.5] luam
%punctul 1, functia fiind crescatoare in 1


%% Exercitiul 7
syms x;
f=8*x^3+4*x-1;
g=diff(f)
% derivata functiei este pozitiva pentru orice x pe R => functia este
% strict crescatoare => solutie unica
f=matlabFunction(f,'vars',x);
g=matlabFunction(g,'vars',x);
y=linspace(-1,1,1000);
figure(7)
plot(y,f(y))
grid on
hold on
xa=MetNR(f,g,1,0.00001)
plot(xa,f(xa),"o")
MetSecantei(f,0,1,0,1,0.00001)
MetPozFalse(f,0,1,0.00001)

%%  Exercitiul 8
syms x;
f=x^3-18*x-10;
f=matlabFunction(f,'vars',x);
x=linspace(-5,5,300);
figure(8)
plot(x,f(x));
grid on
hold on
xs=MetSecantei(f,-5,-3,-4.5,-4,0.00001)
plot(xs,f(xs),"o")
xs=MetSecantei(f,-3,0,-3,0,0.00001)
plot(xs,f(xs),"o");
xs=MetSecantei(f,0,5,1,5,0.00001)
plot(xs,f(xs),"o");

xf=MetPozFalse(f,-5,-3,0.00001)
plot(xf,f(xf),"o")
xf=MetPozFalse(f,-3,0,0.00001)
plot(xf,f(xf),"o")
xf=MetPozFalse(f,0,5,0.00001)
plot(xf,f(xf),"o")

%In metoda Secantei numarul de iteratii pentru fiecare radacina sunt: 6, 8,
%9, iar in metoda Pozitiei false numarul de iteratii este de 10, 6,
%respectiv 9. Se poate observa ca in unele situatii o metoda este mai
%eficienta ca cealalta, in alte situatii cealalta este mai eficienta sau
%ambele sunt la fel de eficiente.



%% Algoritmi metode folosite

function [xaprox] = MetPozFalse(f,a,b,eps)

k=1;
A(1)=a;
B(1)=b;
x(1)=(a*f(b)-b*f(a))/(f(b)-f(a));
k=k+1;
    if f(A(k-1))*f(x(k-1))<0
        A(k)=A(k-1);
        B(k)=x(k-1);
        x(k)=(A(k)*f(B(k))-B(k)*f(A(k)))/(f(B(k))-f(A(k)));
    elseif f(A(k-1))*f(x(k-1))>0
        A(k)=x(k-1);
        B(k)=B(k-1);
        x(k)=(A(k)*f(B(k))-B(k)*f(A(k)))/(f(B(k))-f(A(k)));
    end
while abs(x(k)-x(k-1))/abs(x(k-1)) > eps
   
    k=k+1;
    if f(x(k-1))==0
        x(k)=x(k-1);
        break;
    elseif f(A(k-1))*f(x(k-1))<0
        A(k)=A(k-1);
        B(k)=x(k-1);
        x(k)=(A(k)*f(B(k))-B(k)*f(A(k)))/(f(B(k))-f(A(k)));
    elseif f(A(k-1))*f(x(k-1))>0
        A(k)=x(k-1);
        B(k)=B(k-1);
        x(k)=(A(k)*f(B(k))-B(k)*f(A(k)))/(f(B(k))-f(A(k)));
    end
end
xaprox=x(k);
end

function [xaprox] = MetNR(f,g,x0,eps)

k=1;
x(1)=x0;
k=k+1;
x(2)=x0 - f(x0)/g(x0);

while abs(x(k)-x(k-1))/abs(x(k-1)) > eps
    k=k+1;
    x(k)=x(k-1) - f(x(k-1))/g(x(k-1));
end
xaprox=x(k);
end

function [xaprox] = MetSecantei(f,a,b,x0,x1,eps)
k=2;
x(1)=x0;
x(2)=x1;
while abs(x(k)-x(k-1))/abs(x(k-1)) >= eps

    k=k+1;
    x(k)=(x(k-2)*f(x(k-1))-x(k-1)*f(x(k-2)))/(f(x(k-1)) - f(x(k-2)));
    if x(k)<a || x(k)>b
        msgbox('Introduceti alte valori pentru x0,x1');
        xaprox=x(k);
        return;
    end
end
xaprox=x(k);
end


function [xaprox] = MetBisectiei(f,a,b,eps)

A(1)=a;
B(1)=b;
x(1)=(A(1)+B(1))/2;
n=log2((b-a)/eps);

for k=2:n+1
    if f(x(k-1))==0 
        break;
    elseif f(A(k-1))*f(x(k-1))<0
        A(k)=A(k-1);
        B(k)=x(k-1);
        x(k)=(A(k)+B(k))/2;
    elseif f(A(k-1))f(x(k-1))>0 
        A(k)=x(k-1);
        B(k)=B(k-1);
        x(k)=(A(k)+B(k))/2;
    end
end
xaprox=x(k);
end

        