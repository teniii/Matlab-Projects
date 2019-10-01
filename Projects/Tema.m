%Vectori si matrice
%%Operatii cu vectori si matrice
A=[1,2,3;
   4 5 6
   7 8 9]
B=[9,8,7;6 6 4; 3 2 1]
A/B             %am comparat si am observat ca A/B<=>A*inv(A)
A*inv(B)        
inv(A)*B        %la fel si A\B<=>inv(A)*B
A\B
A.\B
x=[7 4 8 10 12]
length(x)       %am calculat: lungimea vectorului
sum(x)                       %suma elementelor vectorului
prod(x)                      %produsul elementelor
det(A)                       %determinantul matricei A
size(A)                      %dimensiunile matricei A

~A              %operatii logice cu matrice
A&B
A|B
xor(A,B)

sin(A)          %functii matematice cu matrice
cos(A)
log(A)

%%Generarea vectorilor
C=[1:3;2:4;3:5]
D=[1:2:5; 4:2:8]

%%Lucrul cu matrice impartite in blocuri
E=[1:4;5:8;9:-1:6;5:-1:2]
E(1:2,2:3)          %afisam intersectia liniilor 1 2 cu col 2 3
E(2,:)
E(1:3,2:4)

%%Concatenarea matricelor
y=[[1 2 3] [4 5 6] 7]
B=[E [-1 -2 -3 -4]']    %concatenam matricea E cu vect coloana
B=[E ;[-1 -2 -3 -4]]    %concatenam matricea E cu vect linie
v=(1:3)'
m=[v v.^3 3.^v]

%%Functii ce creaza matrice speciale
zeros(2)        %afisam o matrice patratica de 2*2 cu 0
ones(3)         %afisam o matrice patratica de 3*3 cu 1
rand(2,3)       %afisam o matrice de 2*3 cu nr aleatoare
randn(3,2)    %afisam o matrice de 3*2 cu nr aleatoare cu distrib unif
eye(1,2)               %o matrice unitate de 1*2

%%Modificarea dimensiunilor matricelor
a=[1 2 3; 4 5 6]
reshape(a,3,2)
reshape(a,6,1)


%Instructiuni Matlab
%%Instructiunea for
for i=1:5
    a(i)=1/i;   %crearea unui vector cu valorile 1,1/2,...,1/5
end

j=[1:1:5]
for i=j         %acelasi lucru ca in primul for
    a(i)=1/i;
end

a=[[1 2 3]' [-2 -4 -6]']
for i=a
    i*i'
end

D={[1 2] [-1+2j,2+3j]}

for i=1:2
    D{i}        %afisarea componentelor celulelor
end

%%Instructiunea while
i=1;
while i<=5
    a(i)=1/i
    i=i+1;
end

b=[4 5 6; 4+i 5+2i 0]
i=1;
while i<=3 && all(b(:, i))
    b(:, i)     %afisarea coloanelor matricei cu elemente nenule
    i=i+1;
end

 %%Instructiunea if
test=-5;
if test>=0
    'variabila e pozitiva'
else
    'variabila este negativa'
end

%Reprezentari grafice
%%Reprezentarea curbelor plane

x=linspace(0,2*pi,100)  %reprezentarea grafica a functiilor sin si cos
y=cos(x)                %pe intervalul [0,2pi]
plot(x,sin(x),'--r','Linewidth',3,'MarkerFaceColor', 'y','MarkerSize', 10)
hold on
plot(x,y,'--b','Linewidth',3,'MarkerFaceColor', 'b','MarkerSize', 10)
xlabel('time')
grid
title('sin(time), cos(time)')

figure(2)
t = 0:0.01:4*pi;    %functia sin(t)/t pe [0,4pi]
x = sin(t) ./ t;
plot(t, x)
grid
xlabel('\bft[s]')
texstr = '$\frac{sin(t)}{t}$';
text('string',texstr,  'interpreter', 'latex', 'fontsize',40,...
             		'units','norm', 'pos',[.5 .5]);
                
        %functiile sint(t), cos(t), sin(t)^2 si cos(t) in functie de
        %sin(t)
t = 0: pi/20: 2*pi;
x = sin(t);
y = cos(t);
z = x .* x;

subplot(2, 2, 1)
plot(t, x)
grid
xlabel('t[s]')
ylabel('sin(t)')

subplot(2, 2, 2)
plot(t, y)
grid
xlabel('t[s]')
ylabel('cos(t)')

subplot(2, 2, 3)
plot(t, z)
grid
xlabel('t[s]')
ylabel('sin^2(t)')

subplot(2, 2, 4)
plot(x, y)
grid
xlabel('sin(t)')
ylabel('cos(t)')

