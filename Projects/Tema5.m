%% Tema lab05
%% Exercitiul 5

X = [-pi/2,-pi/6,pi/6,pi/2];
f = @(x)sin(x);
fp = @(x)cos(x);
Y = f(X);
Z = fp(X);
x = linspace(-pi/2,pi/2,100);

% 1)
%a)
y = MetNeville(X,Y,x)
%b)
y = MetNDD(X,Y,x)
%c)
[y,z] = MetHermite(X,Y,Z,x)

% 2)
figure(1);
grid on;

subplot(2,2,1);
plot(x, f(x),'y','LineWidth',3);
title('sin(x)');

subplot(2,2,2);
plot(x, MetNeville(X,Y,x),'m','LineWidth',3);
title('Metoda Neville');

subplot(2,2,3);
plot(x, MetNDD(X,Y,x),'c','LineWidth',3);
title('Metoda Newton cu Diferente Divizate');

subplot(2,2,4);
plot(x, y,'r','LineWidth',3);
title('Metoda Hermite');


figure(2);
grid on;

subplot(2,2,1);
plot(x, fp(x),'y','LineWidth',3);
title('cos(x)');

subplot(2,2,2);
plot(x, z,'m','LineWidth',3);
title('Metoda Hermite z');

% 3)
figure(3);
grid on;
plot(x, abs(f(x)-y),'k','LineWidth',3);
title('Eroarea f(x)- Metoda Hermite');


%% Exercitiul 8


f = @(x)sin(x);
fp = @(x)cos(x);
X = [-pi/2,-pi/6,pi/6,pi/2];
Y = f(X);
Z = fp(X);
x = linspace(-pi/2,pi/2,100);

% 1)

[y,z] = MetHermiteDD(X,Y,Z,x);

% 2)
figure(4);
grid on;

subplot(2,2,1);
plot(x, f(x),'y','LineWidth',3);
title('sin(x)');

subplot(2,2,2);
plot(x, y,'m','LineWidth',3);
title('Metoda Hermite cu Diferente Divizate');


figure(5);
grid on;

subplot(2,2,1);
plot(x, fp(x),'y','LineWidth',3);
title('cos(x)');

subplot(2,2,2);
plot(x, z,'m','LineWidth',3);
title('Metoda Hermite cu Diferente Divizate z');

% 3)
figure(6);
grid on;
plot(x, abs(f(x)-y),'k','LineWidth',3);
title('Eroarea f(x)- Metoda Hermite cu Diferente Divizate');
[y,z] = MetHermiteDD(X,Y,Z,pi/2)

%% Algoritmi functii folosite

function [y, z] = MetHermite(X, Y, Z, x)

  n = length(X)-1;
  Her = 0;
  HerD = 0;
  
  for k=1:n+1
      Lpk = zeros(size(x));
      Ld = zeros(size(x));
      L = ones(size(x));
      produs = ones(size(x));
      numitor = ones(size(x));
      asemenea = zeros(size(x));
      
      for i=1:n+1
        if i~=k
            for m=1:length(x)
                if (x(m)-X(i))~=0
                    produs(m) = produs(m) * (x(m)-X(i));
                else
                    asemenea(m) = 1;
                end 
            end
            numitor = numitor .* (X(k)-X(i));
        end
      end
      
      for i=1:n+1
        if i~=k
          L = L .* (x-X(i))./(X(k)-X(i));
          Lpk = Lpk + 1./(X(k)-X(i));
          for m=1:length(x)
              if (x(m)-X(i))~=0
                  if asemenea(m)==0
                      Ld(m) = Ld(m) + produs(m)/((x(m)-X(i))*(numitor(m)));
                  end
              else
                  Ld(m) = Ld(m) + produs(m)/numitor(m);
              end
          end
        end
      end
      
      H = L.*L.*(1-2.*Lpk.*(x-X(k)));
      K = L.*L.*(x-X(k));
      Her = Her + H.*Y(k) + K.*Z(k);
      Hd = 2.*L.*Ld.*(1-2.*Lpk.*(x-X(k))) - L.*L.*2.*Lpk;
      Kd = 2.*L.*Ld.*(x-X(k)) + L.*L;
      HerD = HerD + Hd.*Y(k) + Kd.*Z(k);
      
  end
  
  y = Her;
  z = HerD;
end


function [y,z] = MetHermiteDD(X,Y,Z,x)

    n = length(X)-1;
    
    for i=1:n+1
        XB(2*i-1) = X(i);
        XB(2*i) = X(i);
    end
    
    Q = zeros(2*n+2);
    
    for i=1:n+1
        Q(2*i-1,1) = Y(i);
        Q(2*i,1) = Y(i);
        Q(2*i,2) = Z(i);
        
        if i>=2
            Q(2*i-1,2) = (Q(2*i-1, 1) - Q(2*i-2,1)) / (XB(2*i-1) - XB(2*i-2));
        end
    end
    
    for i=3:2*n+2
        for j=3:i
            Q(i,j) = (Q(i,j-1) - Q(i-1,j-1)) / (XB(i) - XB(i-j+1));
        end
    end
    
    for index=1:length(x)
      y(index) = Q(1,1);
      z(index) = 0;
      
      for k=2:2*n+2
          sumaprod = 0;
          asemenea=0;
          produs = 1;
          produsDiv = 1;
          
          for m=1:k-1
              produs = produs*(x(index)-XB(m));
              if (x(index)-XB(m))~=0
                  produsDiv = produsDiv*(x(index)-XB(m));
              else
                 asemenea= asemenea + 1; 
              end
          end
          
          for m=1:k-1
            if (x(index)-XB(m)) ~= 0
                if asemenea==0
                  sumaprod = sumaprod + produsDiv/(x(index)-XB(m));
                end
            else
                if asemenea==1
                  sumaprod = sumaprod + produsDiv;
                end
            end  
          end
          
          z(index) = z(index) + sumaprod*Q(k,k);
          y(index) = y(index) + Q(k,k)*produs;
      end
    end
end

function [y] = MetNDD(X,Y,x)

  n = length(X)-1;
  Q = zeros(n+1);
  
  for i=1:n+1
    Q(i,1) = Y(i);
  end
  
  for i=2:n+1
      for j=2:i
          Q(i,j) = (Q(i,j-1) - Q(i-1,j-1)) / (X(i)-X(i-j+1));
      end
  end
  
  for index=1:length(x)
      y(index) = Q(1,1);
      for k=2:n+1
          produs = 1;
          for z=1:k-1
              produs = produs*(x(index)-X(z));
          end
          y(index) = y(index) + Q(k,k)*produs;
      end
  end
end

function [y] = MetNeville(X,Y,x)

  n = length(X)-1;
  Q = zeros(n+1);
  
  for index=1:length(x)
      for i=1:n+1
          Q(i,1) = Y(i);
      end
      
      for i=2:n+1
          for j=2:i
              Q(i,j) = ((x(index)-X(i-j+1))*Q(i,j-1)-(x(index)-X(i))*Q(i-1,j-1)) / (X(i) - X(i-j+1));
          end
      end
      
      y(index) = Q(n+1,n+1);
  end
end