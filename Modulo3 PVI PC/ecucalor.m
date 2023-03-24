%____copyright___="Copyright (C) 2022 A. Nolla, A.I. Muñoz, E. Schiavi."
%____license____="GPL-3.0-only"

%Detalles sobre el empleo de los códigos y ejercicios de aplicación pueden
%encontrarse en las direcciones de los siguientes documentos elaborados por
%A. Nolla, A.I. Muñoz, E. Schiavi:
%https://burjcdigital.urjc.es/handle/10115/20132
%https://burjcdigital.urjc.es/handle/10115/20134
%Así como en el fichero Readme

% La mayor parte de los códigos de la colección presentada en esta librería
%son adaptaciones de los publicados en el libro "Cálculo científico con 
% Matlab y Octave" de A. Quarteroni y F. Saliery, 
% que se pueden obtener en https://mox.polimi.it/qs/.



function [x,uf]=ecucalor(C,intespacio,intiempo,pasosespacio, ...
               pasostiempo,theta,u0,cc,f)
%este código resuelve la ecuación del calor en una dimesión espacial
%utilizando un theta-método para la discretización temporal
%du/dt=C (d^2 u/dx^2) +f  (t0,T)x(a,b)
%intespacio=(a,b) intiempo=(t0,T)
%la condición inicial se deota por u0
%la condición de contorno tipo Dirichlet, la denotamos por cc
%theta es el valor de theta empleado en el theta método
%theta=0 es para Euler explícito
%theta=1 es para Euler implícito
%theta=0.5 para Crank-Nicolson
%pasosespacio es el número de pasos en espacio
%pasostiempo es el número de pasos en tiempo
%x es el vector que contiene los nodos en espacio
%uf contiene los valores numéricos obtenidos en t=T

h  = (intespacio(2)-intespacio(1))/pasosespacio;
dt = (intiempo(2)-intiempo(1))/pasostiempo;
N = pasosespacio+1;%número de nodos en espacio
%perparamos la matriz de coeficientes, que variará dependiendo del
%theta-método empleado
e =ones(N,1); D = spdiags([-e 2*e -e],[-1,0,1],N,N);
I = speye(N); M = I+C*dt*theta*D/h^2;
Mn = I-C*dt*(1-theta)*D/h^2;
 M(1,:) = 0;
 M(1,1) = 1; M(N,:) = 0; M(N,N) = 1;
 x =linspace(intespacio(1),intespacio(2),N);
 x = x'; fn = feval(f,intiempo(1),x); 
 un =feval(u0,x);
 [L,U]=lu(M); 
 for t = intiempo(1)+dt:dt:intiempo(2)
    fn1 = feval(f,t,x);
    %lado derecho
    rhs = Mn*un+dt*(theta*fn1+(1-theta)*fn);
    temp = feval(cc,t,[intespacio(1),intespacio(2)]);
    rhs([1,N]) = temp;
    %el sistema se resuleve utilizando una factorización LU
    u = L\rhs;
    u = U\u;
    fn = fn1;
    un = u;
 
 end
uf=u;
return
