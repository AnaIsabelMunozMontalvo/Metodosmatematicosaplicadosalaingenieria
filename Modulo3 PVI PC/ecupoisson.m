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



function [u,x,y,error]=ecupoisson(a,b,c,d,dx,dy,f,...
                       cc,solexacta)
%este código resuelve la ecuación de Poisson bidimensional 
%-Delta(u)=f (x,y)en (a,b)x(c,d)
%con condiciones de frontera Dirichlet, dadas por la función cc(x,y)
%utilizando el esquema de cinco punto en cruz para aproximar el operador
%laplaciano Delta.
%también calcula el error cometido cuando se conoce la solución exacta,
%definida en solexac(x,y)
%en caso de no conocer la solucion exacta, no se mete como argumento
%los tamaños de discretización vienen dados por dx y dy

if nargin == 8
     solexacta=@(x,y) 0.*x;
end

%número de intervalos en x, nx, en y, ny 
nx=(b-a)/dx; ny=(d-c)/dy;
%número de nodos
nx1=nx+1;
dx2=dx^2; dy2=dy^2;
kii=2/dx2+2/dy2;     kix=-1/dx2;  kiy=-1/dy2;
dim=(nx+1)*(ny+1);   K=speye(dim,dim);
rhs=zeros(dim,1);
rhs1=zeros(dim,1);
y = c;

%calculamos la matriz de coeficientes
for m = 2:ny
    x = a; y = y + dy;
    for n = 2:nx
        i = n+(m-1)*(nx+1);
        x = x + dx;
        rhs(i) = feval(f,x,y);
        K(i,i) = kii;         K(i,i-1) = kix;
        K(i,i+1) =kix; K(i,i+nx1)=kiy;
        K(i,i-nx1)=kiy;
      end
    end
    
 %pasamos los valores de la condición de contorno al lado derecho
        
rhs1(1:nx1) = feval(cc,x,c);
rhs1(dim-nx:dim) = feval(cc,x,d);
y = [c:dy:d];
rhs1(1:nx1:dim-nx) = feval(cc,a,y);
rhs1(nx1:nx1:dim) = feval(cc,b,y);
rhs = rhs - K*rhs1;
nbound = [[1:nx1],[dim-nx:dim],...
          [1:nx1:dim-nx],[nx1:nx1:dim]];
ninternal = setdiff([1:dim],nbound);
%resolvemos para los nodos internos
K = K(ninternal,ninternal);
rhs = rhs(ninternal);
utemp = K\rhs;
uh = rhs1;
%imponemos la condición de  contorno en los bordes
uh (ninternal) = utemp;
k = 1; y = c;
for j = 1:ny+1
    x = a;
    for i = 1:nx1
        u(i,j) = uh(k);
        k = k + 1;
        ue(i,j) = feval(solexacta,x,y);
        x = x + dx;
    end
    y = y + dy;
end
x = [a:dx:b];
y = [c:dy:d];
if nargout == 4
  if nargin == 8
     warning('la solución exacta no se conoce');
     error = [ ];
  else
     error = max(max(abs(u-ue)))/max(max(abs(ue)));
  end
end
return



