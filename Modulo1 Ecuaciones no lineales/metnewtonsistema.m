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




function [vectorsol,itera] = metnewtonsistema(fecusistema,jacobiana,...
                                vectorx0,errorper,maxitera)
% este código obtiene una solución aproximada del sistema no lineal
% dado por el vector de funciones fecusistema, F, igualado a cero, F=0
% partiendo de la semilla vectorx0
% jacobiana es la función matriz jacobiana J
% errorper es el máximo error permitido, que se mide tomando el módulo
%del vector -[J^-1(vectorx)]F(vectorx)
% vectorsol es la aproximación del vector raíz que hemos obtenido
% fecusistema y jacobiana deben ser definidas como funciones en scripts .m

itera = 0;
medidaerror = errorper + 1;
 vectorx= vectorx0;
while medidaerror >= errorper && itera < maxitera
    F = fecusistema(vectorx(1),vectorx(2));
    J = jacobiana(vectorx(1),vectorx(2));
    incremento = -inv(J)*F;
    vectorx = vectorx + incremento;
    medidaerror = norm(incremento);
    itera = itera + 1;
end
vectorsol=vectorx;

if (itera==maxitera && medidaerror> errorper)
    fprintf(['El metodo no converge en el número máximo ',...
       'de iteraciones',...
       'por no alcanzarse un error inferior al permitido'],F);

end
return
