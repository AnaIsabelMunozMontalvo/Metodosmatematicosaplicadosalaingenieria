%____copyright___="Copyright (C) 2022 A. Nolla, A.I. Mu�oz, E. Schiavi."
%____license____="GPL-3.0-only"

%Detalles sobre el empleo de los c�digos y ejercicios de aplicaci�n pueden
%encontrarse en las direcciones de los siguientes documentos elaborados por
%A. Nolla, A.I. Mu�oz, E. Schiavi:
%https://burjcdigital.urjc.es/handle/10115/20132
%https://burjcdigital.urjc.es/handle/10115/20134
%As� como en el fichero Readme

% La mayor parte de los c�digos de la colecci�n presentada en esta librer�a
%son adaptaciones de los publicados en el libro "C�lculo cient�fico con 
% Matlab y Octave" de A. Quarteroni y F. Saliery, 
% que se pueden obtener en https://mox.polimi.it/qs/.




function [vectorsol,itera] = metnewtonsistema(fecusistema,jacobiana,...
                                vectorx0,errorper,maxitera)
% este c�digo obtiene una soluci�n aproximada del sistema no lineal
% dado por el vector de funciones fecusistema, F, igualado a cero, F=0
% partiendo de la semilla vectorx0
% jacobiana es la funci�n matriz jacobiana J
% errorper es el m�ximo error permitido, que se mide tomando el m�dulo
%del vector -[J^-1(vectorx)]F(vectorx)
% vectorsol es la aproximaci�n del vector ra�z que hemos obtenido
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
    fprintf(['El metodo no converge en el n�mero m�ximo ',...
       'de iteraciones',...
       'por no alcanzarse un error inferior al permitido'],F);

end
return
