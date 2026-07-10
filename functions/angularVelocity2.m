function    [omega, omegaMagnitude] = angularVelocity2(R,freq)


%----------------------INTRODUZIONE------------------------ 
% Questa funzione permette di calcolare il vettore velocità angolare dalla matrice di posa di un segmento
%----------------------INPUT-------------------------------
% R [3 x 3 x nFrame]
%----------------------OUTPUT------------------------------
% omega [nFrame x 3]
% omegaMagnitude [nFrame x 1]
%----------------AUTHOR------------
% Andrea G. Cutti, DEIS - University of Bologna, INAIL - Prosthesis Centre (c).
% Software provided under Open Source Licence.
%----------------------------------


nFrame=size(R,3);

%Calcolo dei termini di derivata

a11(1,:) = R(1,1,:);
a11 = a11';

a12(1,:) = R(1,2,:);
a12 = a12';

a13(1,:) = R(1,3,:);
a13 = a13';

a21(1,:) = R(2,1,:);
a21 = a21';

a22(1,:) = R(2,2,:);
a22 = a22';

a23(1,:) = R(2,3,:);
a23 = a23';

a31(1,:) = R(3,1,:);
a31 = a31';

a32(1,:) = R(3,2,:);
a32 = a32';

a33(1,:) = R(3,3,:);
a33 = a33';


Rdot(1,1,:) = signalDerivative2( a11, 1, 'cspline',freq);
Rdot(1,2,:) = signalDerivative2( a12, 1, 'cspline',freq);
Rdot(1,3,:) = signalDerivative2( a13, 1, 'cspline',freq);
Rdot(2,1,:) = signalDerivative2( a21, 1, 'cspline',freq);
Rdot(2,2,:) = signalDerivative2( a22, 1, 'cspline',freq);
Rdot(2,3,:) = signalDerivative2( a23, 1, 'cspline',freq);
Rdot(3,1,:) = signalDerivative2( a31, 1, 'cspline',freq);
Rdot(3,2,:) = signalDerivative2( a32, 1, 'cspline',freq);
Rdot(3,3,:) = signalDerivative2( a33, 1, 'cspline',freq);


for i=1:nFrame
    %Modificato 23/11
    %Rdot : glo_loc, R = glo_loc. 
    %Opz 1
    %Rdot*R'= glo_loc*loc_glo = fisso
    S=Rdot(:,:,i)*R(:,:,i)';

    %Opz 2
    %R'*Rdot = loc_glo*glo_loc = mobile
    %Oppure dovrei lasciare opz 1 e poi moltipicare ancora per R(glo_loc)
    %S=R(:,:,i)'*Rdot(:,:,i);
    omega(i,:)=[ S(3,2) S(1,3) S(2,1)]; % espressa in radianti
    omegaMagnitude(i)=norm(omega(i,:))*180/pi; % espressa in gradi
end


