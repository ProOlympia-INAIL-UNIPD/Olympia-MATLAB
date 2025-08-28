function [Fp, Mp] = DistDynCalc(Joint,Fd, Md, g, Fd_appPoint,inertiamask)
% Joint è il giunto sul quale sto calcolando la dinamica inversa
% Fd è forza distale in globale (N) [nF x 3]
% Md è momento distale in globale(Nm) [nF x 3]
% g è il vettore gravità
% Fd_appPoint punto applicazione forza distale [nF x 3]
% inertiamask è un selettore (1/0) per evitare l'uso delle masse e delle inerzie
debugmode=true;
if vecnorm(g)>10
    warning('Gravity vector magnitude (%0.2f) appears to be higher than expected, please check that g is provided in m/s2!',vecnorm(g));
end
JC=Joint.JointCenter;
dist=Joint.Child;
CM=dist.COM;
if isempty(CM)
    CM=JC;
    inertiamask=0;
end

b = Fd_appPoint-CM.Coordinates;
m=dist.Mass*inertiamask;
I=dist.Icm*inertiamask;
T=dist.TransformMat;
I_glo=pagemtimes(T(1:3,1:3,:),pagemtimes(I,'none',T(1:3,1:3,:),'transpose'));
au=dist.AngleUnits;
dist.AngleUnits="rad";
omega=dist.AngularVelocity;
alpha=dist.AngularAcceleration;
dist.AngleUnits=au;
a_com=CM.Acceleration;

b_Fp = JC.Coordinates - CM.Coordinates;
ma=m*a_com;
W=repmat(m*g, [size(T,3) 1]);
Fp = sum(cat(3,-Fd, ma, -W), 3, 'omitnan');
if debugmode==true
    figure()
    for i=1:3
        subplot(2,3,i)
        plot(-Fd(:,i));
        hold on
        plot(ma(:,i));
        plot(-W(:,i));
        plot(Fp(:,i),'k')
        title(sprintf('Force %c',char(87+i)))
        legend('F_D','ma','mg','F_P')
        hold off
    end
end
% Mp_loc+Md_loc + bp_loc X Fp_loc + bd_loc X Fd_loc = Icom * alpha_loc +omega_loc X I*omega_loc
MFd=cross(b, Fd);
MFp=cross(b_Fp, Fp);
wxIw=cross(omega, permute(pagemtimes(I_glo, permute(omega,[2 3 1])),[3 1 2]));
Ia=permute(pagemtimes(I_glo, permute(alpha,[2 3 1])),[3 1 2]);
Mp = sum(cat(3,...
    -MFd, ... %componente forza distale
    -MFp,... %componente forza prossimale (nulla se applicata a centro articolare, con polo in centro articolare)
    -Md, ... %momento distale (esterno se piede)
    wxIw,... %componente legata a quantità di moto
    Ia),...%componente inerziale rot
    3,'omitnan');

if debugmode==true

    for i=1:3
        subplot(2,3,i+3)
        plot(-Md(:,i));
        hold on
        plot(-MFd(:,i));
        plot(-MFp(:,i));
        plot(wxIw(:,i));
        plot(Ia(:,i));
        plot(Mp(:,i),'k')
        title(sprintf('Moment %c',char(87+i)))
        legend('M_D','b_D\timesF_D','b_P\timesF_P','\omega\timesI\omega','I\alpha','M_P')
        hold off
    end
    sgtitle("Inverse Dynamics - "+Joint.Label)
end

end