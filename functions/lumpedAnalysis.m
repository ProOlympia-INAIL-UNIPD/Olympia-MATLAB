function trial=lumpedAnalysis(trial,point)
arguments
    trial Trial
end
arguments (Repeating)
    point string
end
point=string(point);
fp=trial.ForcePlatform.resample(trial.Metadata.POINT.RATE);
fp=fp.combineFP;
R=align2ISB(fp);
GRF=fp.GRF*R';
mask=GRF~=0;
COP=fp.COP*R';

ev=trial.Events.exportEvents('point',false);

for label=point
try
P=trial.Points(matches([trial.Points.Label],label));
PP=P.Coordinates;
% x=1:length(PP(:,1));
% xq=1:fp.NSamples;
% PP=spline(x,PP',xq)';
PP=PP*R';
gr=P.Group;
if gr=="General"
   if extract(P.Label,1)=="R"
      gr="Right";
   elseif extract(P.Label,1)=="L"
      gr="Left";
   end
end
FS=ev.(gr).Foot_Strike;
FO=ev.(gr).Foot_Off;
% GRF(mask)=nan;
% PP(mask)=nan;
% COP(mask)=nan;
%tracks
L=PP-COP;
Lnorm=vecnorm(L,2,2);
Lvers=L./Lnorm;
GRFmod=dot(Lvers,GRF,2);
K=GRFmod./Lnorm;
theta=atan2d(L(:,1),L(:,3));

% for i=length(FS):-1:1
% 
% 
% % L{i}=PP(FS(i):FO(i),:)- COP(FS(i):FO(i),:); %vector PP-COP
% % Lnorm{i}=vecnorm(L{i},2,2);                 % ||PP-COP||
% % Lvers{i}=L{i}./Lnorm{i};                    %versor PP-COP
% % GRF_stance{i}=GRF(FS(i):FO(i),:);           %GRF track during stance
% end
maskLR=mask&false;
%% Parameters Calculation
for i=1:length(FS) %repeat for each step
    if mask(FS(i))==0
       FS(i)=FS(i)+1;
    end
    if mask(FO(i))==0
       FO(i)=FO(i)-1;
    end
    maskLR(FS(i):FO(i))=true;
    %s.L{i}=Lnorm{i};                           % L during stance
    s.L0(i)=Lnorm(FS(i));                       % L0 at FS
    s.Lfin(i)=Lnorm(FO(i));                     % L at FO
    s.ThetaFS(i)=theta(FS(i));
    s.ThetaFO(i)=theta(FO(i));
    s.Lmin(i)=min(Lnorm(FS(i):FO(i)));          % minimum L during step
    s.deltaAPprox(i)=PP(FO(i),1)-PP(FS(i),1);   % forward displacement of proximal point relative to initial position
    % s.deltaL{i}=Lnorm{i}(1)-Lnorm{i};           % dL=L0-L
    s.V0(i)=PP(FS(i),2);                        % vertical position of PP at FS
    s.Vfin(i)=PP(FO(i),2);                      % vertical position of PP at FS
    s.Vmin(i)=min(PP(FS(i):FO(i),2));           % minimum vertical position of PP during stance 
    s.deltaVfin(i)=PP(FS(i),2)-PP(FO(i),2);     % vertical trajectory of PP during stance
    % s.deltaV{i}=PP(FS(i):FO(i),2)-PP(FS(i),2);  % vertical trajectory of PP relative to initial position
    %s.modGRF{i}=dot(Lvers{i},GRF_stance{i},2);  % magnitude of GRF projection on leg versor
    %s.GRFcomp{i}=s.modGRF{i}.*Lvers{i};         % components GRF projection on leg versor
    [GRFmax, Imax]=max(GRFmod(FS(i):FO(i)));
    %s.GRFmax(i)=GRFmax;                         % maximum GRF during stance
    % s.K(i)=GRFmax/s.deltaL{i}(Imax);            % leg stiffness at GRFmax [N/mm]
    % s.Knorm(i)=s.K(i)/(9.81*mass);              % leg stiffness at GRFmax [(N/BW)/mm]
    % s.Kist{i}=(s.modGRF{i})./s.deltaL{i};       % istantaneous leg stiffness [N/mm]
    % s.Kist{1,i}(s.Kist{1,i}==Inf)=nan;

end
maskLR=maskLR&mask;
trial.Scalars(end+1)=Scalar("L_"+label+"-COP",sprintf("Leg length (%s)",P.Units),maskLR.*Lnorm.*[1 0 0]);
trial.Scalars(end+1)=Scalar("GRFmod_"+label+"-COP","Force projection (N)",maskLR.*GRFmod.*[1 0 0]);
trial.Scalars(end+1)=Scalar("K_"+label+"-COP","N"+P.Units,maskLR.*K.*[1 0 0]);
trial.Scalars(end+1)=Scalar("Theta_"+label+"-COP",'degrees',maskLR.*theta.*[1 0 0]);
trial.Metadata.("LUMPED1D_"+label+"_COP")=s;
catch ME
    warning("Lumped Analysis (%s): %s",label,ME.message);
end
end