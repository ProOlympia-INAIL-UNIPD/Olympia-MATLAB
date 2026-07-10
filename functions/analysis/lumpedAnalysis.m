function trial=lumpedAnalysis(trial,point)
arguments
    trial Trial
end
arguments (Repeating)
    point string
end
point=string(point);

fp=trial.ForcePlatform;%.resample(trial.Metadata.POINT.RATE);
fp=fp.combineFP;
sr=trial.NSamples/trial.NFrames;
R=align2ISB(fp); %this align running direction to ISB
AP=1;
V=2;
ML=3;
GRF=fp.GRF(1:sr:end,:)*R';
mask=GRF~=0;
COP=fp.COP(1:sr:end,:)*R';

ev=trial.Events.exportEvents('point',false);

for label=point
    try
        P=trial.Points(matches([trial.Points.Label],label));
        PP=P.Coordinates;
        % x=1:length(PP(:,AP));
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
        theta=atan2d(L(:,AP),L(:,V));
        
        maskLR=mask&false;
        %% Parameters Calculation
        for i=1:length(FO) %repeat for each step
            try
                %To avoid NaN values in theta and L: find the first/last value
                %where the GRF exists
                idxFS = find(mask(FS(i):FO(i)), 1, 'first');
                FS(i) = FS(i) + idxFS - 1;
                idxFO = find(mask(FS(i):FO(i)), 1, 'last');
                FO(i) = FS(i) + idxFO - 1;

                maskLR(FS(i):FO(i))=true;

                s.L0(i)=Lnorm(FS(i));                       % L0 at FS
                s.Lfin(i)=Lnorm(FO(i));                     % L at FO
                s.ThetaFS(i)=theta(FS(i));
                s.ThetaFO(i)=theta(FO(i));

                deltaL=s.L0(i) - Lnorm(FS(i):FO(i));
                [s.GRFp_max(i), Imax]=max(GRFmod(FS(i):FO(i)));
                s.K(i)=s.GRFp_max(i)/deltaL(Imax);             %K=GRFp_max/(L0-L(Imax))
            catch
            end

        end
        maskLR=maskLR&mask;
        trial.Scalars(end+1)=Scalar(label+"-COP_Length",sprintf("Leg length (%s)",P.Units),maskLR.*Lnorm.*[1 0 0]);
        trial.Scalars(end+1)=Scalar(label+"-COP_GRFmod","Force projection (N)",maskLR.*GRFmod.*[1 0 0]);
        %trial.Scalars(end+1)=Scalar("K_"+label+"-COP","N"+P.Units,maskLR.*K.*[1 0 0]);
        trial.Scalars(end+1)=Scalar(label+"-COP_Theta",'degrees',maskLR.*theta.*[1 0 0]);
        trial.Metadata.("LUMPED1D_"+label+"_COP")=s;
        clear s
    catch ME
        warning("Lumped Analysis (%s): %s",label,ME.message);
    end
end