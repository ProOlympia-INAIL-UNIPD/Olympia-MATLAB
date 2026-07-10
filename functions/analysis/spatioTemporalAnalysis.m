function SPATIOTEMP=spatioTemporalAnalysis(trial)
         %SPATIOTEMPORALANALYSIS Computes spatial-temporal parameters
            arguments
            trial (1,1) Trial
            end

        forceplat=trial.ForcePlatform;
        if numel(forceplat)==0
           warning("Trial Contains no Force Platforms!, GRF analysis can't be performed!");
           return
        end   

        
        if getEventCount(trial.Events)==0 %needed because metadata isn't updated (btk bug?)
           warning("Trial Contains no events!, GRF analysis can't be performed!");
           return
        else
           events=trial.Events;
        end

        fp=forceplat.combineFP;
   
        %GRF=sum(cat(3,forceplat.GRF),3)/mass;
        COP=fp.COP;
        R=forceplat.align2ISB();

        %GRF=GRF*R'; %rotate the GRF to the new directions
        COP=COP*R';
        AP=1; % now the first direction has the most variability (i.e., is the running direction)
        %ML=2;
        V=2;  % in a typical flat laboratory, the vertical should change least (i.e., is vertical)
        events=events.selectFootContacts;
    ev_force=events.exportEvents("analog");
    ev_point=events.exportEvents("point");
    ev_time=events.exportEvents("seconds");
side=["Left", "Right"];
ss=["L","R"];
for s=1:2
    lr=side(s);                  %considered side
    cl=side(not(side==side(s))); %controlateral
    %stance time
    try
        SPATIOTEMP.(lr+"ContactTime")=ev_time.(lr).Foot_Off-ev_time.(lr).Foot_Strike;
    catch
    end

    %flight time
    try
        FO_lr=ev_time.(lr).Foot_Off;
        FS_cl=ev_time.(cl).Foot_Strike;
        FS_cl(FS_cl<FO_lr(1))=[];
        ns=min(numel(FO_lr),numel(FS_cl));
        SPATIOTEMP.(lr+"FlightTime")=FS_cl(1:ns)-FO_lr(1:ns);
    catch
    end
    % stance length 
    try
        P=trial.Points.PointStruct.(ss(s)+"GT")*R';
    catch
        P=trial.Points.PointStruct.(ss(s)+"HJC")*R';
    end

    P=P(:,AP);
    try
        SPATIOTEMP.(lr+"ContactLengthGT")=P(ev_point.(lr).Foot_Off)-P(ev_point.(lr).Foot_Strike);
    catch
    end
    

    % stride length
    try
        FS_lr=ev_force.(lr).Foot_Strike;
        FS_cl=ev_force.(cl).Foot_Strike;
        FS_cl(FS_cl<FS_lr(1))=[];
        ns=min(numel(FS_lr),numel(FS_cl));
        COP_lr=COP(FS_lr(1:ns),AP);
        COP_cl=COP(FS_cl(1:ns),AP);

        SPATIOTEMP.(lr+"StrideLength")=COP_cl-COP_lr;
    catch
    end

    try %avgspeed
        v=(trial.Points([trial.Points.Label]=="LGT").Velocity+trial.Points([trial.Points.Label]=="RGT").Velocity)/2;
        v=mean(v,'omitnan');
        SPATIOTEMP.("AvgSpeed")=v;
    catch
    end
    
    end
trial.Metadata.SPATIOTEMP=SPATIOTEMP;
trial.setC3DMetaData;

end
