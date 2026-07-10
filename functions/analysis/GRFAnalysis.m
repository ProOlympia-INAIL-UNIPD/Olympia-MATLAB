        function [trial, grf] = GRFAnalysis(trial,normalize)
         %GRFANALYSIS Computes impulses and GRF parameters
         %  [obj, grf] = GRFAnalysis(obj,normalize) computes impulses and GRF
         %  parameters from the forcePlatform data stored in the acquisition,
         %  assigning it into the Metadata.
         %  Inputs:
         %  trial         - Trial
         %  normalize   - Select what to use for mass normalization
         %  Outputs:
         %  obj - Updated Trial
         %  grf - structure containing the calculated GRF data

            arguments
            trial (1,1) Trial
            normalize string {mustBeMember(normalize,["no","bodymass","bodyweight"])}="bodyweight";
            end
        
        forceplat=trial.ForcePlatform;
        if numel(forceplat)==0
           warning("Trial Contains no Force Platforms!, GRF analysis can't be performed!");
           return
        end   

        
        if btkGetEventNumber(trial.C3DHandle)==0 %needed because metadata isn't updated (btk bug?)
           warning("Trial Contains no events!, GRF analysis can't be performed!");
           return
        else
           events=trial.Events;
        end

 

        try %import mass from trial and adjust units
        mass=trial.Metadata.ANTROPOMETRY.mass;
        if trial.Metadata.ANTROPOMETRY.units=="t"
           mass=mass*1000;
        end
        mustBePositive(mass);
        switch normalize
            case "no"
                mass=1;
                unit="N%s";
            case "bodymass"
                %mass=mass;
                unit="N%s/kg";
            case "bodyweight"
                mass=mass*9.81;
                unit="N%s/N";
        end
        catch
            mass=1;
            unit="N%s";
        end

   
        GRF=sum(cat(3,forceplat.GRF),3)/mass;
        R=forceplat.align2ISB();

        GRF=GRF*R'; %rotate the GRF to the new directions
        AP=1; %in ISB is X, first column
        V=2;  % in ISB is Y, second column
        
        
        f_force=forceplat.SampleRate;
        events=events.selectFootContacts;
        events=exportEvents(events,'analog',false);
        events=rmfield(events,'units');
        
        for ctx=string(fieldnames(events))' %run the analysis for each context of events
            if isempty(fieldnames(events.(ctx)))
                continue
            end
            try
            FS=events.(ctx).Foot_Strike;
            FO=events.(ctx).Foot_Off;
            for j=length(FS):-1:1
            fc=FS(j);
            fo=FO(j);
            
            %%-MAX, MIN, MEAN
            % anteroposterior
            grf.GRF_PARAM.Units=sprintf(unit,"");
            grf.GRF_PARAM.(ctx+"AnteroPosteriorMean")(j)=mean(GRF(fc:fo,AP));
            grf.GRF_PARAM.(ctx+"AnteroPosteriorMin")(j)=min(GRF(fc:fo,AP));
            grf.GRF_PARAM.(ctx+"AnteroPosteriorMax")(j)=max(GRF(fc:fo,AP));
            % vertical
            grf.GRF_PARAM.(ctx+"VerticalMax")(j)=max(GRF(fc:fo,V));
            grf.GRF_PARAM.(ctx+"VerticalMean")(j)=mean(GRF(fc:fo,V));
            %%-IMPULSES
            grf.GRF_IMPULSES.Units=sprintf(unit,"s");
            grf.GRF_IMPULSES.(ctx+"Vertical")(j)=trapz(GRF(fc:fo,V)-1)/f_force; %(Fy-BW)/BW= Fy/BW-1
            %%-Braking/Propulsive impulses X for each step of the same limb
            GRF_AP_e=GRF(fc:fo,AP);
            zerocross=[ne(diff(sign(GRF_AP_e)),0);false]; %find zero crossing
            Ic=cumtrapz(GRF_AP_e)/f_force;
            Ic=[Ic(zerocross)', Ic(end)];
            In=Ic-[0,Ic(1:end-1)];
            
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorStartPropulsion")(j)=sum(In)-(sum(In(In<0))+max(In));
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorBraking")(j)=sum(In(In<0));
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorPropulsive")(j)=max(In);
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorNet")(j)=sum(In);    

            end
            catch ME
                warning(ME.message)
            end
        end
            if exist('grf','var')
                if isfield(grf,"GRF_PARAM")
            trial.Metadata.GRF_PARAM=grf.GRF_PARAM;
                end
                if isfield(grf,"GRF_IMPULSES")
            trial.Metadata.GRF_IMPULSES=grf.GRF_IMPULSES;
                end
            trial.setC3DMetaData;
            else
               warning("GRF analysis not successfull, check events and forceplate data!")
               grf=struct();
            end
        end