function [F,COP,Fc,COPc]=getROS(trial,segmentGroup,segmentLabel)

T=trial.Segments.(segmentGroup).(segmentLabel).TransformMat;

fp=trial.ForcePlatform.combineFP;

fp=fp.resample(trial.Metadata.POINT.RATE);


R=T(1:3,1:3,:);

COP=points2local(fp.COP,T);
F=fp.GRF;
for i=1:size(F,1)
    F(i,:)=R(:,:,i)*F(i,:)';
end

ev=trial.Events.exportEvents('point',false);

FS=ev.(segmentGroup).Foot_Strike;
FO=ev.(segmentGroup).Foot_Off;
ctcmask=false(size(F));
for f=length(FS):-1:1
    Fc(:,:,f)=time2cycle([],F(FS(f):FO(f),:),101);
    COPc(:,:,f)=time2cycle([],COP(FS(f):FO(f),:),101);
    ctcmask(FS(f):FO(f),:)=true;
end
F(ctcmask)=0;
COP(ctcmask)=0;