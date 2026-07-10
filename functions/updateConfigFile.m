function updateConfigFile(configfile,out)

markertemplate="..\Models\Olympia\Olympia.markerlist";
jointtemplate="..\Models\Olympia\Olympia.kinematics";
segmenttemplate="..\Models\Olympia\Olympia.coord";

marker=readstruct(markertemplate,"FileType","xml");
joint=readstruct(jointtemplate,"FileType","xml");
segment=readstruct(segmenttemplate,"FileType","xml");
%configfile="G:\.shortcut-targets-by-id\1BBYI8sbh2BtgOertmVK5pnw2kMLSlbIq\OLIMPIA BIG\PROVE PALAINDOOR\ViconDataNew\TransTibial\Francesco Loragno\S04\TT05_S04_A.xml"

file=readstruct(configfile,"FileType","xml");

%markers
% mkr=[marker.Marker.labelAttribute];
% 
% for i=1:numel(file.MarkerSet.Marker)
%     lab=file.MarkerSet.Marker(i).labelAttribute;
%     idx=find(matches(mkr,lab));    
%     if isscalar(idx)
%       file.MarkerSet.Marker(i)=marker.Marker(idx);
%       fprintf("Marker: %s updated!\n",lab)
%     end
% end

% fix joint and segments showup
f=fieldnames(joint);
i=1;
for ii=1:numel(f)-1
    tmp=joint.(f{ii}).Joint;
    alljoints(i:i+numel(tmp)-1)=tmp;
    i=i+numel(tmp);
end
f=fieldnames(segment);
i=1;
for ii=1:numel(f)
    tmp=segment.(f{ii}).Segment;
    allsegments(i:i+numel(tmp)-1)=tmp;
    i=i+numel(tmp);
end
% joints
jts=[alljoints.labelAttribute];

for i=1:numel(file.KinematicModel.Joint)
    lab=file.KinematicModel.Joint(i).labelAttribute;
    idx=find(matches(jts,lab));    
    if isscalar(idx)
        for f=string(fieldnames(alljoints(idx)))'
            file.KinematicModel.Joint(i).(f)=alljoints(idx).(f);
        end
      fprintf("Joint: %s updated!\n",lab)
    end
end

%Segments
segs=[allsegments.labelAttribute];

for i=1:numel(file.KinematicModel.Segment)
    lab=file.KinematicModel.Segment(i).labelAttribute;
    idx=find(matches(segs,lab));    
    if isscalar(idx)
        for f=string(fieldnames(allsegments(idx)))'
            file.KinematicModel.Segment(i).(f)=allsegments(idx).(f);
        end
      fprintf("Segment: %s updated!\n",lab)
    end
end


% for i=1:numel(f.KinematicModel.Joint)
%     lab=f.KinematicModel.Joint(i).labelAttribute;
%     if ismember(lab,segs)
%       f.KinematicModel.Joint(i).signAttribute=alljoints(segs==lab).signAttribute;
%       f.KinematicModel.Joint(i).aliasAttribute=alljoints(segs==lab).aliasAttribute;
%       fprintf("Joint: %s sign and alias updated!\n",lab)
%     end
% end
% 
% segs=[allsegments.labelAttribute];
% for i=1:numel(f.KinematicModel.Segment)
%     lab=f.KinematicModel.Segment(i).labelAttribute;
%     if ismember(lab,segs)
%       f.KinematicModel.Segment(i).signAttribute=allsegments(segs==lab).signAttribute;
%       f.KinematicModel.Segment(i).aliasAttribute=allsegments(segs==lab).aliasAttribute;
%       fprintf("Segment: %s sign and alias updated!\n",lab)
%     end
% end

if nargin==1
    [p,ff,e]=fileparts(configfile);
else
    [p,ff,e]=fileparts(out);
end

    writestruct(file,fullfile(p,ff+"v2"+e),"FileType","xml")
    fprintf("%s saved in:\n %s\n",ff+"v2"+e,p);