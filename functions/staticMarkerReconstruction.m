function [stat]=staticMarkerReconstruction(stat,configfile)
%[stat]=staticMarkerReconstruction(stat,configfile) uses the information
%stored in the configuration file to add the virtual and calibrated markers
%to an acquisition, returning the modified trial.
p=fileparts(stat.c3dfile);
if isempty(stat.ConfigFile)&&not(exist("configfile","var"))
    error('Unable to configure a trial without an internal or external configfile');
end

if isfield(stat.ConfigFile.MarkerSet,"CalibratedMarkerDef")
mkr=reshape(stat.ConfigFile.MarkerSet.CalibratedMarkerDef,1,[]);
else
    mkr=[];
end
mainpath=mfilename('fullpath');
mainpath=fileparts(mainpath);
wands=readstruct(fullfile(mainpath,'..', 'CalibrationObjects','XML','Wands.calobj'),"FileType","xml");
for cm=mkr
    m=stat.ConfigFile.MarkerSet.Marker(matches([stat.ConfigFile.MarkerSet.Marker.("label"+stat.XMLatt)],cm.("label"+stat.XMLatt)));    
    if length(m)~=1
       warning("Calibrated marker with label: %s is not part of the defined markerset and won't be created!",cm.("label"+stat.XMLatt));
    else
        if isfile(char(fullfile(p,cm.("source"+stat.XMLatt))))
        else
        [cm.("source"+stat.XMLatt),p]=uigetfile(fullfile(p,'*.c3d'),"Select Wand Trial for "+cm.("label"+stat.XMLatt));
        if cm.("source"+stat.XMLatt)==0
            p=fileparts(stat.c3dfile);
            continue
           
        end
        stat.ConfigFile.MarkerSet.CalibratedMarkerDef(matches([stat.ConfigFile.MarkerSet.CalibratedMarkerDef.("label"+stat.XMLatt)],cm.("label"+stat.XMLatt))).("source"+stat.XMLatt)=cm.("source"+stat.XMLatt);
        end
        w_tr=Trial(char(fullfile(p,cm.("source"+stat.XMLatt))));
        try
        wand=wands.(cm.("wand"+stat.XMLatt));
        catch
        error('Wand object not existing!');
        end
        try
        
        stat=stat.wandReconstruct(w_tr,wand,cm.("label"+stat.XMLatt));
        catch ME
            warning("Error occurred during reconstruction of %s",cm.("label"+stat.XMLatt)+ME.message);            
        end
           
    end
end

stat=stat.reconstructVirtualMarkers;

expectedpoints=[stat.ConfigFile.MarkerSet.Marker.("label"+stat.XMLatt)];
existingpoints=[stat.Points.Label];
missingpoints=expectedpoints(~matches(expectedpoints,existingpoints));
while ~isempty(missingpoints)
    answer=questdlg(["The following points defined in the MarkerSet are missing from the current trial, what do you want to?"; missingpoints'],"Points Missing","Solidify from Trial","Remove","Ignore","Ignore");
    if isempty(answer)
       break
    end
    switch answer
        case "Solidify from Trial"
            [p,~,e]=fileparts(stat.c3dfile);
            [f,p]=uigetfile(fullfile(p,e));
            if f==0
               break
            end
            tr=Trial(fullfile(p,f),stat.ConfigFile);
            
            tr=tr.mean;
            pp=listdlg("ListString",missingpoints,"PromptString","Select point to retrieve from current trial","SelectionMode","single");
            pp=missingpoints(pp);
            try
            cs=tr.Points(matches([tr.Points.Label],pp)).Segment;
            catch
               warning(sprintf("%s unavailable in %s",pp,f))
               break
            end
            tr.Points(not(matches([tr.Points.Segment],cs)))=[];
            stat=stat.segmentSolidification(tr,'ignore');
        case "Remove"
            rmv=listdlg("PromptString","Select points to remove from configuration:","SelectionMode","multiple","ListString",missingpoints);
            stat.ConfigFile.MarkerSet.Marker(matches(expectedpoints,missingpoints(rmv)))=[];
        case "Ignore"
            break
    end
expectedpoints=[stat.ConfigFile.MarkerSet.Marker.("label"+stat.XMLatt)];
existingpoints=[stat.Points.Label];
missingpoints=expectedpoints(~matches(expectedpoints,existingpoints));
end
end