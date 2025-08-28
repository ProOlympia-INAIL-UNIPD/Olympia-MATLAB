function addRecord(rootdir,trial,static,configuration)
file=fullfile(rootdir,'ProcessedTrialList.xlsx');
trial=string(trial);
exp="(\w+)(\d\d)_S(\d\d)_T(\d\d)_";

code=regexp(trial,exp,'tokens');
code=string(code{:});

code=[code(1),code(1)+code(2)];
varname=["Key","Classification","Subject","FootModel","FootCat","SocketName","SocketAlignment","Trial","Static","Configuration"]
ID=height(file)+1;
prompt=["Classification","Subject","FootModel","FootCat","SocketName","SocketAlignment"];
def=[code(1),code(2),"",nan,"",nan];
answer=inputdlg(prompt,"Enter Trial details",[1 30],def);
answer={ID,answer{:},trial,static,configuration};
T=cell2table(answer);
T.Properties.VariableNames=varname;
writetable(T,file,"WriteMode","append");
