function [Options,OnsetTimes,Labels]=LoadSymptomsTable(SxOrAsx,Virus,OldOrNew,Symptoms,Options)
N=length(Symptoms.Virus);
Options.SxSubjectIDs=[];
Found=0;
for i=1:N
    if strcmp(Symptoms.Virus{i},Virus)
        if strcmp(Symptoms.OldOrNew{i},OldOrNew)
            if ~isempty(Symptoms.SxSubjects{i})
                if Symptoms.SxOrAsx(i)==SxOrAsx
                    Found=1;
                    Options.SxSubjectIDs=Symptoms.SxSubjects{i}(:,1);
                    OnsetTimes=Symptoms.SxSubjects{i}(:,2);
                    Labels=ones(Symptoms.TotalSubjects{i},1);
                    Labels(Options.SxSubjectIDs)=0;
                end
            end
        end
    end
end

Options.S=length(Options.SxSubjectIDs);
if ~Found
   disp('ERROR: Couldn t find Symptoms for the choice of (SxOrASx,Virus,OldOrNew,)') 
end

