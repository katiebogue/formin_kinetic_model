function output = getstat(obj,NameValueArgs)
% GETSTAT retrieve polymer stat value at any level from a lookuptable
% 
% output = LOOKUPTABLE.GETSTAT(NameValueArgs)
%
%   Inputs (NameValueArgs): 
%       N       : (double) FH1 length
%       type    : (string) filament type ('single','double','dimer','ratio')
%       Stat    : (string) polymer stat
%       iSite   : (double) PRM location
%       Fil     : (string) filament, a or b
%       NName   : (string) N in string format, with "N" at front, reduces
%       runtime if added
% 
%   Outputs in order of precedence: 
%       If type is not specified, output is FilType 
%       If N is not specified, output is struct of N__
%       If Stat is not specified, output is struct of polymer stats
%       If Fil is not specified, output is Filament
%       If iSite is not specified, output is numerical array
% 
%   If N is beyond the maximum N values for the relevent filament types,
%   extrapolation is used.
% 
% See also LOOKUPTABLE, LOOKUPTABLE/EXTRAPOLATE.
    arguments
        obj Lookuptable
        NameValueArgs.N double =[]
        NameValueArgs.type string {mustBeMember(NameValueArgs.type,{'single','double','dimer','ratio'})}
        NameValueArgs.Stat string
        NameValueArgs.iSite double=[]
        NameValueArgs.Fil string {mustBeMember(NameValueArgs.Fil,{'a','b'})}
        NameValueArgs.NName string
    end
    N=NameValueArgs.N;
    if isfield(NameValueArgs,"type")
        type=NameValueArgs.type;
    else
        type=[];
    end
    if isfield(NameValueArgs,"Stat")
        Stat=NameValueArgs.Stat;
    else
        Stat=[];
    end
    iSite=NameValueArgs.iSite;
    if isfield(NameValueArgs,"Fil")
        Fil=NameValueArgs.Fil;
    else
        Fil=[];
    end

    obj.holdratio=true;
    
    tempNNames=obj.NNames;
    tempAddedNs=obj.AddedNs;
    tempStatNames=obj.StatNames;
    tempAddedStats=obj.AddedStats;

    if ~isempty(N)
        if ~isempty(iSite) && iSite>N
            error("iSite is larger than FH1 length")
        end
        if isfield(NameValueArgs,"NName")
            NName=NameValueArgs.NName;
        else
            NName=strcat("N",num2str(N));
        end
        if isempty(Stat)
            objcopy=obj.copytable;
            objcopy.holdratio=true;
            
            if ismember(NName,tempNNames)
                if ~isempty(type)
                    if isempty(objcopy.missingNs.(type))
                        tempout=tempAddedNs.(NName).(type);
                        if ~isempty(Fil) && type~="single"
                            for i=1:length(tempStatNames)
                                if class(tempout.(tempStatNames(i)))=="Filament"
                                    tempout.(tempStatNames(i))=tempout.(tempStatNames(i)).(Fil);
                                end
                            end
                        end
                    else
                        for i=1:length(tempStatNames)
                            extrap=obj.extrapolate(tempStatNames(i));
                            tempout.(tempStatNames(i))=generateextrapolation(type); %function already accounts for Fil
                        end
                    end
                else
                    tempout=tempAddedNs.(NName);
                    if ~isempty(objcopy.missingNs.single)
                        for i=1:length(tempStatNames)
                            extrap=obj.extrapolate(tempStatNames(i));
                            tempout.single.(tempStatNames(i))=generateextrapolation("single");
                        end
                    end
                    if ~isempty(objcopy.missingNs.double)
                        for i=1:length(tempStatNames)
                            extrap=obj.extrapolate(tempStatNames(i));
                            tempout.double.(tempStatNames(i))=generateextrapolation("double");
                        end
                    elseif ~isempty(Fil)
                        for i=1:length(tempStatNames)
                            if class(tempout.double.(tempStatNames(i)))=="Filament"
                                tempout.double.(tempStatNames(i))=tempout.double.(tempStatNames(i)).(Fil);
                            end
                        end
                    end
                    if ~isempty(objcopy.missingNs.dimer)
                        for i=1:length(tempStatNames)
                            extrap=obj.extrapolate(tempStatNames(i));
                            tempout.dimer.(tempStatNames(i))=generateextrapolation("dimer");
                        end
                    elseif ~isempty(Fil)
                        for i=1:length(tempStatNames)
                            if class(tempout.dimer.(tempStatNames(i)))=="Filament"
                                tempout.dimer.(tempStatNames(i))=tempout.dimer.(tempStatNames(i)).(Fil);
                            end
                        end
                    end
                end
            else
                for i=1:length(tempStatNames)
                    extrap=obj.extrapolate(tempStatNames(i));
                    if ~isempty(type)
                        tempout.(tempStatNames(i))=generateextrapolation(type);
                    else
                        if i==1
                            tempout=FilType;
                        end
                        tempout.single.(tempStatNames(i))=generateextrapolation("single");
                        tempout.double.(tempStatNames(i))=generateextrapolation("double");
                        tempout.dimer.(tempStatNames(i))=generateextrapolation("dimer");
                    end
                end
            end
        else
            if ismember(NName,tempNNames)
                if isempty(type)
                    tempout=FilType;
                    if ~isempty(obj.missingNs.single) && ismember(NName,obj.missingNs.single)
                        extrap=obj.extrapolate(Stat);
                        tempout.single=generateextrapolation("single");
                    else
                        tempout.single=tempAddedStats.(Stat).single.(NName);
                    end

                    if ~isempty(obj.missingNs.double) && ismember(NName,obj.missingNs.double)
                        extrap=obj.extrapolate(Stat);
                        tempout.double=generateextrapolation("double");
                    else
                        tempout.double=tempAddedStats.(Stat).double.(NName);
                        if ~isempty(Fil)
                            if class(tempout.double)=="Filament" 
                                tempout.double=tempout.double.(Fil);
                            end
                        end
                    end

                    if ~isempty(obj.missingNs.dimer) && ismember(NName,obj.missingNs.dimer)
                        extrap=obj.extrapolate(Stat);
                        tempout.dimer=generateextrapolation("dimer");
                    else
                        tempout.dimer=tempAddedStats.(Stat).dimer.(NName);
                        if ~isempty(Fil)
                            if class(tempout.dimer)=="Filament"
                                tempout.dimer=tempout.dimer.(Fil);
                            end
                        end
                    end
                else
                    if ~isempty(obj.missingNs.(type))
                        tempout=tempAddedStats.(Stat).(type).(NName);
                        if ~isempty(Fil)
                            if class(tempout)=="Filament"
                                tempout=tempout.(Fil);
                            end
                        end
                    else
                        extrap=obj.extrapolate(Stat);
                        tempout=generateextrapolation(type);
                    end
                end
            else
                extrap=obj.extrapolate(Stat);
                if ~isempty(type)
                    tempout=generateextrapolation(type);
                else
                    tempout=FilType;
                    tempout.single=generateextrapolation("single");
                    tempout.double=generateextrapolation("double");
                    tempout.dimer=generateextrapolation("dimer");
                end
            end
        end
    else
        objcopy=obj.copytable;
        objcopy.holdratio=true;

        if isempty(Stat)
            if isempty(type)
                tempout=objcopy;
                if ~isempty(Fil)
                    tempout.double=getonefil("double");
                    tempout.dimer=getonefil("dimer");
                    tempout.updateAddedNs();
                    tempout.updateAddedStats();
                end
            else
                if isempty(Fil) || type=="single"
                    tempout=objcopy.(type);
                else
                    tempout=getonefil(type);
                end
            end
        else
            tempout=tempAddedStats.(Stat);
            if ~isempty(type)
                if isempty(Fil) || type=="single"
                    tempout=tempout.(type);
                else
                    tempout=getonefil(type,Stat);
                end
            elseif ~isempty(Fil)
                 tempout.double=getonefil("double",Stat);
                 tempout.dimer=getonefil("dimer",Stat);
            end
        end
    end

    if isempty(iSite)
        output=tempout;
    else
        output=assigniSite(tempout);
    end

    obj.holdratio=false;

    function out=assigniSite(input)
        inclass=class(input);
        if inclass=="struct" 
            out=input;
            fields=fieldnames(out);
            for ii=1:length(fields)
                out.(fields{ii})=assigniSite(input.(fields{ii}));
                if fields{ii}(1)=="N"
                    [fieldnum,tf] = str2num(fields{ii}(2:end));
                    if tf
                        if fieldnum<iSite
                            out=rmfield(out,fields{ii});
                        end
                    end
                end
            end

        elseif inclass=="FilType" 
            out=input;
            out.single=assigniSite(input.single);
            out.double=assigniSite(input.double);
            out.dimer=assigniSite(input.dimer);
        elseif inclass=="Filament"
            out=input;
            out.a=assigniSite(input.a);
            out.b=assigniSite(input.b);
        elseif inclass=="Lookuptable"
            out=input;
            fields=fieldnames(out);
            props=metaclass(input).PropertyList;
            for ii=1:length(fields)
                if ~props(ii).Dependent
                    out.(fields{ii})=assigniSite(input.(fields{ii}));
                end
                if fields{ii}(1)=="N"
                    [fieldnum,tf] = str2num(fields{ii}(2:end));
                    if tf
                        if fieldnum<iSite
                            out=rmfield(out,fields{ii});
                        end
                    end
                end
            end
        else
            size1=size(input,1);
            size2=size(input,2);
            if size1==1 && size2==1
                out=input;
            elseif size1==1 && size2>1
                if size2<iSite
                    out=input;
                else
                    out=input(:,iSite);
                end
            elseif class(input)=="string"
                out=input;
            else
                error("Entries do not match iSite format")
            end
        end
    end
    

    function out= getonefil(typ,instat)
        out=struct;
        Names=objcopy.NList.(typ);
        for ii=1:length(Names)
            Nnum=Names(ii);
            Nname=strcat("N",num2str(Nnum));
            if nargin>1
                value=objcopy.(typ).(Nname).(instat);
                if class(value)=="Filament"
                        out.(Nname)=value.(Fil);
                else
                    out.(Nname)=value;
                end
            else
                stats=fieldnames(objcopy.(typ).(Nname));
                for jj=1:length(stats)
                    value=objcopy.(typ).(Nname).(stats{jj});
                    if class(value)=="Filament"
                        out.(Nname).(stats{jj})=value.(Fil);
                    else
                        out.(Nname).(stats{jj})=value;
                    end
                end
            end
        end
    end

    function out= generateextrapolation(typ)
        exttype=extrap.(typ);
        if class(exttype)=="scatteredInterpolant"
            out=zeros(1,N);
            for jj=1:N
                out(1,jj)=exttype(N,jj);
            end
        elseif class(exttype)=="Filament"
            if isempty(Fil)
                if class(exttype.a)=="scatteredInterpolant"
                    out=Filament(zeros(1,N),zeros(1,N));
                    for jj=1:N
                        out.a(1,jj)=exttype.a(N,jj);
                        out.b(1,jj)=exttype.b(N,jj);
                    end
                else
                    out=Filament(exttype.a(N),exttype.b(N));
                end
            else
                if class(exttype.a)=="scatteredInterpolant"
                    out=zeros(1,N);
                    for jj=1:N
                        out(1,jj)=exttype.(Fil)(N,jj);
                    end
                else
                    out=exttype.(Fil)(N);
                end
            end
        else
             out=exttype(N);
        end
    end
end