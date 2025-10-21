classdef Lookuptable <handle 
    %LOOKUPTABLE Holds formatted polymer statistics from polymer-c   
    %   
    %   Construction:
    %       obj = LOOKUPTABLE(ltvar) 
    %           Inputs:
    %               ltvar : (struct) contains single, double, and/or dimer
    %               fields full of N struct objects to be added to the
    %               lookuptable; struct should be created from
    %               makeLookupMat and reading in polymer-c outputs
    % 
    % See also MAKELOOKUPMAT. 

    properties (GetAccess=public, SetAccess=protected)
        single % struct with structs for single filament simulations
        dimer % struct with structs for dimer filament simulations
        double % struct with structs for double filament simulations
    end

    properties (SetAccess=private, SetObservable, AbortSet)
        StatNames % names of polymer stats/ outputs from simulation
        NNames % list of Ns that have been added
        AddedStats struct % each stat/ output as a FilType
        AddedNs struct % each N as a FilType
        NList FilType % array of N values as FilType
        max_length FilType % Max filament length before extrapolation
        missingNs FilType % List of missing Ns below the max_length 
        ratio struct % ratio of dimer/ double values
    end

    properties (SetAccess=protected)
        interpolant struct % contains extrapolation functions for stats (indexed by stat)
    end

    properties (SetAccess=protected,GetAccess=private)
        holdratio logical=false % if true, do not update ratio when called
        ratiostore % holds the ratio value
    end

    methods 
        function obj = Lookuptable(ltvar)
            %LOOKUPTABLE Construct an instance of Lookuptable class
            %   
            % obj = LOOKUPTABLE(ltvar)
            %
            %   Converts lookuptable structure object into a class object.
            %   Modifies double and dimer input properties to be of
            %   Filament class if applicable.
            % 
            %   Inputs:
            %       ltvar : (struct) contains single, double, and/or dimer
            %       fields full of N struct objects to be added to the
            %       lookuptable; struct should be created from
            %       makeLookupMat and reading in polymer-c outputs
            %
            %   Runs updateAddedNs and updateAddedStats.
            % 
            % See also LOOKUPTABLE/UPDATEADDEDNS, LOOKUPTABLE/UPDATEADDEDSTATS, MAKELOOKUPMAT.
            % any fields in the input struct that are arrays must be 1x# or, only if there are different values for different filaments, 2xN
            if nargin>0
                if class(ltvar)~="struct"
                    error("Invalid argument at position 1. Value must be of type struct or be convertible to struct.")
                end
                fields=fieldnames(ltvar);
                if nargin == 1
                    if length(fields)>3
                        error("Too many fields in input lookuptable variable")
                    elseif isempty(fields)
                        error("No fields in input lookuptable variable")
                    else
                       for i=1:length(fields)
                           field=fields{i};
                           if ismember(field,{'single','double','dimer'})
                               obj.(field)=ltvar.(field);
                               if ismember(field,{'double','dimer'})
                                   fields2=fieldnames(obj.(field));
                                   for j=1:length(fields2)
                                       Nval=obj.(field).(fields2{j}).N(1);
                                       obj.(field).(fields2{j}).N=Nval;
                                       fields3=fieldnames(obj.(field).(fields2{j}));
                                       for k=1:length(fields3)
                                           value=obj.(field).(fields2{j}).(fields3{k});
                                           if size(value,1)==2
                                               a=value(1,:);
                                               b=value(2,:);
                                               obj.(field).(fields2{j}).(fields3{k})=Filament(a,b);
                                           elseif size(value,2)==2
                                               a=value(:,1);
                                               b=value(:,2);
                                               obj.(field).(fields2{j}).(fields3{k})=Filament(a,b);
                                           elseif size(value,2)==(2*Nval)+1
                                               a=value(:,2:(Nval+1));
                                               b=value(:,(Nval+2):end);
                                               obj.(field).(fields2{j}).(fields3{k})=Filament(a,b);
                                           end
                                       end
                                   end
                               end
                           else
                               error("Input lookuptable variable contains invalid field name: %s",field)
                           end
                       end
                    end
                    obj.updateAddedNs();
                    obj.updateAddedStats();
                end
            end
        end

        function value= get.max_length(obj)
            % compute the maximum N value present for each filament type,
            % returns a FilType object
            value=FilType(find_maxN("single"),find_maxN("double"),find_maxN("dimer"));
            function max=find_maxN(type)
                if isempty(obj.(type))
                    max=zeros(1,0);
                else
                    Ns=fieldnames(obj.(type));
                    max=0;
                    for i=1:length(Ns)
                        N=str2double(erase(Ns{i},"N"));
                        if N>max
                            max=N;
                        end
                    end
                end
            end
        end

        function obj = getmissingNs(obj)
            %GETMISSINGNS updates the obj.missingNs property to contain any
            %N values that are missing for each filament type (single,
            %double, dimer)
            %
            %   out = LOOKUPTABLE.GETMISSINGNS
            %
            %   Scans through every entry and adds the N values
            %   accordingly. An N value is considered missing if it is not
            %   present in one of the types but is less than the maximum N
            %   value for that type. 
            %
            % See also LOOKUPTABLE.
            addNs("single");
            addNs("double");
            addNs("dimer");
            function addNs(type)
                obj.missingNs(1).(type)=[];
                for i=1:obj.max_length.(type)
                    Nlabel=strcat("N",num2str(i));
                    if ~isfield(obj.(type),Nlabel)
                        obj.missingNs.(type)=[obj.missingNs.(type);Nlabel];
                    end
                end
            end
        end

        function updateNList(obj)
            %UPDATENLIST creates a FilType object for the NList property
            %that contains the N values present in single, double, and
            %dimer
            %
            %   out = LOOKUPTABLE.UPDATENLIST
            %
            %   Scans through every entry and adds the N values
            %   accordingly.
            %
            %   Runs getmissingNs.
            %
            % See also LOOKUPTABLE, FILTYPE, LOOKUPTABLE/GETMISSINGNS.
            obj.NList=FilType(find_N("single"),find_N("double"),find_N("dimer"));
            obj.NList.intersectratio=true;
            obj.getmissingNs;
            function Nlist=find_N(type)
                if isempty(obj.(type))
                    Nlist=zeros(1,0);
                else
                    Ns=fieldnames(obj.(type));
                    Nlist=zeros(1,length(Ns));
                    for i=1:length(Ns)
                        N=str2double(erase(Ns{i},"N"));
                        Nlist(i)=N;
                    end
                end
            end
        end

        function updateAddedStats(obj)
            %UPDATEADDEDSTATS populates the AddedStats and StatNames properties with
            %entries from single, double, and dimer.
            %
            %   out = LOOKUPTABLE.UPDATEADDEDSTATS
            %
            %   Does not modify data, but rearranges how the information in
            %   single, double, and dimer is stored/ accessed for ease of
            %   access. Now indexed by Lookuptable.(stat).(type).(N)
            %
            % See also LOOKUPTABLE.
            value=struct;
            addstats("single")
            addstats("double")
            addstats("dimer")

            obj.AddedStats= value;
            obj.StatNames=string(fieldnames(obj.AddedStats));
            function addstats(type)
                if ~isempty(obj.(type))
                    Ns=fieldnames(obj.(type));
                    for i=1:length(Ns)
                        vars=fieldnames(obj.(type).(Ns{i}));
                        for j=1:length(vars)
                            var=vars{j};
                            if ~isfield(value,var)
                                value.(var)=FilType();
                            end
                            value.(var).(type).(Ns{i})=obj.(type).(Ns{i}).(var);
                        end
                    end
                end
            end
        end

        function updateAddedNs(obj)
            %UPDATEADDEDNS populates the AddedNs and NNames properties with
            %entries from single, double, and dimer.
            %
            %   out = LOOKUPTABLE.UPDATEADDEDNS
            %       
            %   Runs updateNList.
            %
            %   Does not modify data, but rearranges how the information in
            %   single, double, and dimer is stored/ accessed for ease of
            %   access. Now indexed by Lookuptable.(N).(type).(stat)
            %
            % See also LOOKUPTABLE, LOOKUPTABLE/UPDATENLIST.
            value=struct;
            addNs("single")
            addNs("double")
            addNs("dimer")

            obj.AddedNs=value;
            obj.NNames=string(fieldnames(obj.AddedNs));
            obj.updateNList();
            function addNs(type)
                if length(obj.(type))>0
                    Ns=fieldnames(obj.(type));
                    for i=1:length(Ns)
                        N=Ns{i};
                        if ~isfield(value,N)
                            value.(N)=FilType();
                        end
                        value.(N).(type)=obj.(type).(N);
                    end
                end
            end
        end
        
        function addN(obj,Nstruct)
            %ADDN adds the input N struct, overwrites any existing struct
            %for that N and type entry
            %
            %   out = LOOKUPTABLE.ADDN(Nstruct)
            %
            %   Inputs:
            %       Nstruct : (struct) a singular N struct with
            %       corresponding polymer stat values. Must have a 'type'
            %       and 'N' property.
            %       
            %   Runs updateAddedNs and updateAddedStats.
            %
            %   For any entries with 2 rows/ columns, splits them into
            %   filament objects
            %
            % See also LOOKUPTABLE, LOOKUPTABLE/UPDATEADDEDNS, LOOKUPTABLE/UPDATEADDEDSTATS.
            Ntype=Nstruct.type;
            N=num2str(Nstruct.N(1));
            obj.(Ntype).(strcat("N",N))=Nstruct;
            if ismember(Ntype,{'double','dimer'})
               fields2=strcat("N",N);
               obj.(Ntype).(fields2).N=str2double(N);
               fields3=fieldnames(obj.(Ntype).(fields2));
               for k=1:length(fields3)
                   value=obj.(Ntype).(fields2).(fields3{k});
                   if size(value,1)==2
                       a=value(1,:);
                       b=value(2,:);
                       obj.(Ntype).(fields2).(fields3{k})=Filament(a,b);
                   elseif size(value,2)==2
                       a=value(:,1);
                       b=value(:,2);
                       obj.(Ntype).(fields2).(fields3{k})=Filament(a,b);
                   elseif size(value,2)==(2*(Nstruct.N(1)))+1
                       a=value(:,2:(Nstruct.N(1)+1));
                       b=value(:,(Nstruct.N(1)+2):end);
                       obj.(Ntype).(fields2).(fields3{k})=Filament(a,b);
                   end
               end
            end
            obj.updateAddedNs();
            obj.updateAddedStats();
        end

        function delN(obj,N,type)
            %DELN removes the entry for the specified N value (FH1 length) and type
            %
            %   out = LOOKUPTABLE.DELN(N,type)
            %
            %   Inputs:
            %       N     : (double) N value to delete
            %       type  : (String) type to delete (single, double, or dimer)
            %       
            %   Runs updateAddedNs and updateAddedStats.
            %
            % See also LOOKUPTABLE, LOOKUPTABLE/UPDATEADDEDNS, LOOKUPTABLE/UPDATEADDEDSTATS.
            obj.(type)=rmfield(obj.(type),(strcat("N",num2str(N))));
            obj.updateAddedNs();
            obj.updateAddedStats();
        end

        function out=copytable(obj)
            %COPYTABLE create new Lookuptable with the same properties but
            %not the same object (since this is a handle class)
            %
            %   out = LOOKUPTABLE.COPYTABLE
            %
            %   Copies the following properties: holdratio, ratio, single,
            %   double, dimer, StatNames, NNames, AddedStats, AddedNs,
            %   NList, interpolant, missingNs
            %
            % See also LOOKUPTABLE.
            out=Lookuptable();
            out.holdratio=obj.holdratio;
            out.ratio=obj.ratio;
            out.single=obj.single;
            out.double=obj.double;
            out.dimer=obj.dimer;
            out.StatNames=obj.StatNames;
            out.NNames=obj.NNames;
            out.AddedStats=obj.AddedStats;
            out.AddedNs=obj.AddedNs;
            out.NList=obj.NList;
            out.interpolant=obj.interpolant;
            out.missingNs=obj.missingNs;
        end

        function out=stattable(obj,stat,type)
            %STATTABLE create matrix of lookuptable stats for a specific type 
            %
            %   out = LOOKUPTABLE.STATTABLE(stat,type)
            %
            %   Inputs:
            %       stat  : (String) lookup table stat to retrieve
            %       type  : (String) single, double, dimer, or ratio
            %       
            %   Output is table with 3 columns: [FH1 size, PRM location,
            %   stat value]
            %
            %   Uses getstat, AddedStats, and NList
            %
            % See also LOOKUPTABLE, LOOKUPTABLE/GETSTAT, LOOKUPTABLE/UPDATEADDEDSTATS.
            obj.holdratio=true;
            stats=struct2cell(obj.AddedStats.(stat).(type));
            if class(stats{1})=="Filament"
                mata=makemat(struct2cell(obj.getstat(Fil="a",Stat=stat,type=type)));
                matb=makemat(struct2cell(obj.getstat(Fil="b",Stat=stat,type=type)));
                out=Filament(mata,matb);
            else
                out=makemat(stats);
            end

            obj.holdratio=false;
            
            function mat=makemat(stats2)
                if type=="ratio"
                    numrows=min(obj.max_length.dimer,obj.max_length.double)*size(stats2,1);
                else
                    numrows=obj.max_length.(type)*size(stats2,1);
                end
                mat=zeros(numrows,3);
                for i=1:size(stats2,1)
                    N=obj.NList.(type)(i);
                    if length(stats2{i})>N
                        istart=(N+1)*N/2;
                        looplen=length(stats2{i});
                    else
                        istart=(N-1)*N/2;
                        looplen=N;
                    end
                    
                    for j=1:looplen
                        mat(istart+j,1)=N;
                        mat(istart+j,2)=j;
                        if length(stats2{i})==1
                            mat(istart+j,3)=stats2{i}(1);
                        else
                            mat(istart+j,3)=stats2{i}(j);
                        end
                    end
                end
                mat=mat(any(mat,2),:);
            end
        end

        function value = get.ratio(obj)
            % compute (or retrive) element by element ratio dimer/double
            if isempty(obj.ratiostore)
                value=getratio(obj);
                obj.ratiostore=value;
            elseif obj.holdratio
                value=obj.ratiostore; % use the stored ratio value (instead of re calculating) if holdratio is true
            else
                value=getratio(obj);
                obj.ratiostore=value; 
            end
            function out= getratio(input)
                % compute element by element ratio dimer/double
                if isempty(input.dimer)
                    out=[];
                elseif class(input.dimer)=="struct"
                    fieldsdimer=fieldnames(input.dimer);
                    fieldsdouble=fieldnames(input.double);
                    fields=intersect(fieldsdimer,fieldsdouble);
                    out=struct;
                    for i=1:length(fields)
                        temp1=FilType(0,input.double.(fields{i}),input.dimer.(fields{i}));
                        temp=getratio(temp1);                                                                                                                                                                                                                                
                        out.(fields{i})=temp;
                    end
                elseif class(input.dimer)=="Filament"
                    outa=getratio(FilType(0,input.double.a,input.dimer.a));
                    outb=getratio(FilType(0,input.double.b,input.dimer.b));
                    out=Filament(outa,outb);
                elseif class(input.dimer)=="string"
                    out="n/a";
                else
                    if input.double==0
                        out=NaN;
                    elseif input.dimer==0
                        out=0;
                    else
                        out=input.dimer./input.double;
                    end
                end
            end
        end
    end

    methods
      output = getstat(obj,NameValueArgs)
    end

    methods
      output=extrapolate(obj,stat)
    end


end