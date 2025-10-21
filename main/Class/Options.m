classdef Options <handle
    %OPTIONS Contains options for reading in formins, calculating rates,
    %including all of the rate constants
    %
    %   Construction:
    %       obj = OPTIONS(lookup,python_path,kpoly_type,resultsdir,k_cap,k_del,r_cap,r_del,k_rel)
    %           lookup      : (Lookuptable) refrence lookup table
    %                          object
    %           python_path : (String) path to python files
    %           kpoly_type  : (String) kpoly model type (capture, 3st,
    %                          or 4st)
    %           resultsdir  : (String) where to save results
    %           k_cap       : (double) capture rate constant
    %           k_del       : (double) delivery rate constant
    %           r_cap       : (double) reverse capture rate constant
    %           r_del       : (double) revsere delivery rate constant
    %           k_rel       : (double) release rate constant
    %           r_cap_exp   : (double) constant in front of PRM size in r_cap (if applicable, default is 1)
    % 
    % 
    %   See also FORMIN.

    properties (SetObservable, AbortSet)
        python_path string % location of python files
        kpoly_type string {mustBeMember(kpoly_type,{'capture','3st','4st'})}='3st' % model for kpoly
        equations struct % equations used to calculate rates
        lookup Lookuptable % lookup table with polymer stats

        % reading sequence options (see get_formin_info.py):
        min_lenPRM double % minimum PRM length (including interruptions)
        nInt double % number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions)
        max_lenInt double % maximum interruption length
        min_nP double % minimum number of Ps 
        NTopt double % FH1 NT definition options 
        CTopt double % FH1 CT definition options
        PRMscanopt double % PRM scan options

        colors (:,1) % list of hexcodes to be used for plotting
        shapes (:,1) % list of shape symbols to be used for plotting
        resultsdir % location of folder to save results folder to
        resultsfolder % results folder name

        equationstext struct % string for each equation 

        notes string % notes to add to output

        % parameters
        k_cap double % rate constant for PRM + profilin-actin binding (capture) | μM^(-1)s^(-1)
        k_del double % rate constant for barbed end + PRM-profilin-actin binding (delivery) | μM^(-1)s^(-1)
        r_cap double % rate constant for PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        r_del double % rate constant for barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        k_rel double % rate constant for PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
        r_cap_exp double % constant in front of PRM size in r_cap (if applicable)

        % delivery location, only used if added to kpoly calcs
        del_x double=0 % x location (parallel to FH2) for delivery, only used if in kpoly calcs
        del_y double=0 % y location (perpendicular to FH2) for delivery, only used if in kpoly calcs
    end

    methods
        function obj = Options(lookup,python_path,kpoly_type,resultsdir,k_cap,k_del,r_cap,r_del,k_rel,r_cap_exp)
            %OPTIONS Construct an instance of options class
            %    
            %   obj = OPTIONS(lookup,python_path,kpoly_type,resultsdir,k_cap,k_del,r_cap,r_del,k_rel)
            % 
            %   Inputs:
            %       lookup      : (Lookuptable) refrence lookup table
            %                      object
            %       python_path : (String) path to python files
            %       kpoly_type  : (String) kpoly model type (capture, 3st,
            %                     or 4st)
            %       resultsdir  : (String) where to save results
            %       k_cap       : (double) capture rate constant
            %       k_del       : (double) delivery rate constant
            %       r_cap       : (double) reverse capture rate constant
            %       r_del       : (double) revsere delivery rate constant
            %       k_rel       : (double) release rate constant
            %       r_cap_exp   : (double) constant in front of PRM size in r_cap (if applicable, default is 1)
            % 
            %   Sets input variables, update_results_folder, and set_fh1
            % 
            % See also OPTIONS/SET_FH1, LOOKUPTABLE, OPTIONS/UPDATE_RESULTS_FOLDER.
            arguments
              lookup Lookuptable
              python_path string
              kpoly_type string
              resultsdir string
              k_cap double
              k_del double
              r_cap double
              r_del double
              k_rel double
              r_cap_exp double =1
            end
            if nargin>0
            obj.lookup=lookup;
            obj.python_path=python_path;
            obj.kpoly_type=kpoly_type;
            obj.resultsdir=resultsdir;
            obj.update_results_folder;
            obj.k_cap = k_cap;
            obj.k_del = k_del;
            obj.r_cap = r_cap;
            obj.r_del = r_del;
            obj.k_rel = k_rel;
            obj.r_cap_exp=r_cap_exp;
            obj.set_FH1();
            end
        
        end

        function update_results_folder(obj)
            %UPDATE_RESULTS_FOLDER update results folder to current time
            % 
            % OPTIONS.UPDATE_RESULTS_FOLDER
            % 
            %   Format is "RESULTS_yyyy-MM-dd HH-mm"
            % 
            % See also OPTIONS.m.
            time= datetime('now', 'Format','yyyy-MM-dd HH-mm');
            time= string(time);
            obj.resultsfolder= 'RESULTS_' + time;
        end

        function out=getconsts(obj)
            %GETCONSTS get array of rate constants
            %
            %   OPTIONS.GETCONSTS
            %
            %   output is this array: [k_cap,k_del,r_cap,r_del,k_rel,r_cap_exp]
            %
            % See also OPTIONS.m.
            out=[obj.k_cap,obj.k_del,obj.r_cap,obj.r_del,obj.k_rel,obj.r_cap_exp];
        end

        function set_FH1(obj,NameValueArgs)
            % SET_FH1 set options for reading in FH1 sequences
            %
            %   OPTIONS.SET_FH1 sets FH1 options properties to the default
            %   values 
            % 
            %   OPTIONS.SET_FH1(NameValueArgs) sets FH1 options properties
            %   to the values specified plus other default values
            %   
            %   NameValueArgs:
            %       min_lenPRM : (double) minimum PRM length (including interruptions), use instead of opts property (default is 4)
            %       nInt       : (double) number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions), use instead of opts property (default is 1)
            %       max_lenInt : (double) maximum interruption length, use instead of opts property (default is 1)
            %       min_nP     : (double) minimum number of Ps in PRM, use instead of opts property (default is 4)
            %       NTopt      : (double) FH1 NT definition options , use instead of opts property (default is 1)
            %                   1 - first instance of PRM with at least 4 Ps with max 1 interruption of length 1 (in sequence of at least 3 PRMs no father than 100 amino acids apart) (default)
            %                   2 - first instance of PRM (as defined by args 3 and 4) (uin sequence of at least 3 PRMs no father than 100 amino acids apart)
            %                   3 - Uniprot defined FH1 start (or option 1 if no FH1)
            %                   4 - Uniprot defined FH1 start (or option 2 if no FH1)
            %                   5 - Start of sequence (for input sequences)
            %       CTopt      : (double) FH1 CT definition options, use instead of opts property (default is 1)
            %                   1 - Uniprot defined FH2 start (default)
            %                   2 - Uniprot defined FH1 end (or option 1 if no FH1)
            %                   3 - End of sequence (for input sequences)
            %       PRMscanopt : (double) PRM scan options, use instead of opts property (default is 1)
            %                   1 - Search for PRMs starting from the FH2 (CT) (default)
            %                   2 - Search for PRMs starting from the FH1 NT
            %
            % See also OPTIONS.m, FORMIN.
            arguments
                obj Options
                NameValueArgs.min_lenPRM double=4 % minimum PRM length (including interruptions)
                NameValueArgs.nInt double=1 % number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions)
                NameValueArgs.max_lenInt double=1 % maximum interruption length
                NameValueArgs.min_nP double=4 % minimum number of Ps 
                NameValueArgs.NTopt double {mustBeInRange(NameValueArgs.NTopt,1,5)}=1 % FH1 NT definition options 
                NameValueArgs.CTopt double {mustBeInRange(NameValueArgs.CTopt,1,3)}=1 % FH1 CT definition options
                NameValueArgs.PRMscanopt double {mustBeInRange(NameValueArgs.PRMscanopt,1,2)}=1 % PRM scan options
            end
            obj.min_lenPRM=NameValueArgs.min_lenPRM;
            obj.nInt=NameValueArgs.nInt; 
            obj.max_lenInt=NameValueArgs.max_lenInt; 
            obj.min_nP=NameValueArgs.min_nP; 
            obj.NTopt=NameValueArgs.NTopt;
            obj.CTopt=NameValueArgs.CTopt;
            obj.PRMscanopt=NameValueArgs.PRMscanopt;
        end

        function set_equation(obj,preset,step,vars)
            %SET_EQUATION set equation(s) for rate calculations
            % 
            % OPTIONS.SET_EQUATION(preset) set equations based on preset
            %   value
            %       preset 1 : 
            %           kcap= (1-POcclude)*c_PA
            %           kdel= (1-POcclude_base)*(1.0e33*Prvec0/27*6.022e23)*gating
            %           rcap= e^-size*r_cap_exp
            %       preset 2 : 
            %           kcap= (1-POcclude)*c_PA
            %           kdel= (1-POcclude_base)*gating
            %           rcap= e^-size*r_cap_exp
            % 
            % OPTIONS.SET_EQUATION(0,step,{var1,type1,var2,type2,...}) set equation for
            %   step as equal to the product of each of the input vars
            %   implemented accoridng to the following types. Can have as
            %   many repeating var,type as you want. 
            %
            % OPTIONS.SET_EQUATION(0,step1,{var1,type1,var2,type2,...},step2, {var,type...},...) 
            %   set equations for each step as equal to the product of each 
            %   of the input vars implemented accoridng to the following types. 
            %   Can have as many repeating var,type as you want. Can have as many
            %   repeating step,{var,type,...} as you want.
            % 
            % Note: all equations are still multiplied by the corresponding
            % rate constant
            %
            % If a variable is set for r_cap using negexp, the resulting
            % equation is e^-var*r_cap_exp
            % 
            % Inputs:
            %       preset : (double) preset to use (1,2,3) or, if 0, don't
            %               use a preset and instead use input values
            %       step   : (String) step to set equation to (repeating)
            %               options- 'kcap','kdel','rcap','rdel','krel'
            %       vars   : (cell) of format {var1,type1,var2,type2,...}
            %               indicating the variables and how to implement 
            %               them into the equations
            %           var must be a lookuptable property or 'dist_FH2','dist_NT', 'size','gating','c_PA','dist_FH2_start'
            %           type must be one of the following:
            %               "linear"    : var
            %               "negexp"    : e^-var (or e^-var*r_cap_exp if the step is r_cap)
            %               "exp"       : e^var
            %               "1-"        : (1-var)
            %               "amino"     : (1.0e33*var/27*6.022e23)
            %               "1-base"    : (1-var_base) (var at PRM loc 0)
            %
            %   Sets equations and equationstext for all of the inputs, and
            %   sets equations equal to 1 for any rate equation not yet
            %   specified (if an eqaution already exists but is not
            %   specified, no changes are made)
            %
            %
            % See also OPTIONS.m, LOOKUPTABLE.
            arguments
                obj
                preset double {mustBeMember(preset,[0,1,2,3])}
            end
            arguments (Repeating)
                step string {mustBeMember(step,{'kcap','kdel','rcap','rdel','krel'})}
                vars (1,:) cell 
            end
            
            if isempty(obj.equations)
                obj.equations=struct;
            end
            if isempty(obj.equationstext)
                obj.equationstext=struct;
            end
            if preset==0
                vars=vars{:};
                if logical(mod(length(vars),2))
                    if vars{end}==""
                        vars=vars(1:end-1);
                    else
                        eid = 'validvars:sizeNotEven';
                        msg = 'An even number of inputs must be provided to vars.';
                        throwAsCaller(MException(eid,msg))
                    end
                end
                valid_vars=obj.lookup.StatNames;
                scale_types={'linear','negexp','exp','1-','amino','1-base'};
                for k=1:2:length(vars)
                    if ~ismember(vars{k},valid_vars) && ~ismember(vars{k},{'dist_FH2','dist_NT', 'size','gating','c_PA','dist_FH2_start'})
                        eid = 'validvars:varNotValid';
                        msg = 'Variables to include in equations must either be in the reference lookuptable or "dist_FH2","dist_NT","size".';
                        throwAsCaller(MException(eid,msg))
                    end
                    if ~ismember(vars{k+1},scale_types)
                        valid_scale_string="";
                        for j=1:length(scale_types)
                            if j==length(scale_types)
                                valid_scale_string=strcat(valid_scale_string,scale_types{k});
                            else
                                valid_scale_string=strcat(valid_scale_string,scale_types{k},", ");
                            end
                        end
    
                        eid = 'validvars:scaleTypeNotValid';
                        msg = strcat('Sacle type must be one of the following: ',valid_scale_string,".");
                        throwAsCaller(MException(eid,msg))
                    end
                end
                invars=vars;
                if invars{end}==""
                    invars=invars(1:end-1);
                end
                obj.equations.(step{:})=makeeq();
                obj.equationstext.(strcat(step{:},"_eq"))=[sprintf("%s(%s),",invars{1:end})];
            elseif preset==1
                invars={"POcclude","1-","c_PA","linear"};
                obj.equations.kcap=makeeq();
                obj.equationstext.kcap_eq=[sprintf("%s(%s),",invars{1:end})];
                invars={"POcclude","1-base","Prvec0","amino","gating","linear"};
                obj.equations.kdel=makeeq();
                obj.equationstext.kdel_eq=[sprintf("%s(%s),",invars{1:end})];
                invars={"size","negexp"};
                step="rcap";
                obj.equations.rcap=makeeq();
                obj.equationstext.rcap_eq=[sprintf("%s(%s),",invars{1:end})];
            elseif preset==2
                invars={"POcclude","1-","c_PA","linear"};
                obj.equations.kcap=makeeq();
                obj.equationstext.kcap_eq=[sprintf("%s(%s),",invars{1:end})];
                invars={"POcclude","1-base","gating","linear"};
                obj.equations.kdel=makeeq();
                obj.equationstext.kdel_eq=[sprintf("%s(%s),",invars{1:end})];
                invars={"size","negexp"};
                step="rcap";
                obj.equations.rcap=makeeq();
                obj.equationstext.rcap_eq=[sprintf("%s(%s),",invars{1:end})];
            end
            
            steps={"kcap","kdel","rcap","rdel","krel"};
            for j=1:length(steps)
                if ~isfield(obj.equations,steps{j})
                    obj.equations.(steps{j})=@(PRM) 1;
                end
            end
            
            function fxn = makeeq()
                fxn=@(PRM) 1;
                for i=1:2:length(invars)
                    if invars{i+1}=="linear"
                        fxn=@(PRM) fxn(PRM)*PRM.(invars{i});
                    elseif invars{i+1}=="negexp"
                        if step=="rcap"
                            fxn=@(PRM) fxn(PRM)*exp(-1*PRM.(invars{i})*PRM.formin.opts.r_cap_exp);
                        else
                            fxn=@(PRM) fxn(PRM)*exp(-1*PRM.(invars{i}));
                        end
                    elseif invars{i+1}=="exp"
                        fxn=@(PRM) fxn(PRM)*exp(PRM.(invars{i}));
                    elseif invars{i+1}=="1-"
                        fxn=@(PRM) fxn(PRM)*(1-PRM.(invars{i}));
                    elseif invars{i+1}=="amino"
                        fxn=@(PRM) fxn(PRM)*(1.0e33*(PRM.(invars{i}))/(27*6.022e23));
                    elseif invars{i+1}=="1-base"
                        fxn=@(PRM) fxn(PRM)*(1-PRM.(strcat(invars{i},"_Base")));
                    end
                end
            end
        end

        function cleareq(obj)
            %CLEAREQ remove all rate equations and replace them with 1
            %
            %   OPTIONS.CLEAREQ
            %
            %   Updates obj.equations and removed all fields from obj.equationstext
            %
            % See also OPTIONS.m.
            eqfields=fieldnames(obj.equationstext);
            for i=1:length(eqfields)
                obj.equationstext=rmfield(obj.equationstext,eqfields{i});
            end
            eqfields=fieldnames(obj.equations);
            for i=1:length(eqfields)
                obj.equations.(eqfields{i})=@(PRM) 1;
            end
        end

        function tab=optionstable(obj)
            % OPTIONSTABLE create a table with properties of the options
            % object
            %
            %   tab = OPTIONS.OPTIONSTABLE
            %
            %   Table entry for every property that is a double or string
            %   (of length 1) as well as every entry in obj.equationstext
            %
            % See also OPTIONS.m.
            tabstruct=struct;

            addel(obj);
            addel(obj.equationstext);

            Rowtitles=fieldnames(tabstruct);
            Datavars=struct2cell(tabstruct);
            tab=cell2table([Datavars]);
            tab.Properties.RowNames=Rowtitles;
            tab.Properties.VariableNames={' '};

            function addel(object)
                optfields=fieldnames(object);
                for i=1:length(optfields)
                    if class(object.(optfields{i}))=="double" 
                        tabstruct.(optfields{i})=num2str(object.(optfields{i}));
                    elseif class(object.(optfields{i}))=="string" && length(object.(optfields{i}))==1
                        tabstruct.(optfields{i})=object.(optfields{i});
                    end
                end
            end
        end

        function applytable(obj,row)
            % APPLYTABLE modify object properties to allign with those
            % listed in a table
            %
            %   OPTIONS.APPLYTABLE(tab) modify object properties to allign
            %   with those in the table tab
            %
            %   Searches for table variables with the same name as an
            %   Options property and applies them.
            %   Searches for table variables corresponding to any of the
            %   rate equations and runs obj.set_equation accordingly. These
            %   entries are assumed to have the name "rate_eq" and be of
            %   the format of obj.equationstext entries.
            %
            % See also OPTIONS.m, OPTIONS/OPTIONSTABLE, OPTIONS/SET_EQUATION.
            arguments
                obj Options
                row table
            end
            vars=row.Properties.VariableNames;
            obj.cleareq;
            for i=1:length(vars)
                value=row.(vars{i});
                if ismember(vars{i},fieldnames(obj))
                    if class(obj.(vars{i}))=="double" && class(value)=="string"
                        value=str2double(value);
                    end
                    if obj.(vars{i})~=value
                        obj.(vars{i})=value;
                    end
                    continue
                end
                if ismember(vars{i},strcat([fieldnames(obj.equations)],"_eq"))
                    obj.set_equation(0,vars{i}(1:end-3),cellfun(@string,strsplit(char(value),{',','(',')'}),"UniformOutput",false))
                    continue
                end
            end
        end
    end
end