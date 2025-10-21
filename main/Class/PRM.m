classdef PRM < handle & dynamicprops
%PRM Contains information and functions about one specific PRM in a formin.
    %
    %   Construction:
    %       obj = PRM(formin, dist_FH2, dist_NT, size, dist_FH2_start)
    %
    %   See also FORMIN, OPTIONS.

    properties(SetAccess={?Formin})
        formin Formin % parent formin for the PRM
        dist_FH2 double % distance from center of PRM to FH2
        dist_NT double % distance from center of PRM to NT
        size double % number of prolines
        dist_FH2_start double % distance from CT P of PRM to FH2
        stat_props struct=struct % contains meta.DynamicProperty objects, from lookuptable
    end

    properties(Dependent)
        gating double % gating factor
        c_PA double % concentration of profilin-actin | μM
        kpoly FilType % rate of polymerization by FilType
        fh1length double % length of parent FH1 domain
        FH2dist_frac double % fractional distance from center of PRM to FH2
        
        % Parameters all by filament and FilType:
        kcap FilType % rate of PRM + profilin-actin binding (capture)| s^(-1)
        kdel FilType % rate of barbed end + PRM-profilin-actin binding (delivery) | s^(-1)
        rcap FilType % rate of PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        rdel FilType % rate of barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        krel FilType % rate of PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
    end

    methods
        function obj = PRM(formin, dist_FH2, dist_NT, size, dist_FH2_start)
            %PRM Construct an instance of PRM
            %  Inputs:
            %       formin  : parent formin
            %       dist_FH2: distance from the PRM to the FH2
            %       dist_NT : distance from the PRM to the N-terminus
            %       size    : number of amino acid in the PRM
            %       dist_FH2_start: distance from the C-terminal proline of
            %                       the PRM to the FH2
            if nargin>0
                obj.formin=formin;
                obj.dist_FH2=dist_FH2;
                obj.dist_NT=dist_NT;
                obj.size=size;
                obj.dist_FH2_start=dist_FH2_start;

                statlist=obj.formin.lookup.StatNames;
                for i=1:length(statlist)
                    obj.addStat(statlist(i));
                end
                addlistener(obj.formin.lookup,'StatNames','PostSet',@obj.updateStats);
            end
        end

        function value=get.gating(obj)
            % obtain gating factor from parent formin
            value=obj.formin.gating;
        end

        function value=get.c_PA(obj)
            % obtain concentration of profilin-actin from parent formin
            value=obj.formin.c_PA;
        end

        function value=get.kcap(obj)
            % compute the rate of capture, using k_cap and equation from options 
            value=obj.formin.opts.k_cap*obj.formin.opts.equations.kcap(obj);
        end

        function value=get.kdel(obj)
            % compute the rate of delivery, using k_del and equation from options 
            value=obj.formin.opts.k_del*obj.formin.opts.equations.kdel(obj);
        end

        function value=get.rcap(obj)
            % compute the rate of reverse capture, using r_cap and equation from options 
            value=obj.formin.opts.r_cap*obj.formin.opts.equations.rcap(obj);
        end

        function value=get.rdel(obj)
            % compute the rate of reverse delivery, using r_del and equation from options 
            value=obj.formin.opts.r_del*obj.formin.opts.equations.rdel(obj);
        end

        function value=get.krel(obj)
            % compute the rate of release, using k_rel and equation from options 
            value=obj.formin.opts.k_rel*obj.formin.opts.equations.krel(obj);
        end

        function value=get.kpoly(obj)
            %compute the overall polymerization rate, integrating the
            % rates of each of the steps with kpolymerization
            kpoly_type = obj.formin.opts.kpoly_type;
            value=kpolymerization(kpoly_type,obj.kcap,obj.kdel,obj.rcap,obj.rdel,obj.krel);
        end

        function r=calculate_kpoly(obj,kpoly_type)
            % CALCULATE_KPOLY calculate kpoly for a specific kpoly type
            % (i.e. 4st vs. 3st)
            %
            %   kpoly= PRM.CALCULATE_KPOLY(type) calculates kpoly for the
            %   kpoly type type
            %
            %   kpoly= PRM.CALCULATE_KPOLY calculates kpoly for the
            %   kpoly type specified in the parent formin's opts
            %
            % See also KPOLYMERIZATION.
            arguments
                obj PRM
                kpoly_type string = obj.formin.opts.kpoly_type
            end
            r=kpolymerization(kpoly_type,obj.kcap,obj.kdel,obj.rcap,obj.rdel,obj.krel);
        end

        function value=get.fh1length(obj)
            %  obtain formin length from parent formin
            value=obj.formin.length;
        end

        function value=get.FH2dist_frac(obj)
            % compute fractional distance of PRM from FH2 wrt total formin
            % length
            value=obj.dist_FH2./obj.formin.length;
        end

        function obj = addStat(obj,stat)
            % ADDSTAT add a polymer statistic to the object from the
            % lookuptable
            %
            %   obj= PRM.ADDSTAT(stat) adds the stat from the parent 
            %   formin's lookuptable entry corresponding to the PRM's location
            %
            %   Input:
            %       stat : (string) the name of the polymer statistic to
            %       add, must be a property of the relevent lookuptable
            %
            %   Stats are added as dependent properties and an entry in
            %   obj.stat_props is updated
            %
            % See also FORMIN, LOOKUPTABLE/GETSTAT.
            if ~isprop(obj,stat)
                prop=obj.addprop(stat);
                prop.Dependent=true;
                prop.GetMethod = @(obj) obj.formin.lookup.getstat(iSite=obj.dist_FH2,N=obj.formin.length,Stat=stat,NName=obj.formin.NName);
                obj.stat_props.(stat)=prop;
            end
            obj.addBaseStat(stat);
        end

        function obj = addBaseStat(obj,stat)
            % ADDBASESTAT add a polymer statistic to the object from the
            % lookuptable, for the amino acid at the FH2 domain
            %
            %   obj= PRM.ADDBASESTAT(stat) adds the stat from the parent 
            %   formin's lookuptable entry corresponding to location 0 of
            %   the corresponding FH1 length
            %
            %   Input:
            %       stat : (string) the name of the polymer statistic to
            %       add, must be a property of the relevent lookuptable
            %
            %   Stats are added as dependent properties and an entry in
            %   obj.stat_props is updated
            %
            % See also FORMIN, LOOKUPTABLE/GETSTAT.
            if ~isprop(obj,strcat(stat,"_Base"))
                prop=obj.addprop(strcat(stat,"_Base"));
                prop.Dependent=true;
                prop.GetMethod = @(obj) obj.formin.lookup.getstat(iSite=1,N=obj.formin.length,Stat=stat,NName=obj.formin.NName);
                obj.stat_props.(strcat(stat,"_Base"))=prop;
            end
        end

        function obj = updateStats(obj,src,evnt)
            % UPDATESTATS refresh polymer statistics for the PRM
            %
            %   obj= PRM.UPDATESTATS updates the polymer statistics for the
            %   PRM
            %
            %
            %   All values in obj.stat_props are deleted and new stats are
            %   added via PRM.addStat for each stat in the parent formin's
            %   StatNames property.
            %
            % See also PRM/ADDSTAT.
           fields=fieldnames(obj.stat_props);
           for i=1:length(fields)
               delete(obj.stat_props.(fields{i}))
           end
           obj.stat_props=struct;
           statlist=obj.formin.lookup.StatNames;
           for i=1:length(statlist)
               obj.addStat(statlist(i));
           end
        end

        function names=getfields(obj)
            % GETFIELDS obtain list of fieldnames for the PRM or array of
            % PRMs
            %
            %   names= PRM.GETFIELDS obtain list of fieldnames for the PRM or array of
            %   PRMs
            %   
            %   Input PRM can be a single PRM object or an array of PRMs.
            %   If an input array, the output is an array of fieldnames.
            %
            %
            %   Output fieldnames are strings.
            if length(obj)==1
                names=fieldnames(obj);
            else
                fields1=fieldnames(obj);
                fields2=fieldnames([obj.stat_props]);
                names=[fields1;fields2];
            end
        end

        function out=getprop(obj,prop)
            % GETPROP obtain property of PRM or array of PRMs
            %
            %   out= PRM.GETPROP(prop) obtain property prop for the PRM or array of
            %   PRMs
            %   
            %   Input PRM can be a single PRM object or an array of PRMs.
            %   If an input array, the output is an array of property values.

            if length(obj)==1
                out=obj.(prop);
            else
                out=[];
                for i=1:length(obj)
                    out=[out,obj(i).(prop)];
                end
            end
        end


    end

    methods (Static)
        function obj = loadobj(s)
            %LOADOBJ run tasks on load of a PRM object
            %
            % On load, run PRM.updateStats and add new listeners for the
            % parent formin's corresponding lookuptable.
            %
            % See also PRM/UPDATESTATS.
         obj=s;
         obj.updateStats;
         addlistener(obj.formin.lookup,'StatNames','PostSet',@obj.updateStats);
      end
    end
    

end