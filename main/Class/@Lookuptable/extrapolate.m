function output=extrapolate(obj,stat)
% EXTRAPOLATE generate extrapolation function for the specified polymer
% statistic
% 
% output = LOOKUPTABLE.EXTRAPOLATE(stat)
%
%   Inputs: 
%       stat : (string) lookuptable property to extrapolate
%
%   Output is a FilType with an extrapolation function for each filament
%   entry (as well as filament a and b if dimer or double)
%
%   If the lookuptable property is a string, the functions will be a string
%   that says "unable to extrapolate," unless the property is 'type,' in
%   which case the function will return the type.
%
%   If the stat is dependent on PRM location, uses scatteredInterpolant
%   ('linear','nearest' setting)and output functions can be accessed with
%   fxn(FH1 size, PRM loc). 
%   
%   If the stat is the same for the same FH1 size, uses polyfit, and the
%   output functions can be accessed with fxn(FH1 size).
% 
%   Adds the output to obj.interpolant. If obj.interpolant.(stat) already
%   exists, the output will be the preexisting value and new calculations
%   will not be done. 
%
% See also LOOKUPTABLE, POLYFIT, SCATTEREDINTERPOLANT.
    tempNList=obj.NList;
    if isfield(obj.interpolant,stat)
        output=obj.interpolant.(stat);
    else
        statobj=obj.AddedStats.(stat);
        if stat=="type"
            singlefxn= @(n) "single";
            doubleafxn= @(n) "double";
            dimerafxn= @(n) "dimer";
            output=FilType(singlefxn,doubleafxn,dimerafxn);
        elseif class(statobj.single.(obj.NNames(1)))=="string"
            singlefxn= @(n) "unable to extrapolate";
            doubleafxn= @(n) "unable to extrapolate";
            dimerafxn= @(n) "unable to extrapolate";
            output=FilType(singlefxn,doubleafxn,dimerafxn);
        else
            singlefxn=getfxn("single",1);
            doubleafxn=getfxn("double",1);
            dimerafxn=getfxn("dimer",1);
            if class(statobj.dimer.(strcat("N",num2str(tempNList.dimer(1)))))=="Filament"
                doublebfxn=getfxn("double",2);
                dimerbfxn=getfxn("dimer",2);
                output=FilType(singlefxn,doubleafxn,dimerafxn,doublebfxn,dimerbfxn);
            else
                output=FilType(singlefxn,doubleafxn,dimerafxn);
            end
        end
        
        if isempty(obj.interpolant)
            x.(stat)=output;
            obj.interpolant=x;
        else
            obj.interpolant.(stat)=output;
        end
    end
    
    function fxn=getfxn(type,fil)
        if fil==1
            filname="a";
        elseif fil==2
            filname="b";
        end
        stats=struct2cell(statobj.(type));
        if class(stats{1})=="Filament"
            for i=1:size(stats,1)
                stat_=stats{i};
                stats{i}=stat_.(filname);
            end

            if size(stats{2},2)==1 && size(stats{3},2)==1
                stats=cell2mat(stats);
                [stats,sortIn]=sort(stats);
                fit = polyfit(tempNList.(type)(sortIn(:)),stats(:),1);
                fxn=@(n) fit(1)*n+fit(2);
            else
                mat=obj.stattable(stat,type);
                mat=mat.(filname);
                fxn = scatteredInterpolant(mat(:,1),mat(:,2),mat(:,3),'linear','nearest');
            end
        elseif size(stats{1},1)==1
            if size(stats{1},2)==1 && size(stats{2},2)==1
                stats=cell2mat(stats);
                [stats,sortIn]=sort(stats);
                fit = polyfit(tempNList.(type)(sortIn),stats(:)',1);
                fxn=@(n) fit(1)*n+fit(2);
            else
                mat=obj.stattable(stat,type);
                fxn = scatteredInterpolant(mat(:,1),mat(:,2),mat(:,3),'linear','nearest');
            end
        end
    end
end