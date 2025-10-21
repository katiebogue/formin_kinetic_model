function [fig,h]=kpolyheatmap(options,col,smooth,smoothfac,save)
% KPOLYHEATMAP  creates heatmap of kpoly ratios for a fictitious formin
% construct swept across PRM location (NT and CT dist)
    %
    %   [fig,h]=KPOLYHEATMAP(options,col) creates heatmap of NT and CT dist sweep
    %   with the specified colormap options, uses
    %   smoothing factor of 0.58
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,0) creates heatmap of NT and CT dist sweep
    %   with the specified colormap options, with
    %   no smoothing
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,1,smoothfac) creates
    %   heatmap of NT and CT dist sweep with the specified colormap options, 
    %   uses specified smoothing factor 
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,0,x,1) creates and saves
    %   heatmap of NT and CT dist sweep with the specified
    %   colormap options, with no smoothing
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,1,smoothfac,1) creates and
    %   saves heatmap of NT and CT dist sweep with the specified
    %   colormap options, uses specified smoothing factor
    %
    %   Inputs:
    %       options    : Options object with parameters and rate equations
    %       col        : (double) colormap type
    %                       1- use jet
    %                       2- use parula, if smooth not true, set all increasing value to red
    %                       3- load 'customcolorbar_red_blue.mat' and set color limits to min values
    %                       [x, x] - use 'customcolorbar_red_blue.mat' and
    %                               set colorlimits to [x, x]
    %       smooth     : (logical) whether or not to use smoothing on the
    %                     heatmap (default is true)
    %       smoothfac  : (double) smoothing factor, applied if smooth=true
    %                    (default is 0.58)
    %       save       : (logical) whether or not to save the resulting
    %                    heatmap (default is false)
    %
    %   Outputs:
    %       fig : Figure with heatmap
    %       h   : heatmap object
    %
    %
    %   Loads customcolorbar_red_blue.mat
    % 
    %   See also FORMIN, OPTIONS, PRM, KPOLYMERIZATION.

arguments
    options Options
    col double
    smooth logical=true
    smoothfac double =0.58
    save logical=false
end
fig=figure;

if save
    set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs
end


kpoly_lab="log_{2}(k_{poly} dimerized/k_{poly} non-dimerized)";

x_label="Distance from PRM to FH2";
y_label="Distance from PRM to N-term";
scatn=0;
PRM_size=10;

PRM_lab=strcat("PRM size: ",num2str(PRM_size));

yvals=1:300;
xvals=1:300;

LMAX=200;

kpolyratios=zeros(length(yvals)*length(xvals),1);
yvals_matrix=kpolyratios;
xvals_matrix=kpolyratios;

index=0;

for i=1:length(xvals)
    PRM_loc=xvals(i);

    formin1=Formin("formin1",options,c_PA=1,gating=1,length=PRM_loc,PRMloc=PRM_loc,PRMsize=PRM_size);
    for n=1:length(yvals)
        FH1_length=PRM_loc+yvals(n);
        if FH1_length>LMAX
            index=index+1;
            kpolyratios(index,1)=nan;
            yvals_matrix(index,1)=yvals(n);
            xvals_matrix(index,1)=xvals(i);
        else
            addlen=FH1_length-formin1.length;
            formin1.add_length(addlen);
            index=index+1;
            kpolyratios(index,1)=formin1.kpoly.ratio;
            yvals_matrix(index,1)=yvals(n);
            xvals_matrix(index,1)=xvals(i);
        end
    end
end

if smooth
    for n=1:length(yvals)
        index=yvals_matrix==yvals(n);
        vals=kpolyratios(index);
        smootheddata=smoothdata(vals,"lowess","SmoothingFactor",smoothfac);
        smootheddata(smootheddata<0)=vals(smootheddata<0);
        kpolyratios(index)=smootheddata;
    end
    kpoly_lab=strcat(kpoly_lab," (smoothed:",num2str(smoothfac),")");
end

tbl=table(xvals_matrix,yvals_matrix,log2(kpolyratios));

h = heatmap(tbl,'xvals_matrix','yvals_matrix','ColorVariable','Var3');
h.XLabel = x_label;
h.YLabel = y_label;
h.ColorMethod = 'none';
h.GridVisible="off";
h.NodeChildren(3).YDir='normal'; 
if col==1
    h.Colormap=[jet];
elseif col==2
    h.Colormap=[parula];
    
    if not(smooth)
        red = [1 0 0];
        h.Colormap=[parula;red];
        h.ColorLimits = [min(log2(kpolyratios)) 0];
    end
elseif col==3
    load('customcolorbar_red_blue.mat');
    h.Colormap=CustomColormap;
    minn=min(log2(kpolyratios(kpolyratios~=0)));
    h.ColorLimits=[minn abs(minn)];
elseif length(col)==2
    load('customcolorbar_red_blue.mat');
    h.Colormap=CustomColormap;
    h.ColorLimits=col;
end


prvec=getmeanstat(options.lookup,'Prvec0',PRM_loc);
pocc=getmeanstat(options.lookup,'POcclude',PRM_loc);
pocc0=getmeanstat(options.lookup,'POcclude',1);

kcapavg=log10(options.k_cap*(1-pocc));
kdelavg=log10(options.k_del*(1.0e33*(prvec)/(27*6.022e23))*(1-pocc0));

lab_avg=strcat("kcap avg: ",num2str(kcapavg)," kdel avg: ",num2str(kdelavg));

lab=strcat("k_{cap}: ",num2str(options.k_cap)," k_{del}: ",num2str(options.k_del)," r_{cap}: ",num2str(options.r_cap));

h.Title = {kpoly_lab,PRM_lab,lab,lab_avg};

% Convert each number in the array into a string
CustomXLabels = string(xvals);
CustomYLabels = string(yvals);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(xvals,20) ~= 0) = " ";
CustomYLabels(mod(yvals,20) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;
s = struct(h); 
s.XAxis.TickLabelRotation = 0;   % horizontal

if save
    fname=append("heatmapsweep_bounded_",sweep_type,"kcap ", num2str(options.k_cap)," kdel ", num2str(options.k_del)," rcap ", num2str(options.r_cap),'.png');
    saveas(gcf,fname);
    
    load('customcolorbar_red_blue.mat');
    h.Colormap=CustomColormap;
    minn=min(log2(kpolyratios(kpolyratios~=0)));
    h.ColorLimits=[minn abs(minn)];

    fname=append("heatmapsweep_",sweep_type,"kcap ", num2str(options.k_cap)," kdel ", num2str(options.k_del)," rcap ", num2str(options.r_cap),'.png');
    saveas(gcf,fname);
end
end

function value=getmeanstat(lookup,stat,isite)
    a=lookup.getstat(type='double',Stat=stat,iSite=isite,fil='a');
    a=structfun(@mean,a);
    a=mean(a);

    b=lookup.getstat(type='double',Stat=stat,iSite=isite,fil='b');
    b=structfun(@mean,b);
    b=mean(b);

    dob=mean([a,b]);

    a=lookup.getstat(type='dimer',Stat=stat,iSite=isite,fil='a');
    a=structfun(@mean,a);
    a=mean(a);

    b=lookup.getstat(type='dimer',Stat=stat,iSite=isite,fil='b');
    b=structfun(@mean,b);
    b=mean(b);

    dim=mean([a,b]);

    value=mean([dim,dob]);
end