function fullfig = makeKineticHeatmaps()
% MAKEKINETICHEATMAPS generates heatmaps of polymerization rate ratios (dimer/nondimer)
%for sweeps across possible probability density and occlusion probability
%ratios for 3 overall parameter regimes.
    % 
    % Makes heatmaps using the default parameters in POCCPRSWEEP plus the
    % following 3 capture and delivery rate constant combinations:
    %
    %   1) k_cap=1; k_del=10^3
    %   2) k_cap=10^4; k_del=1
    %   3) k_cap=10^7; k_del=1
    %
    % See also POCCPRSWEEP, KPOLY_IND.

    fullfig = figure('units','centimeters','position',[0,5,55,18]);hold on;
    tiles = tiledlayout(1,3,'TileSpacing','tight','Padding','none');
    
    k_cap=1;
    k_del=1000;
    nexttile(1)
    h1 = poccPrSweep(k_cap,k_del);
    
    k_cap=10000;
    k_del=1;
    nexttile(2)
    h2 = poccPrSweep(k_cap,k_del);
    
    k_cap=10000000;
    k_del=1;
    nexttile(3)
    h3 = poccPrSweep(k_cap,k_del);
    
end
