function generateFigs()
% GENERATEFIGS creates heatmaps of kpoly ratios for fictitious formin
% constructs swept across PRM location (NT and CT dist) as well as heatmaps
% of kpoly ratios for a fictitious formin construct swept across input
% probability density and occlusion probability ratios
    %
    %   Generates and saves two figures:
    %   1) heatmaps of polymerization rate ratios (dimer/nondimer)
    %      for sweeps across possible probability density and occlusion probability
    %      ratios for 3 parameter regimes:
    %           1) k_cap=1; k_del=10^3
    %           2) k_cap=10^4; k_del=1
    %           3) k_cap=10^7; k_del=1
    %      See MAKEKINETICHEATMAPS for more information.
    %      The 1x3 heatmap figure is saved to runs/kineticHeatmaps.png
    %   2) heatmap of kpoly ratios for a fictitious formin construct swept 
    %      across PRM location (NT and CT dist) for 3 parameter regimes:
    %           1) k_cap=1; k_del=10^5
    %           2) k_cap=1.5x10^4; k_del=1
    %           3) k_cap=10^9; k_del=1
    %      See MAKECOMBOSWEEPS for more information.
    %      The 1x3 heatmap figure is saved to runs/kpolyHeatmaps.png
    %   
    %   Loads lookup_35_stripped from lookup_35_stripped.mat
    %   Runs makeKineticHeatmaps and makecombosweeps
    %
    %   Must be run from the formin_kinetic_model/analysis folder
    % 
    %   See also MAKEKINETICHEATMAPS, POCCPRSWEEP, KPOLY_IND, MAKECOMBOSWEEPS.

    currentFolder = pwd;
    targetFolder = 'formin_kinetic_model/analysis'; 
    isMatch = endsWith(currentFolder, targetFolder);
    if ~isMatch
        error(['Must run from ''', targetFolder, ''''])
    end
    cd("../..")
    p=genpath("formin_kinetic_model");
    addpath(p)
    cd(currentFolder)

    fullfig = makeKineticHeatmaps();
    exportgraphics(fullfig,"../runs/kineticHeatmaps.png")

    load("lookup_35_stripped.mat","lookup_35_stripped")
    lt=Lookuptable(lookup_35_stripped);

    fullfig = makecombosweeps(lt);
    exportgraphics(fullfig,"../runs/kpolyHeatmaps.png")

end