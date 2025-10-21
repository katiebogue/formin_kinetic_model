function generateFigs()

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