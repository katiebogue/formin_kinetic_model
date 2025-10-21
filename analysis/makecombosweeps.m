function fullfig = makecombosweeps(lookuptab)

    fullfig = figure('units','centimeters','position',[0,5,55,18]);hold on;
    tiles = tiledlayout(1,3,'TileSpacing','tight','Padding','none');

    pythonpath=""; % path to python files -- not functional in this repo version
    resultsloc=""; % not saving anything
    
    opts_1=Options(lookuptab,pythonpath,...
        "3st",...       % kpoly type
        resultsloc,...
        10^3.2778,...          % k_cap
        10^-2.8678,... % k_del
        10^3.15,...    % r_cap
        1,...    % r_del
        1);     % k_rel
    
    opts_1.r_cap_exp=0.57121;
    opts_1.set_equation(1); % using preset #1 (see Options class)
    opts_1.NTopt=2; 
    
    opts_1.resultsfolder=strcat(opts_1.resultsfolder,"combosweep_1");
    
    opts_1.k_cap=1;
    opts_1.k_del=1;
    opts_1.r_cap=10^4;
    opts_1.r_cap_exp=0.8;
    
    opts_1.k_cap=1;
    opts_1.k_del=10^5;
    
    [fig,h]=kpolyheatmap(opts_1,[-3 3]);

    cmc1 = fig.Children;
    h1 = copyobj(cmc1,tiles);
    h1.Layout.Tile = 1;
    h1.NodeChildren(3).YDir='normal';
    close(fig)
    
    opts_1.k_cap=15000;
    opts_1.k_del=1;

    [fig,h]=kpolyheatmap(opts_1,[-3 3]);

    cmc2 = fig.Children;
    h2 = copyobj(cmc2,tiles);
    h2.Layout.Tile = 2;
    h2.NodeChildren(3).YDir='normal';
    close(fig)

    opts_1.k_cap=10^9;
    opts_1.k_del=1;

    [fig,h]=kpolyheatmap(opts_1,[-3 3]);

    cmc3 = fig.Children;
    h3 = copyobj(cmc3,tiles);
    h3.Layout.Tile = 3;
    h3.NodeChildren(3).YDir='normal';
    close(fig)
end

