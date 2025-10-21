function h = poccPrSweep(k_cap,k_del,r_cap,r_cap_exp,c_PA,G,prm_size,pocc_base_nondimer,pocc_base_dimer,p_occ_nondimer,p_r_nondimer)
%POCCPRSWEEP make heatmaps of polymerization rate ratios (dimer/nondimer)
%for sweeps across possible probability density and occlusion probability
%ratios
%
% Inputs:
%    k_cap : (double) capture rate constant
%    k_del : (double) delivery rate constant
%    r_cap : (double) reverse capture rate constant; default is 10000
%    r_cap_exp : (double) reverse capture exponent constant; default is 0.8
%    c_PA : (double) concentration of profilin-actin; default is 1
%    G : (double) gating factor; default is 1
%    prm_size : (double) size (number of amino acids) of the PRM; default is 10
%    pocc_base_nondimer : (double) occlusion probability of the delivery site for the non-dimerized formin; default is 0.75
%    pocc_base_dimer : (double) occlusion probability of the delivery site for the dimerized formin; default is 0.75
%    p_occ_nondimer : (double) occlusion probability of the PRM for the non-dimerized formin; default is 0.5
%    p_r_nondimer : (double) probability density of the PRM at the delivery site for the non-dimerized formin; default is 10^5
%
%   Generates heatmap of polymerization rate ratios for input ratios
%   ranging from:
%       2^-5 to 2^5 for probability density
%       50% to 100% for 1-occlusion probability
%
%   The polymerization rate for the non-dimerized formin is held constant,
%   while the polymerization rate for the dimerized formin varies according
%   to the probability density and occlusion probability ratios.
%
%   Note that accessibility probability (1-pocc) is plotted instead of occlusion
%   probability for interpretability
%       
% 
% See also KPOLY_IND.
arguments
    k_cap double
    k_del double
    r_cap double =10000
    r_cap_exp double =0.8 
    c_PA double =1
    G double =1
    prm_size double =10
    pocc_base_nondimer double =0.75
    pocc_base_dimer double =0.75
    p_occ_nondimer double =0.5
    p_r_nondimer double =10^5
end

% compute static polymerization rate for the non-dimerized formin
kp_nondimer = kpoly_ind(p_occ_nondimer,p_r_nondimer,pocc_base_nondimer,c_PA,G,k_cap,k_del,r_cap,r_cap_exp,prm_size);

% set range for pr and pocc ratio sweeps
p_r_ratio_sweep_vals=-5:0.05:5;
p_occ_comp_ratio_sweep_vals=[-1:0.005:-0.0000000001 0:0.005:log2(1/(1-p_occ_nondimer))];

% Generate all combinations of i and j using ndgrid
[I, J] = ndgrid(p_occ_comp_ratio_sweep_vals, p_r_ratio_sweep_vals);

% Flatten the grids
I_flat = I(:);
J_flat = J(:);

% Preallocate kpval array
Kpvals = arrayfun(@(x, y) kpolyratio(x, y), I_flat, J_flat);

% Combine into result matrix
vals = [I_flat, J_flat, Kpvals];
x=array2table(vals);

% generate heatmap
h=heatmap(x,"vals1","vals2",'ColorVariable',"vals3","ColorMethod","none","GridVisible","off");

% formating
CustomColormap = [];
load('customcolorbar_red_blue.mat','CustomColormap');
h.Colormap=CustomColormap;
h.ColorLimits=[-5 5];
h.NodeChildren(3).YDir='normal';
% Convert each number in the array into a string
CustomXLabels = string(p_occ_comp_ratio_sweep_vals);
CustomYLabels = string(p_r_ratio_sweep_vals);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(p_occ_comp_ratio_sweep_vals,0.5) ~= 0) = " ";
CustomYLabels(mod(p_r_ratio_sweep_vals,1) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;
s = struct(h); 
s.XAxis.TickLabelRotation = 0;   % horizontal

% Labeling
kpoly_lab="log_{2}(k_{poly} dimerized/k_{poly} non-dimerized)";
lab_rates=strcat("kcap: ",num2str(k_cap)," kdel: ",num2str(k_del)," rcap: ",num2str(r_cap)," rcap_{exp}: ",num2str(r_cap_exp));
lab_formin=strcat("c_{PA}: ",num2str(c_PA)," G: ",num2str(G)," PRM size: ",num2str(prm_size));
lab_sim1=strcat("p_{occ} non-dimerized: ",num2str(p_occ_nondimer)," p_{r} non-dimerized: ",num2str(p_r_nondimer));
lab_sim2=strcat(" p_{occ}^{0} non-dimerized: ",num2str(pocc_base_nondimer)," p_{occ}^{0} dimerized: ",num2str(pocc_base_dimer));

h.Title = {kpoly_lab,lab_rates,lab_formin,lab_sim1,lab_sim2};
h.YLabel = "log_{2}(p_{r} dimerized/p_{r} non-dimerized)";
h.XLabel = "log_{2}(1-p_{occ} dimerized/1-p_{occ} non-dimerized)";

    function kpratio = kpolyratio(p_occ_comp_ratio, p_r_ratio)
        p_r_dimer=p_r_nondimer.*(2.^p_r_ratio);

        p_occ_comp_ratio_nonlog=2.^p_occ_comp_ratio;
        p_occ_dimer=1-((1-p_occ_nondimer).*(p_occ_comp_ratio_nonlog));

        kp_dimer = kpoly_ind(p_occ_dimer,p_r_dimer,pocc_base_dimer,c_PA,G,k_cap,k_del,r_cap,r_cap_exp,prm_size);

        kpratio = log2(kp_dimer./kp_nondimer);
    end
end