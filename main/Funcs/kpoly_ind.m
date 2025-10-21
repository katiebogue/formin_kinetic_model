function kp = kpoly_ind(p_occ,p_r,pocc_base,c_PA,G,k_cap,k_del,r_cap,r_cap_exp,prm_size)
%KPOLY_IND Computes the rate of polymerization for a single binding site (PRM)
    % kp=
    %   KPOLY_IND(p_occ,p_r,pocc_base,c_PA,k_cap,G,k_del,r_cap,r_cap_exp,prm_size) computes the rate of
    %   elongation
    %
    %   Inputs:
    %         p_occ : (double) occlusion probability of the PRM
    %         p_r : (double) probability density of the PRM at the delivery
    %         site
    %         pocc_base : (double) occlusion probability of the delivery
    %         site
    %         c_PA : (double) concentration of profilin-actin
    %         G : (double) gating factor
    %         k_cap : (double) capture rate constant
    %         k_del : (double) delivery rate constant
    %         r_cap : (double) reverse capture rate constant
    %         r_cap_exp : (double) reverse capture exponential constant
    %         prm_size : (double) size (number of amino acids) of the PRM
    %
    %   Output is the computed polymerization rate (in array form if inputs
    %   are arrays)
    %
    %   Note-- all calculations are done with matrix multiplication or
    %   division

% rate of capture
kcap = (1-p_occ).*c_PA.*k_cap;

% rate of delivery
kdel = (1-pocc_base).*p_r.*k_del.*G;

% rate of reverse capture
rcap = r_cap.*exp(-prm_size.*r_cap_exp);

% overall polymerization rate
kp =  1./((1./kdel) + ((kdel + rcap)./(kdel.*kcap)));

end