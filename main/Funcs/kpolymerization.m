function output = kpolymerization(type,kcap,kdel,rcap,rdel,krel)
%KPOLYMERIZATION Computes the rate of elongation for a single PRM
    % kpoly=
    %   KPOLYMERIZATION(type,kcap,kdel,rcap,rdel,krel) computes the rate of
    %   elongation
    %
    %   Inputs:
    %         type : (string) the elongation rate equation type; either
    %                 "capture", "3st", or "4st"
    %         kcap : (double or FilType) rate of capture for the PRM
    %         kdel : (double or FilType) rate of delivery for the PRM
    %         rcap : (double or FilType) reverse rate of capture for the 
    %                PRM
    %         rdel : (double or FilType) reverse rate of delivery for the 
    %                 PRM
    %         krel : (double or FilType) rate of release for the PRM
    %
    %   Output is either double or FilType (depending on inputs)
    %
    %   Note-- all calculations are done with matrix multiplication or
    %   division if the input is not of class "FilType"
    %
    %   See also KPOLY, PRM, FORMIN, FILTYPE.
arguments
    type (1,1) string {mustBeMember(type,{'capture','3st','4st'})}
    kcap {mustBeA(kcap,["double","FilType"])}
    kdel {mustBeA(kdel,["double","FilType"])}
    rcap {mustBeA(rcap,["double","FilType"])}
    rdel {mustBeA(rdel,["double","FilType"])}
    krel {mustBeA(krel,["double","FilType"])}
end
    if type=="capture"
        output=kcap;
    elseif type=="3st"
        if class(kcap)=="FilType" || class(kdel)=="FilType" || class(rcap)=="FilType"
            output=1/((1/kdel) + ((kdel + rcap)/(kdel*kcap)));
        else
            output=1./((1./kdel) + ((kdel + rcap)./(kdel.*kcap)));
        end
    elseif type=="4st"
        if class(kcap)=="FilType" || class(kdel)=="FilType" || class(rcap)=="FilType" || class(krel)=="FilType" || class(rdel)=="FilType" 
            output=1/((1/krel) + ((rdel + krel)/(kdel * krel)) + (((rcap * rdel) + (rcap * krel) + (kdel * krel))/(kcap * kdel * krel)));
        else
            output=1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel)));
        end
        %sum(1./((1./r_PF_rev)+((r_paf_rev+r_PF_rev)./((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))).*r_PF_rev))+((((k_paf_rev.*exp(-1.*pp_length_vec)).*r_paf_rev)+((k_paf_rev.*exp(-1.*pp_length_vec)).*r_PF_rev)+((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))).*r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2a(:))).*(k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))).*r_PF_rev))))+sum(1./((1./r_PF_rev)+((r_paf_rev+r_PF_rev)./((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))).*r_PF_rev))+((((k_paf_rev.*exp(-1.*pp_length_vec)).*r_paf_rev)+((k_paf_rev.*exp(-1.*pp_length_vec)).*r_PF_rev)+((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))).*r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2b(:))).*(k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))).*r_PF_rev))))
    end
end