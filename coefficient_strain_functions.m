
%Equilibrium (∞) Coefficient:
function [ a_ui ] = alpha_u(i)
a_ui = 0.058*log(i)+ 0.0034;
end

%Fast (f) Coefficient:
function [ a_fi ] = alpha_f(i)
a_fi = 0.71 - (0.016 * i);
end

%Modified Fast (f) Coefficient:
function [ a_fi ] = alpha_f_mod(i, S)
a_fi = 0.71 - ((0.71/S) * i);
end

%Slow (s) Coefficient:
function [ a_si ] = alpha_s(a_ui, a_fi)
a_si = 1 - a_ui - a_fi;
end

%Fast (f) Time Constant:
Constant (0.97)

%Slow (s) Time Constant:
function [ t_si ] = tau_s(i)
t_si = 7.42*i + 60.2;
end

%Maximum E-Modulus:
function [ E_0i ] = E_modulus_0(i)
E_0i = 81.3*i - 78.1;
end


%Calculation of Longitudinal Distraction Osteogenesis (LDO) Strain:
function [ e_neu ] = strain_L(l_delta, l_0, i, e_i, a_ui)
e_neu = (l_delta/(l_0 + l_delta*(i-1))) + a_ui * e_i;
% Calculation of the current stress in the callus tissue 
% considering the residual stress in the tissue
end

Calculation of Radial Distraction Osteogenesis (RDO) Strain:
function [ e_neu ] = strain_C(l_delta, l_0, i, e_i, a_ui)
e_neu = (l_delta/(l_0 + l_delta*(i-1))) + a_ui * e_i;
end

%Calculation of LDO/RDO E-Modulus:
function [ E_ti ] = E_modulus(E_0i, a_ui, a_fi, a_si, t_fi, t_si, t, k)
E_ti = (E_0i/k) * (a_ui + ((a_fi * exp(-t/t_fi)) + (a_si * exp(-t/t_si))));
% Calculation of the time-dependent E-Modulus at the time points 
% of distraction
end