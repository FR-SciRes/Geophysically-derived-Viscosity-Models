function [SR_Total, SR_Dis, SR_Dif, SR_DisGBS, SR_peierls] = TotalStrainRate_HK03_fH2O(stress, d, COH, P, T, melt) 

%The inputs are: stress in MPa, d in mm, COH in ppm H/Si, T in kelvin, P in
%GPa, melt in percent
%The stress and d can be an array. Other parameters are in value.
%The function outputs the ff: total strain rate, dislocation, diffusion,
%DisGBS, and peierls strain rates which could be in matrix form. The rows
%correspond to d elements while columns correspond to stress elements.


%The following are the flow law paramaters for different deformation
%mechanisms at dry and wet conditions, and the respective strain rate
%calculations.
P = P*10^9;         %in Pa
d = d*10^-3;        %in meter, note that in Hirth & Kohlstedt, the unit is in micrometer so I multiplied this new d with 10^6 (for diffusion and dislocation creep only)
melt = melt/100;    %in decimal form

%A. For DRY condition
%A.1 Diffusion Mechanism
A_dryDif = 1.5e9;       %MPa^(-n-r)micrometer^(p)s^-1; constant that carries oxygen fugacity and Si activity
n_dryDif = 1;           %stress exponent
p_dryDif = 3;           %grain size exponent
r_dryDif = 0;           %water fugacity exponent
a_dryDif = 30;          %melt constant
E_dryDif = 375*10^3;    %activation energy in J/mol (error of 50*10^3)
V_dryDif = 6*10^-6;  %activation volume in m^3/mol (range is 2-10*10-6)

%A.2 Dislocation Mechanism
A_dryDis = 1.1e5;    %constant that carries oxygen fugacity and Si activity
n_dryDis = 3.5;         %stress exponent (error of 0.3)
p_dryDis = 0;           %grain size exponent
r_dryDis = 0;           %water fugacity exponent
a_dryDis = 30;          %melt constant (30-45)
E_dryDis = 530*10^3;    %activation energy in J/mol (error of 4*10^3)
V_dryDis = 15*10^-6;    %activation volume in m^3/mol based on recovery experiment for P=0.0001-10 GPa

%A.3 DisGBS Mechanism
A_dryGBS = 10^-4.89;    %constant that carries oxygen fugacity and Si activity
n_dryGBS = 3.0;            %stress exponent (error of 0.3)
p_dryGBS = 1;               %grain size exponent
r_dryGBS = 1.25;               %water fugacity exponent
a_dryGBS = 30;             %melt constant (30-45)
E_dryGBS = 423*10^3;      %activation energy in J/mol (error of 75*10^3)
V_dryGBS = 17.6*10^-6;         %activation volume in m^3/mol , assumed to be the same in dislocation creep, from Kawazoe et al. (2009)

%B. For WET condition
%B.1 Diffusion Mechanism
A_wetDif = 2.5e7;    %constant that carries oxygen fugacity and Si activity 4.0*10^5
n_wetDif = 1;           %stress exponent
p_wetDif = 3;           %grain size exponent
r_wetDif = 1;%0.8;          %water fugacity exponent
a_wetDif = 30;          %melt constant
E_wetDif = 375*10^3;    %activation energy in J/mol (error of 75*10^3)335
V_wetDif = 14*10^-6;     %activation volume in m^3/mol [4*10^-6] from Ohuchi et al. (2012)
%using Zhao et al. 2004 fH2O-COH conversion
A_wetDif = 3.5^r_wetDif*A_wetDif;

%B.2 Dislocation Mechanism
A_wetDis = 1.6e3;        %constant that carries oxygen fugacity and Si activity
n_wetDis = 3.5;         %stress exponent (error of 0.3)
p_wetDis = 0;           %grain size exponent
r_wetDis = 1.2;         %water fugacity exponent
a_wetDis = 30;          %melt constant (30-45)
E_wetDis = 520*10^3;    %activation energy in J/mol (error of 40*10^3)
V_wetDis = 22*10^-6;    %activation volume in m^3/mol (error of 11*10^-6)--used by Ohuchi et al (2015)
%using Zhao et al. 2004 fH2O-COH conversion
A_wetDis = 3.5^r_wetDis*A_wetDis;

% %% From Karato and Jung (2003) open system
% A_wetDis = 10^2.9;        %constant that carries oxygen fugacity and Si activity
% n_wetDis = 3.0;         %stress exponent (error of 0.3)
% p_wetDis = 0;           %grain size exponent
% r_wetDis = 1.2;         %water fugacity exponent
% a_wetDis = 30;          %melt constant (30-45)
% E_wetDis = 470*10^3;    %activation energy in J/mol (error of 40*10^3)
% V_wetDis = 24*10^-6;    %activation volume in m^3/mol (error of 11*10^-6)--used by Ohuchi et al (2015)
% % using Zhao et al. 2004 fH2O-COH conversion
% A_wetDis = 3.5^r_wetDis*A_wetDis;


%B.3 DisGBS Mechanism (Dislocation accommodated grain boundary sliding
A_wetGBS = 10^-4.89;     %constant that carries oxygen fugacity and Si activity
n_wetGBS = 3.0;           %stress exponent (error of 0.3)
p_wetGBS = 1;             %grain size exponent
r_wetGBS = 1.25;          %water fugacity exponent
a_wetGBS = 30;            %melt constant (30-45)?
E_wetGBS = 423*10^3;      %activation energy in J/mol (error of 56*10^3)
V_wetGBS = 17.6*10^-6;    %activation volume in m^3/mol (error of 0.8*10^-6)

% A_wetGBS = 3.5^r_wetGBS*A_wetGBS;
%D. Peierls mechanism (Evan and Goeze, 1979)
A_p = 10^12.10;
n_p = 1;
E_p = 564*10^3;
V_p = 30*10^-6;           %Kawazoe et al. (2010)
stress_p = 9.1*10^3;        %in MPa

%E. Total strain rate calculation
for i=1:length(d);
for j=1:length(stress);
    
    fwater = fH2O(COH,P,T);
    SR_peierls(i,j) = A_p*exp(-((E_p+P*V_p)/(8.3144621*T))*(1-(stress(j)/stress_p))^2);
    
    if COH <= 100;
        SR_Dif(i,j) = StrainRate(A_dryDif, stress(j), n_dryDif, d(i)*10^6, p_dryDif, fwater, r_dryDif, a_dryDif, melt, E_dryDif, P, V_dryDif, T);
        SR_Dis(i,j) = StrainRate(A_dryDis, stress(j), n_dryDis, d(i)*10^6, p_dryDis, fwater, r_dryDis, a_dryDis, melt, E_dryDis, P, V_dryDis, T);
        SR_DisGBS(i,j) = StrainRate(A_dryGBS, stress(j), n_dryGBS, d(i), p_dryGBS, fwater, r_dryGBS, a_dryGBS, melt, E_dryGBS, P, V_dryGBS, T);
    else
        SR_Dif(i,j) = StrainRate(A_wetDif, stress(j), n_wetDif, d(i)*10^6, p_wetDif, fwater, r_wetDif, a_wetDif, melt, E_wetDif, P, V_wetDif, T);
        SR_Dis(i,j) = StrainRate(A_wetDis, stress(j), n_wetDis, d(i)*10^6, p_wetDis, fwater, r_wetDis, a_wetDis, melt, E_wetDis, P, V_wetDis, T);
        SR_DisGBS(i,j) = StrainRate(A_wetGBS, stress(j), n_wetGBS, d(i), p_wetGBS, fwater, r_wetGBS, a_wetGBS, melt, E_wetGBS, P, V_wetGBS, T);
    end
end
end
SR_Total = SR_Dif + SR_Dis + SR_DisGBS;% + SR_peierls;

%Used sub-functions:
function [f] = StrainRate(A, stress, n, d, p, fwater, r, alpha, melt, E, P, V, T)
%The function computes the strain rate based on power law rheology for certain deformation mechanism as
%described by f:
R = 8.3144621; %gas constant in J/mol.K
f = A*(stress^n)*(d^(-p))*(fwater^r)*exp(alpha*melt)*exp(-(E+P*V)/(R*T));
end

function [f] = fH2O(COH,P,T) %in MPa
%convert water in ppm H/Si to f(H2O) using Keppler and Bolfan-Cassanova(2006)
Ac= 90;%87.81;          %ppm H/Si/MPa (derived by using A(T)=1.1 ppm H/Si from Kohlstedt et al.(1997)
Vc=10*10^-6;      %m^3/mol (Zhao et al., 2004)
R = 8.3144621;      %gas constant in J/mol.K
Ec = 50000;         %J/mol (Zhao et al., 2004)
alpha = 97000;      %J/mol; +\- 4000
XFa = 0.10;         %Fo90 that means that fayalite (iron-bearing) is 10% 

f=COH/(Ac*exp(-(Ec+(P*Vc))/(R*T))*exp((alpha*XFa)/(R*T)));
end

% function [f] = fH2O(COH,P,T) %in MPa from HK03 Eq. 6
% %convert water in ppm H/Si to f(H2O) using Keppler and Bolfan-Cassanova(2006)
% Ac= 26;%87.81;          %ppm H/Si/MPa (derived by using A(T)=1.1 ppm H/Si from Kohlstedt et al.(1997)
% Vc=10*10^-6;      %m^3/mol (Zhao et al., 2004)
% R = 8.3144621;      %gas constant in J/mol.K
% Ec = 40000;         %J/mol (Zhao et al., 2004)
% 
% f=COH/(Ac*exp(-(Ec+(P*Vc))/(R*T)));
% end

end


