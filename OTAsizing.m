clc; 
close all; 
clear all; 
format short;

%% Library extraction
%-----------------------------------------------------
% you can access to : 
% Vgs, Vdsat, Vea [V]
% In [A]
% gm/Id [1/V]
% Cgs [F/m2]
% Cgso, Cgdo, Cbd [F/m]
%-----------------------------------------------------
% To access the gm/Id of an LVT nMOS for example : nlvtlp.GMID
%-----------------------------------------------------
rowHeadings = {'VGS','GMID','IN','CGS','CGSO','CGDO','CBD','VEA','VDSAT'};
nlvtlp = cell2struct(textscan(fopen('./data/nlvtlp.txt'), '%f %f %f %f %f %f %f %f %f ', 'headerlines', 2, 'Delimiter', 'whitespace'),rowHeadings,2);
nsvtlp = cell2struct(textscan(fopen('./data/nsvtlp.txt'), '%f %f %f %f %f %f %f %f %f ', 'headerlines', 2, 'Delimiter', 'whitespace'),rowHeadings,2);
nhvtlp = cell2struct(textscan(fopen('./data/nhvtlp.txt'), '%f %f %f %f %f %f %f %f %f ', 'headerlines', 2, 'Delimiter', 'whitespace'),rowHeadings,2);
plvtlp = cell2struct(textscan(fopen('./data/plvtlp.txt'), '%f %f %f %f %f %f %f %f %f ', 'headerlines', 2, 'Delimiter', 'whitespace'),rowHeadings,2);
psvtlp = cell2struct(textscan(fopen('./data/psvtlp.txt'), '%f %f %f %f %f %f %f %f %f ', 'headerlines', 2, 'Delimiter', 'whitespace'),rowHeadings,2);
phvtlp = cell2struct(textscan(fopen('./data/phvtlp.txt'), '%f %f %f %f %f %f %f %f %f ', 'headerlines', 2, 'Delimiter', 'whitespace'),rowHeadings,2);
fclose('all');

%% Data modification: suppress first points to obtain a monotonic function (for interpolations)
M = max(nlvtlp.GMID);   Inl = find(nlvtlp.GMID==M, 1, 'last');
M = max(nsvtlp.GMID);   Ins = find(nsvtlp.GMID==M, 1, 'last');
M = max(nhvtlp.GMID);   Inh = find(nhvtlp.GMID==M, 1, 'last');
M = max(plvtlp.GMID);   Ipl = find(plvtlp.GMID==M, 1, 'last');
M = max(psvtlp.GMID);   Ips = find(psvtlp.GMID==M, 1, 'last');
M = max(phvtlp.GMID);   Iph = find(phvtlp.GMID==M, 1, 'last');

nlvtlp.VGS = nlvtlp.VGS(Inl:end);   nlvtlp.GMID = nlvtlp.GMID(Inl:end);   nlvtlp.IN = nlvtlp.IN(Inl:end);   nlvtlp.CGS = nlvtlp.CGS(Inl:end);   nlvtlp.CGSO = nlvtlp.CGSO(Inl:end);   nlvtlp.CGDO = nlvtlp.CGDO(Inl:end);   nlvtlp.CBD = nlvtlp.CBD(Inl:end);   nlvtlp.VEA = nlvtlp.VEA(Inl:end);     nlvtlp.VDSAT = nlvtlp.VDSAT(Inl:end);
nsvtlp.VGS = nsvtlp.VGS(Ins:end);   nsvtlp.GMID = nsvtlp.GMID(Ins:end);   nsvtlp.IN = nsvtlp.IN(Ins:end);   nsvtlp.CGS = nsvtlp.CGS(Ins:end);   nsvtlp.CGSO = nsvtlp.CGSO(Ins:end);   nsvtlp.CGDO = nsvtlp.CGDO(Ins:end);   nsvtlp.CBD = nsvtlp.CBD(Ins:end);   nsvtlp.VEA = nsvtlp.VEA(Ins:end);     nsvtlp.VDSAT = nsvtlp.VDSAT(Ins:end);
nhvtlp.VGS = nhvtlp.VGS(Inh:end);   nhvtlp.GMID = nhvtlp.GMID(Inh:end);   nhvtlp.IN = nhvtlp.IN(Inh:end);   nhvtlp.CGS = nhvtlp.CGS(Inh:end);   nhvtlp.CGSO = nhvtlp.CGSO(Inh:end);   nhvtlp.CGDO = nhvtlp.CGDO(Inh:end);   nhvtlp.CBD = nhvtlp.CBD(Inh:end);   nhvtlp.VEA = nhvtlp.VEA(Inh:end);     nhvtlp.VDSAT = nhvtlp.VDSAT(Inh:end);
plvtlp.VGS = plvtlp.VGS(Ipl:end);   plvtlp.GMID = plvtlp.GMID(Ipl:end);   plvtlp.IN = plvtlp.IN(Ipl:end);   plvtlp.CGS = plvtlp.CGS(Ipl:end);   plvtlp.CGSO = plvtlp.CGSO(Ipl:end);   plvtlp.CGDO = plvtlp.CGDO(Ipl:end);   plvtlp.CBD = plvtlp.CBD(Ipl:end);   plvtlp.VEA = plvtlp.VEA(Ipl:end);     plvtlp.VDSAT = plvtlp.VDSAT(Ipl:end);
psvtlp.VGS = psvtlp.VGS(Ips:end);   psvtlp.GMID = psvtlp.GMID(Ips:end);   psvtlp.IN = psvtlp.IN(Ips:end);   psvtlp.CGS = psvtlp.CGS(Ips:end);   psvtlp.CGSO = psvtlp.CGSO(Ips:end);   psvtlp.CGDO = psvtlp.CGDO(Ips:end);   psvtlp.CBD = psvtlp.CBD(Ips:end);   psvtlp.VEA = psvtlp.VEA(Ips:end);     psvtlp.VDSAT = psvtlp.VDSAT(Ips:end);
phvtlp.VGS = phvtlp.VGS(Iph:end);   phvtlp.GMID = phvtlp.GMID(Iph:end);   phvtlp.IN = phvtlp.IN(Iph:end);   phvtlp.CGS = phvtlp.CGS(Iph:end);   phvtlp.CGSO = phvtlp.CGSO(Iph:end);   phvtlp.CGDO = phvtlp.CGDO(Iph:end);   phvtlp.CBD = phvtlp.CBD(Iph:end);   phvtlp.VEA = phvtlp.VEA(Iph:end);     phvtlp.VDSAT = phvtlp.VDSAT(Iph:end);

%% Plot technology curves
% obtain_techno_curves({nlvtlp nsvtlp nhvtlp}, {'nlvtlp', 'nsvtlp', 'nhvtlp'});
% obtain_techno_curves({plvtlp psvtlp phvtlp}, {'plvtlp', 'psvtlp', 'phvtlp'});

%% Specifications
Vdd = 1.2;
Vss = -1.2;
CL = 5e-12; %[F]
fT = 5e6;
omega_u = 2*pi*fT;
Pm = 80;
SR = 5e6; %[V/s]
Cc = 2.5e-12;
L = 2e-6;
Vout_max = 2;
u_p = 0.01156; %[m^2/V*s]
u_n = 0.0506; %[m^2/V*s]
V_CM = 1;
V_CMHR = Vdd-V_CM;

L1 = L; L2 = L1;
L3 = L; L4 = L3; 
L5 = L; L6 = L5;
L7 = L; L8 = L7;
L9 = L; L10= L9;
L9A = L; L9B = L;
L9C = L; L9D = L;

%% Design choices
%gmid6 = 10; must be determine from the value of the tension
%gmid7 = gmid6;


M1 = nlvtlp;
M2 = nlvtlp;
M3 = plvtlp;
M4 = plvtlp;
M5 = nlvtlp;
M6 = plvtlp;
M7 = nlvtlp;
M8 = nlvtlp;
M9 = nhvtlp;
M10 = nhvtlp;
M9A = nhvtlp;
M9B = nhvtlp;
M9C = phvtlp;
M9D = phvtlp;

%% Design algorithm without gm/Id
% L6 = sqrt((3*u_p*V_HR*Cc)/(2*fT*2*pi*(Cc+CL)*tand(Pm)));
% W6 = (2*SR*(CL+Cc)*L6)/(u_p*Cox*(V_HR)^2);
% WL6 = W6/L6
%Vsg6 = V_HR+abs(V_tp)
%test6 = interp1(M6.VGS, M6.GMID, V_HR+abs(V_tp),'spline')
% 

% WL1 = (omega_u^2*Cc)/(u_n*Cox*SR)
% WL5 = (2*SR*Cc)/(u_n*Cox*(V_CMHR-V_tn-SR/omega_u)^2)
% WL7 = ((Cc+CL)/Cc)*WL5
% WL3 = ((WL6)/2*WL7)*WL5
% 
%WL9 =((2*Cc*SR)/(u_p*Cox*V_HR*(Vdd-V_HR-2*abs(V_tp))))



%% Design algorithm with gm/Id
% in1 = interp1(M1.GMID, M1.IN, gmid1,'spline');
% %peut d�finir le gain! Veaeq?
% %Veaeq = 1/(1/interp1(M6.GMID, M6.VEA, gmid6,'spline') + 1/interp1(M7.GMID, M7.VEA, gmid7,'spline')); 
% %Gain = gmid1*Veaeq*gmid6*Veaeq;
% gmid1 = 16.5;
gmid1 = 16.5;
in1 = interp1(M1.GMID, M1.IN, gmid1,'spline');
gm1 = 2*pi*fT*Cc;
id1 = gm1/gmid1;
W1 = id1/in1*L1;
WL1 = W1/L1
W2 = W1;
gmid2 = gmid1;
id2 = id1;
% gmid6 = 30;
% Cox = interp1(M6.GMID, M6.CGS, gmid6,'spline');
% % W2 = W1;
% % 
% L6 = sqrt((3*u_p*V_HR*Cc)/(2*fT*2*pi*(Cc+CL)*tand(Pm)));
% W6 = (2*SR*(Cc+CL)*L6)/(u_p*Cox*(V_HR)^2);
% WL_6 = W6/L6
% 
gmid5 = 9;
id5 = Cc*SR;

in5 = interp1(M5.GMID, M5.IN, gmid5,'spline');
W5 = id5/in5*L5;
WL5 = W5/L5
W8 = W5;
WL8 = WL5;
in8 = in5;
id8 = in8*WL8;
Vgs8 = interp1(M8.GMID, M8.VGS, gmid5, 'spline');
Vg8 = Vgs8+Vss;
R_bias = (Vdd-Vg8)/id8;

WL7 = ((Cc+CL)/(Cc))*WL5;
W7 = WL7*L7;

%gmid7 = 5;
gmid7 = 10;
in7 = interp1(M7.GMID, M7.IN, gmid7,'spline');
id7 = in7*WL7;
id6 = id7;

gmid6 = 12.5;
gm6 = gmid6*id6;
omega_T6 = tand(Pm)*omega_u*((Cc+CL)/Cc); %formula 30
W6L6 = gm6/(omega_T6*2/3*interp1(M6.GMID, M6.CGS, gmid6,'spline')); %formula 31
in6 = interp1(M6.GMID, M6.IN, gmid6,'spline');
W6sL6 = id6/in6;
L6 = sqrt(W6L6/W6sL6);
W6 = W6L6/L6;
WL6 = W6/L6


WL3 = WL6/(2*WL7) * WL5; 
W3 = WL3*L3;
WL4 = WL3;

 
Vsg6 = interp1(M6.GMID, M6.VGS, gmid6,'spline');
V_tp = 0.243;
V_HR = Vsg6-abs(V_tp);
tox = 2.3926e-9;
E0 = 8.854e-12; %[F/m]
Er = 3.9; %for SiO2;
Cox = E0*Er/tox;



W9 = L9*((2*Cc*SR)/(u_p*Cox*V_HR*(Vdd-V_HR-2*abs(V_tp))));
WL9 = W9/L9;

gd2 = id2/interp1(M2.GMID, M2.VEA, gmid2,'spline');
id4 = id1;
in4 = id4/WL4;
gmid4 = interp1(M4.IN, M4.GMID, in4, 'spline');
gd4 = id4/interp1(M4.GMID, M4.VEA, gmid4,'spline');
gd6 = id6/interp1(M6.GMID, M6.VEA, gmid6,'spline');
gd7 = id7/interp1(M7.GMID, M7.VEA, gmid7,'spline');
gmid7test = [0:1:15];
in7test = interp1(M7.GMID, M7.IN, gmid7test,'spline');
id7test = in7test*WL7;
gd7test = id7test./interp1(M7.GMID, M7.VEA, gmid7test,'spline');
%semilogy(gmid7test,gd7test);
RA = 1/(gd2+gd4);
RB = 1/(gd6+gd7);
A0 = gm1*gm6*RA*RB;
A0_dB = 20*log10(A0);

Pole_dp = 1/(gm6*RA*RB*Cc);
numerateur = gm6*Cc;
denominateur = (CL*2/3*W6*L6*interp1(M6.GMID, M6.CGS, gmid6,'spline'));
Pole_dn = gm6*Cc/(CL*2/3*W6*L6*interp1(M6.GMID, M6.CGS, gmid6,'spline'));

num = A0;
den=conv([1/Pole_dn 1],[1/Pole_dp 1]);
sys = tf(num,den);
sys2 = feedback(sys,1);
vgs1 = interp1(M1.GMID, M1.VGS, gmid1,'spline');
Vdsat5 = interp1(M5.GMID, M5.VDSAT, gmid5, 'slpine');
Vin = vgs1+Vdsat5-Vdd;
%!!!!!!!!-Vdd!!!!!!!!!!!

figure;
bode(sys,sys2,{1e0,1e12})
grid on
axes_handles = findall(gcf, 'type', 'axes');
legend(axes_handles(3),'Open-loop gain','Closed-loop gain')


figure;
margin(sys)
xlim([1 1e12])
grid on

[Gm,Pm,Wg,Wp] = margin(sys);
Pm

    


%---------------------------------------------
% one = W1/L1
% two = W2/L2
% three = W3/L3
% four = W4/L4
% five = W5/L5
% six = W6/L6
% seven = W7/L7
% eight = W8/L8
% nine = W9/L9
%---------------------------------------------
% Gain, poles and zeros
%---------------------------------------------
%---------------------------------------------

%% Plots

% num=Gain*[1/Zeron 1];
% den=conv([1/Polep 1],conv([1/Polen 1],[1/Poled 1]));
% 
% sys = tf(num,den);
% sys2 = feedback(sys,1);
% 
% figure;
% bode(sys,sys2,{1e0,1e12})
% grid on
% axes_handles = findall(gcf, 'type', 'axes');
% legend(axes_handles(3),'Open-loop gain','Closed-loop gain')
% 
% figure;
% margin(sys)
% xlim([1 1e12])
% grid on
% 
% [Gm,Pm,Wg,Wp] = margin(sys);
% 
% %% Algorithm results
% 
 fid = fopen('OTA.txt','wt');
 fprintf(fid,'------------------------------------------------------------------\n');
 fprintf(fid,'Performances\n');
 fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'Gain: %.4g [dB]\n', 20*log10(Gain));
 fprintf(fid,'Transition frequency: %.4g [MHz]\n',Wp/(2*pi)*1e-6);
% fprintf(fid,'Dominant pole: %.4g [kHz]\n',Poled/(2*pi)*1e-3);
% fprintf(fid,'Position of polep: %.5g [MHz]\n',Polep/(2*pi)*1e-6);
% fprintf(fid,'Position of polen: %.5g [MHz]\n',Polen/(2*pi)*1e-6);
% fprintf(fid, 'Phase margin: %.4g [degrees]\n',Pm);
% fprintf(fid,'Bias current: %.4g [uA]\n',ibias*1e6); 
% fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'DC point\n');
% fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'Vgs1: %3.2g [V]\n', vgs1);
% fprintf(fid,'Vsg4: %3.2g [V]\n', vsg4);
% fprintf(fid,'Vsg6: %3.2g [V]\n', vsg6);
% fprintf(fid,'Vgs7: %3.2g [V]\n', vgs7);
% fprintf(fid,'Vgs9: %3.2g [V]\n', vgs9);
% fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'Dimensions\n');
% fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'W1: %3.2g [um]\n', W1*1e6);
% fprintf(fid,'W4: %3.2g [um]\n', W4*1e6);
% fprintf(fid,'W6: %3.2g [um]\n', W6*1e6);
% fprintf(fid,'W7: %3.2g [um]\n', W7*1e6);
% fprintf(fid,'W9: %3.2g [um]\n', W9*1e6);
% fclose(fid);