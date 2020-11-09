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
Vdd = 2.5;
Vss = -2.5;
CL = 5e-12; %[F]
fT = 5e6;
Pm = 65;
SR = 5e6; %[V/s]
Cc = 2.5e-12;
L = 2e-6;
Vout_max = 2;
u_p = 0.01156; %[m^2/V*s]
u_n = 0.0506; %[m^2/V*s]
V_HR = Vdd-Vout_max;
E0 = 8.854e-12; %[F/m]
Er = 3.9; %for SiO2;
tox = 9.6e-9; %[m]
Cox = E0*Er/tox;
V_tp = -0.901;

L1 = L; L2 = L1;
L3 = L; L4 = L3; 
L5 = L; L6 = L5;
L7 = L; L8 = L7;
L9 = L; L10= L9;
L9A = L; L9B = L;
L9C = L; L9D = L;

%% Design choices
gmid1 = 15;
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

%% Design algorithm
in1 = interp1(M1.GMID, M1.IN, gmid1,'spline');
%peut définir le gain! Veaeq?
%Veaeq = 1/(1/interp1(M6.GMID, M6.VEA, gmid6,'spline') + 1/interp1(M7.GMID, M7.VEA, gmid7,'spline')); 
%Gain = gmid1*Veaeq*gmid6*Veaeq;
gm1 = 2*pi*fT*Cc;
id1 = gm1/gmid1;
W1 = id1/in1*L1;
W2 = W1;

L6 = sqrt((3*u_p*V_HR*Cc)/(2*fT*2*pi*(Cc+CL)*tand(Pm)));
W6 = (2*SR*(Cc+CL)*L6)/(u_p*Cox*(V_HR)^2);
W_over_L_6 = W6/L6;

id5 = Cc*SR;

%%atention
in5 = interp1(M5.GMID, M5.IN, gmid1,'spline');
W5 = id5/in5*L5;
W8 = W5;

W7 = ((Cc+CL)/(Cc))*W5;

W3 = L3*( W_over_L_6/(2*W7) * W5 );
W4 = W3;

W9 = L9*((2*Cc*SR)/(u_p*Cox*V_HR*(Vdd-V_HR-2*abs(V_tp))));

one = W1/L1
two = W2/L2
three = W3/L3
four = W4/L4
five = W5/L5
six = W6/L6
seven = W7/L7
eight = W8/L8
nine = W9/L9

    


%---------------------------------------------

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
% fid = fopen('OTA.txt','wt');
% fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'Performances\n');
% fprintf(fid,'------------------------------------------------------------------\n');
% fprintf(fid,'Gain: %.4g [dB]\n', 20*log10(Gain));
% fprintf(fid,'Transition frequency: %.4g [MHz]\n',Wp/(2*pi)*1e-6);
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