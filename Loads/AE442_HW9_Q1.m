v_e = linspace(0,1000);
W = 14000; %lbs
S = 200; %ft^2
rho = 0.0008; %slugs/ft^2
h = 30000;
clmax = 2;
clmin = -.8;
cla = 4.7;
c = 8;
g = 32.17;
Tmax = 4250;
Re = 1496111;

CL_max = 0.9*clmax;
CN_max = 1.1*CL_max;

CL_maxneg = 0.9*clmin;
CN_maxneg = 1.1*CL_maxneg;

n_pos = 7.33; %top horizontal
n_neg = -3; %bottom horizontal

syms nz
V_s = sqrt((2*(W/S))/(rho*CN_max))*0.592484;
V_sn = V_s*sqrt(nz); %left top sloped line. nZ chosen to be 1 from Roskam

V_sneg = sqrt((2*-1*(W/S))/(rho*CN_maxneg))*0.592484;
V_snneg = V_sneg*sqrt(nz);

kc = 36;
V_C = 425;
V_D = V_C*1.25;
V_H = V_C*1.4;%in this case max cruise speed
V_L = 1.25*V_H;

n_pos_range = [linspace(0,1,12) linspace(1.1016,n_pos,23)]; %Vstall at index 12.
n_neg_range = [linspace(0,1,12) linspace(1.1016,-n_neg,23)];
pos_stall = double(subs(V_sn, nz, n_pos_range));
neg_stall = double(subs(V_snneg, nz, n_neg_range));


V_stall = pos_stall(12);
V_dive  = V_D;
V_cruise = V_C;
V_maneuv = pos_stall(end);

%%
rhos = 0.0008; %slugs/ft^2
miug = (2*W/S)/(rhos*g*c*cla);
Kg = miug^1.03/(6.9+miug^1.03);
syms V_B
UdeB =39;
UdeC = 29; %for VC
UdeD = 15;
nlimgpC = 1+Kg*UdeC*V_C*cla/(498*W/S);
nlimgnC = 1-Kg*UdeC*V_C*cla/(498*W/S);
nlimgpD = 1+Kg*UdeD*V_D*cla/(498*W/S);
nlimgnD = 1-Kg*UdeD*V_D*cla/(498*W/S);
nlimgpB = 1+Kg*UdeB*V_B*cla/(498*W/S);
nlimgnB = 1-Kg*UdeB*V_B*cla/(498*W/S);
fvb = V_s*sqrt(nlimgpB)==V_B;
VBs = double(solve(fvb,V_B));
nlimgpB = double(subs(nlimgpB,V_B,VBs));
nlimgnB = double(subs(nlimgnB,V_B,VBs));
n_pos_rangeg = linspace(0,nlimgpB,30); %Vstall at index 12.
pos_stallg = double(subs(V_sn, nz, n_pos_rangeg));
%%
close all;
figure(1)
hold on;
grid on;
line([0,V_L],[0,0],'Color','black','LineStyle','-')
% draw 1g stall
line([pos_stall(12) pos_stall(12)], [0 1],'Color','red')
line([0 pos_stall(12)],[1 1],'Color','red','LineStyle','-.')
line([neg_stall(12) neg_stall(12)], [0 -1],'Color','red')
line([0 neg_stall(12)],[-1 -1],'Color','red','LineStyle','-.')
% draw design cruising speed
line([V_H V_H],[n_neg 0],'Color','blue','LineWidth',1);

% draw design manuevering speed 
line([pos_stall(end) pos_stall(end)],[0 n_pos],'Color','blue','LineWidth',1);

% draw dive speed line
line([V_L V_L],[-1 n_pos],'Color','black','LineWidth',3);

% draw neg load limit from Vc to Vd
line([V_H V_L],[n_neg, -1],'Color','black','LineWidth',3);

%draw n_pos and n_neg
line([pos_stall(end) V_L],[n_pos n_pos],'Color','black','LineWidth',3);
line([neg_stall(end) V_H],[n_neg n_neg],'Color','black','LineWidth',3);

% plot lift curves
plot(pos_stall,n_pos_range,'k-','LineWidth',3);
plot(neg_stall,-n_neg_range,'k-','LineWidth',3);

xlim([0,1.1*V_L]);
ylim([floor(n_neg),ceil(n_pos)]);
xlabel('Equivalent Airspeed (kts)');
ylabel('Load Factor (g)');
title('MILITARY AIRCRAFT V-N MANEUVER DIAGRAM')

figure(2)
hold on;
grid on;
line([0,V_D],[0,0],'Color','black','LineStyle','-')
% draw 1g stall
line([pos_stall(12) pos_stall(12)], [0 1],'Color','red')
line([0 pos_stall(12)],[1 1],'Color','red','LineStyle','-.')
% draw design cruising speed
line([VBs VBs],[nlimgpB 0],'Color','blue','LineWidth',1);
line([VBs VBs],[nlimgnB 0],'Color','blue','LineWidth',1);


line([V_C V_C],[nlimgpC 0],'Color','blue','LineWidth',1);
line([V_C V_C],[nlimgnC 0],'Color','blue','LineWidth',1);

line([V_D V_D],[nlimgpD 0],'Color','blue','LineWidth',1);
line([V_D V_D],[nlimgnD 0],'Color','blue','LineWidth',1);

% draw design manuevering speed 
line([0 VBs],[1 nlimgpB],'Color','black','LineWidth',1.5,'LineStyle','--');
line([0 VBs],[1 nlimgnB],'Color','black','LineWidth',1.5,'LineStyle','--');

line([0 V_C],[1 nlimgpC],'Color','black','LineWidth',1.5,'LineStyle','--');
line([0 V_C],[1 nlimgnC],'Color','black','LineWidth',1.5,'LineStyle','--');

line([0 V_D],[1 nlimgpD],'Color','black','LineWidth',1.5,'LineStyle','--');
line([0 V_D],[1 nlimgnD],'Color','black','LineWidth',1.5,'LineStyle','--');



% draw dive speed line

line([VBs V_C],[ nlimgpB nlimgpC],'Color','black','LineWidth',1.5,'LineStyle','--');
line([VBs V_C],[ nlimgnB nlimgnC],'Color','black','LineWidth',1.5,'LineStyle','--');
line([V_C V_D],[ nlimgpC nlimgpD],'Color','black','LineWidth',1.5,'LineStyle','--');
line([V_C V_D],[ nlimgnC nlimgnD],'Color','black','LineWidth',1.5,'LineStyle','--');
% draw neg load limit from Vc to Vd


%draw n_pos and n_neg


%draw cruise FAR 23

% plot lift curves
plot(pos_stallg,n_pos_rangeg,'k-','LineWidth',3);
xlim([0,1.1*V_D]);
ylim([floor(nlimgnB),ceil(nlimgpB)]);
xlabel('Equivalent Airspeed (kts)');
ylabel('Load Factor (g)');
title('FAR 25 AIRCRAFT V-N GUST DIAGRAM')
