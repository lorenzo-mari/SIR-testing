clearvars
close all
clc

%% parameters
% free parameters
N=1e6; % baseline population size
R0=3; % basic reproduction number
gamma=1/7; % recovery rate (1/day)
kappa=1/3; % testing rate for infected (1/day)
epsilon=1/60; % testing rate for non-infected (1/day)
theta=1/10; % isolation release rate (1/day)
tpr=0.85; % test sensitivity
a=2.681; % susceptible-infected average separation in test results
b=1; % susceptible-infected variance ratio in test results

% other parameters
beta=R0*gamma; % transmission rate (1/person/day)
tnr=normcdf(a-b*norminv(tpr)); % test specificity
fnr=1-tpr; % false negative rate
fpr=1-tnr; % false positive rate

% simulation parameters
ndays=180; % simulation timespan (days)
I0=1; % initial number of infected people
ctrl_thresh=100/1e5; % threshold for control (new cases in the past 7 days)
twL=90; % time window length for intervention assessment (days)

%% control reproduction number
Rc=R0*(gamma*(gamma+theta*fnr)+kappa*tpr*theta*fnr)/(gamma+kappa*tpr)/(gamma+theta*fnr);
display(Rc)

%% simulation without controls
tspan_noctrl=0:ndays;
x0=[N-I0 I0 0 I0];
ode=@(t,x) SIR_ode(t,x,beta,gamma);
[~,x_noctrl]=ode45(ode,tspan_noctrl,x0);

%% two-step simulation (without/with controls)
Delta=tspan_noctrl(find(movsum([I0 diff(x_noctrl(:,4))'],[7 0])/N>=ctrl_thresh,1,'first')); % time lag for control deployment (days)

% before detection
tspan1=0:1:Delta;
ode=@(t,x) SIR_ode(t,x,beta,gamma);
[~,x1_temp]=ode45(ode,tspan1,x0);

% after detection
x0=zeros(1,12);
x0([1 4 7 10])=x1_temp(end,:);
tspan2=Delta:ndays;
ode=@(t,x) SIRtest_ode(t,x,beta,gamma,epsilon,kappa,theta,tpr,tnr,fpr,fnr);
[~,x2]=ode45(ode,tspan2,x0);

% assembling the simulation
tspan=[tspan1,tspan2(2:end)];
x1=zeros(Delta+1,12);
x1(:,[1 4 7 10])=x1_temp;
x=[x1;x2(2:end,:)];

%% post-processing
Sn=x(:,1); % susceptible individuals who never tested positive
Sp=x(:,2); % susceptible individuals who tested positive and are currently isolated
Sc=x(:,3); % susceptible individuals who have been released from isolation after their recovery has been certified via testing
In=x(:,4); % infected individuals who never tested positive
Ip=x(:,5); % infected individuals who tested positive and are currently isolated
Ic=x(:,6); % infected individuals who have been released from isolation after their recovery has been certified via testing
Rn=x(:,7); % recovered individuals who never tested positive
Rp=x(:,8); % recovered individuals who tested positive and are currently isolated
Rc=x(:,9); % recovered individuals who have been released from isolation after their recovery has been certified via testing
cumC=x(:,10); % cumulated cases
cumP_inf=x(:,11); % cumulated isolation mandates (infected individuals)
cumP_noninf=x(:,12); % cumulated isolation mandates (non-infected individuals)
cumP_noctrl=x_noctrl(:,4); % cumulated cases (without controls)

% daily incidence values
dailyC=[I0;diff(cumC)]; % new cases
dailyP_inf=[0;diff(cumP_inf)]; % new isolation mandates (infected individuals)
dailyP_noninf=[0;diff(cumP_noninf)]; % new isolation mandates (non-infected individuals)
dailyP=dailyP_inf+dailyP_noninf; % new isolation mandates (total)

% prevalence values
noninf_comm=Sn+Sc+Rn+Rc; % non-infected individuals in the community
inf_comm=In+Ic; % infected individuals in the community
tot_comm=noninf_comm+inf_comm; % total abundance of individuals in the community
noninf_isol=Sp+Rp; % non-infected individuals in isolation
inf_isol=Ip; % infected individuals in isolation
tot_isol=noninf_isol+inf_isol; % total abundance of individuals in isolation
inf_tot=inf_comm+inf_isol; % total abundance of infected individuals

% case count and isolation person-time
tw=1+[Delta Delta+twL];
tc=cumC(tw(2)-tw(1))/N*100; % total case count (% population)
display(tc)
cr=(1-(cumC(tw(2))-cumC(tw(1)))/(cumP_noctrl(tw(2))-cumP_noctrl(tw(1))))*100; % case-count reduction (% with respect to no control)
display(cr)
tq=trapz(tspan(tw(1):tw(2)),tot_isol(tw(1):tw(2)))/(N*twL)*100; % total isolation person-time (% total person-time)
display(tq)
Sp=trapz(tspan(tw(1):tw(2)),noninf_isol(tw(1):tw(2)))/(N*twL)*100; % superfluous isolation person-time (% total person-time)
display(Sp)

%% effective reproduction number
Rt=NaN(ndays+1,1);
for t=1:Delta
    Rt(t)=beta/gamma*Sn(t)*(Sn(t)+Rn(t))/(Sn(t)+In(t)+Rn(t))^2;
end
for t=Delta+1:ndays+1 
    Rt(t)=beta/gamma*(Sn(t)+Sc(t)+Rn(t)+Rc(t))/(Sn(t)+Sc(t)+In(t)+Ic(t)+Rn(t)+Rc(t))^2*...
        (Sn(t)*(gamma*(gamma+theta*fnr)+kappa*tpr*theta*fnr)+Sc(t)*(gamma+kappa*tpr)*(gamma+theta*fnr))/...
        (gamma+kappa*tpr)/(gamma+theta*fnr);
end

%% plots
dailyP_inf(1:Delta+1)=NaN;
dailyP_noninf(1:Delta+1)=NaN;
tot_isol(1:Delta+1)=NaN;
inf_isol(1:Delta+1)=NaN;
noninf_isol(1:Delta+1)=NaN;

figure('Renderer','painters')
tiledlayout(2,2,'TileSpacing','compact')
labels={'(a)','(b)','(c)','(d)'};
j=0;

j=j+1; 
nexttile;
yyaxis left
plot(tspan,dailyC./N*100,tspan,dailyP./N*100,'linewidth',1)
box off
set(gca,'XLim',[0 ndays],'TickDir','out','TickLabelInterpreter','latex','FontSize',8)
xlabel('Time $t$ (days)','Interpreter','latex','FontSize',10)
ylabel({'Daily incidence';'(\% population/day)'},'Interpreter','latex','FontSize',10)
yyaxis right
plot(tspan,Rt,'linewidth',1)
box off
ylabel({'Effective reproduction';'number, $\mathcal{R}_t$'},'Interpreter','latex','FontSize',10)
patch([0 Delta Delta 0],[min(get(gca,'YLim')) min(get(gca,'YLim')) max(get(gca,'YLim')) max(get(gca,'YLim'))],'k','EdgeColor','n','FaceAlpha',1/5)
hl=legend('Infected','Isolated');
set(hl,'Box','on','Interpreter','latex','FontSize',8,'Location','northwest','EdgeColor','none')
hl.ItemTokenSize=[15,18];
hl.BoxFace.ColorType='truecoloralpha';
hl.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
title(labels{j},'Interpreter','latex','FontSize',10)

j=j+1; 
nexttile;
plot(tspan,inf_tot./N*100,tspan,inf_comm./N*100,tspan,inf_isol./N*100,'linewidth',1)
box off
set(gca,'XLim',[0 ndays],'TickDir','out','TickLabelInterpreter','latex','FontSize',8) 
xlabel('Time $t$ (days)','Interpreter','latex','FontSize',10)
ylabel({'Infection prevalence';'(\% population)'},'Interpreter','latex','FontSize',10)
patch([0 Delta Delta 0],[min(get(gca,'YLim')) min(get(gca,'YLim')) max(get(gca,'YLim')) max(get(gca,'YLim'))],'k','EdgeColor','n','FaceAlpha',1/5)
hl=legend('Total','Community','Isolation');
set(hl,'Box','on','Interpreter','latex','FontSize',8,'Location','northwest','EdgeColor','none')
hl.ItemTokenSize=[15,18];
hl.BoxFace.ColorType='truecoloralpha';
hl.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
title(labels{j},'Interpreter','latex','FontSize',10)

j=j+1; 
nexttile;
yyaxis left
plot(tspan,inf_comm./tot_comm*100,'linewidth',1)
set(gca,'XLim',[0 ndays],'TickDir','out','TickLabelInterpreter','latex','FontSize',8) 
xlabel('Time $t$ (days)','Interpreter','latex','FontSize',10)
ylabel({'Infection prevalence';' (\% community)'},'Interpreter','latex','FontSize',10)
title(labels{j},'Interpreter','latex','FontSize',10)
yyaxis right
plot(tspan,inf_isol./tot_isol*100,'linewidth',1)
ylabel({'Infection prevalence';'(\% isolated)'},'Interpreter','latex','FontSize',10)
box off
patch([0 Delta Delta 0],[min(get(gca,'YLim')) min(get(gca,'YLim')) max(get(gca,'YLim')) max(get(gca,'YLim'))],'k','EdgeColor','n','FaceAlpha',1/5)

j=j+1; 
nexttile;
plot(tspan,tot_isol./N*100,tspan,inf_isol./N*100,tspan,noninf_isol./N*100,'linewidth',1)
box off
set(gca,'XLim',[0 ndays],'TickDir','out','TickLabelInterpreter','latex','FontSize',8)
patch([0 Delta Delta 0],[min(get(gca,'YLim')) min(get(gca,'YLim')) max(get(gca,'YLim')) max(get(gca,'YLim'))],'k','EdgeColor','n','FaceAlpha',1/5)
hl=legend('Total','Infected','Non infected');
set(hl,'Box','on','Interpreter','latex','FontSize',8,'Location','northwest','EdgeColor','none')
hl.ItemTokenSize=[15,18];
hl.BoxFace.ColorType='truecoloralpha';
hl.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
xlabel('Time $t$ (days)','Interpreter','latex','FontSize',10)
ylabel({'Isolation prevalence';'(\% population)'},'Interpreter','latex','FontSize',10)
title(labels{j},'Interpreter','latex','FontSize',10)

%% ODE functions
function d_x_dt=SIR_ode(~,x,beta,gamma)
    
    % state variables
    S=x(1); % susceptible individuals
    I=x(2); % infected individuals
    R=x(3); % recovered individuals

    % force of infection
    lambda=beta*I/(S+I+R);

    % ODEs
    d_S_dt=-lambda*S; % susceptible individuals 
    d_I_dt=lambda*S-gamma*I; % infected individuals
    d_R_dt=gamma*I; % recovered individuals
    d_C_dt=lambda*S; % cumulated cases

    % output
    d_x_dt=[d_S_dt; d_I_dt; d_R_dt; d_C_dt];
end

function d_x_dt=SIRtest_ode(~,x,beta,gamma,epsilon,kappa,theta,tpr,tnr,fpr,fnr)
    
    % state variables
    Sn=x(1); % susceptible individuals who never tested positive
    Sp=x(2); % susceptible individuals who tested positive and are currently isolated
    Sc=x(3); % susceptible individuals who have been released from isolation after their recovery has been certified via testing
    In=x(4); % infected individuals who never tested positive
    Ip=x(5); % infected individuals who tested positive and are currently isolated
    Ic=x(6); % infected individuals who have been released from isolation after their recovery has been certified via testing
    Rn=x(7); % recovered individuals who never tested positive 
    Rp=x(8); % recovered individuals who tested positive and are currently isolated
    Rc=x(9); % recovered individuals who have been released from isolation after their recovery has been certified via testing

    % force of infection
    lambda=beta*(In+Ic)/(Sn+Sc+In+Ic+Rn+Rc);

    % ODEs
    d_Sn_dt=-(lambda+epsilon*fpr)*Sn; % susceptible individuals who never tested positive
    d_Sp_dt=epsilon*fpr*Sn-theta*tnr*Sp; % susceptible individuals who tested positive and are currently isolated
    d_Sc_dt=theta*tnr*Sp-lambda*Sc; % susceptible individuals who have been released from isolation after their recovery has been certified via testing
    d_In_dt=lambda*Sn-(gamma+kappa*tpr)*In; % infected individuals who never tested positive
    d_Ip_dt=kappa*tpr*In-(gamma+theta*fnr)*Ip; % infected individuals who tested positive and are currently isolated
    d_Ic_dt=lambda*Sc+theta*fnr*Ip-gamma*Ic; % infected individuals who have been released from isolation after their recovery has been certified via testing
    d_Rn_dt=gamma*In-epsilon*fpr*Rn; % recovered individuals who never tested positive
    d_Rp_dt=gamma*Ip+epsilon*fpr*Rn-theta*tnr*Rp; % recovered individuals who tested positive and are currently isolated
    d_Rc_dt=gamma*Ic+theta*tnr*Rp; % recovered individuals who have been released from isolation after their recovery has been certified via testing
    d_C_dt=lambda*(Sn+Sc); % cumulated cases
    d_P_inf_dt=kappa*tpr*In; % cumulated isolation mandates (infected individuals)
    d_P_noninf_dt=epsilon*fpr*(Sn+Rn); % cumulated isolation mandates (non-infected individuals)

    % output
    d_x_dt=[d_Sn_dt; d_Sp_dt; d_Sc_dt; d_In_dt; d_Ip_dt; d_Ic_dt;d_Rn_dt; d_Rp_dt; d_Rc_dt; d_C_dt; d_P_inf_dt; d_P_noninf_dt];
end