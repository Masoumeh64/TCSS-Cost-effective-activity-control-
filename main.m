clear all
clc
global N  activation_rate_asym  inactivation_rate_asym activation_rate_infected inactivation_rate_infected
global p_link_const beta_c beta_I delta eta_prime eta G
global act_rate_asym inact_rate_asym  act_rate_I  inact_rate_I p_link d
global  not_Nei J1 J2 I1_s I2_s Nei_s Nei_d
global l1 l2 kappa
N = 500;
act_rate_asym = 0.2;
inact_rate_asym = 0.2;
act_rate_I = 0.2;
inact_rate_I = 0.2;
p_link = 0.2;
activation_rate_asym = act_rate_asym*ones(1,N);
inactivation_rate_asym = inact_rate_asym*ones(1,N);
activation_rate_infected = act_rate_I*ones(1,N);
inactivation_rate_infected = inact_rate_I*ones(1,N);
p_link_const = p_link*ones(1,N);
beta_c = 0.1;
beta_I = 0.1;
delta = 1;
eta_prime = 0.8;
eta = 0.56;
d = 4;
l1 = 4;
l2 = 50;
kappa_vec = 0:0.1:1;
R_vec = [];
R0_1_vec = [];
sum_R = [];
for kappa = kappa_vec
[R,Rt] = reproduction_no_ext
R_vec = [R_vec R(1)];
R0_1_vec = [R0_1_vec Rt(1,1)];
sum_R = [sum_R sum(Rt(1,:))];
end
cnt = 0;
r_vec = 0;
% beta_c_vec = 0:0.005:0.15;
% beta_I_vec = 2*beta_c_vec;



% A = RandomRegularGraph(N, d);
% G1 = graph(A);
% A1 = ones(N)-eye(N)-A;
% G2 = graph(max(A1,0));
% save('graph_st_dy.mat','G1','G2');
load('graph_st_dy','G1','G2');

% Deriving neighbours


Nei_s =[];
not_Nei = [];
cnt = 1;
cnt1 = 1;
for i =1:N
    set_neigh = neighbors(G1,i);
    Nei_s = [Nei_s; neighbors(G1,i) ];
    e = neighbors(G2,i);
    not_Nei = [not_Nei; e];
    I1_s(i) = cnt;
    I2_s(i) = cnt+length(set_neigh)-1;
    J1(i) = cnt1;
    J2(i) = cnt1+length(e)-1;
    d(i) = length(set_neigh);
    cnt = cnt+length(set_neigh);
    cnt1 = cnt1+length(e);
end
Nei_d = not_Nei;

%________________________________________________________________________


% G = graph(A);
% plot(G)
% save('graph_static.mat','A')
% prevalence_MFE = zeros(1,length(r_vec));
% prevalence_sim = zeros(1,length(beta_c_vec));
% prevalence_ODE = zeros(1,length(beta_c_vec));
% fid_static = fopen('plot1.m','w');


vec = [0:0.001:0.2]
cnt = 0;
for act_rate_asym = vec
    cnt = cnt+1;
    [R,R0_1,R0_2,discrim,root1,root2] = reproduction_no;
    R_vec(cnt) = R;
    
end

cnt = 0;


% for beta_c = beta_c_vec
    
    cnt = cnt+1;
    
%     beta_I = 2*beta_c;
    
%     [beta_c R_vec(cnt)]
     no_run = 1000;
%     load('graph_static','A')
%     G = graph(A);
%     plot(G)
    %___________________________________________________________
    %                  RESULTS OF ANALYSIS
    %___________________________________________________________
    z = [];
    y0 = repmat([1;0; 0;0; 0;0; 0;0],N-1,1);
    y0 = [[0;0; 1;0; 0;0; 0;0];y0];
    [t,y] = ode45(@ODE_SCIR_two_layer,[0 80],y0);
    z(:,1) = sum(y(:,8*(0:N-1)+1),2);
    z(:,2) = sum(y(:,8*(0:N-1)+2),2);
    z(:,3) = sum(y(:,8*(0:N-1)+3),2);
    z(:,4) = sum(y(:,8*(0:N-1)+4),2);
    z(:,5) = sum(y(:,8*(0:N-1)+5),2);
    z(:,6) = sum(y(:,8*(0:N-1)+6),2);
    z(:,7) = sum(y(:,8*(0:N-1)+7),2);
    z(:,8) = sum(y(:,8*(0:N-1)+8),2);
    plot(t,z(:,1)+z(:,2),'-o')
    hold on
    plot(t,z(:,3)+z(:,4),'-x')
    hold on
    plot(t,z(:,5)+z(:,6),'-*')
    hold on
    plot(t,z(:,7)+z(:,8),'->')
    prevalence_ODE(cnt) = (z(end,7)+z(end,8))/N;
    
    %___________________________________________________________
    %                  SIMULATION RESULTS
    %___________________________________________________________
    
    y0 = repmat([0;1;0;0;0;0;0;0],N-50,1);
    y0 = [repmat([0;0;0;1;0;0;0;0],50,1);y0];
    no_run = 1000;
   [prevalence,z1_mean,z2_mean,z3_mean,z4_mean] = sim_SCIR_two_layer_v2(no_run,y0);
    
    %____________________________________________________________
    fprintf(fid_static,'%.4f  %.4f  %.4f  \n',[beta_c  R_vec(cnt) prevalence_ODE(cnt) ]);
    save('plot1_ODE.mat','gamma1_vec','gamma2_vec','prevalence_ODE')
    
% end
plot(beta_c_vec,prevalence_sim,'bo')
hold on
plot(beta_c_vec,prevalence_ODE,'g-')
hold on
vec = [];
for R0 = R1_vec
    fun = @(x) log((N-1)/x)-R0*(1-x/N);
    if R0 <= 1
        x0 = N;
    else
        x0 = N/2;
    end
    y = fsolve(fun,x0);
    vec = [vec 1-y/N];
end
plot(beta_c_vec,vec,'r--')


% load('result.mat','beta_c_vec','prevalence_MFE','prevalence_sim')
% beta_c_vec1 = beta_c_vec;
% prevalence_MFE1 = prevalence_MFE;
% prevalence_sim1 = prevalence_sim;
% load('result1.mat','beta_c_vec','prevalence_MFE','prevalence_sim')
% plot(beta_c_vec1,[prevalence_MFE1 prevalence_MFE])



