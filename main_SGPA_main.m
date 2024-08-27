clear all
clc
global act_rate_I inact_rate_I
global inact_rate_asym
global beta_c beta_I delta eta_prime eta
global N PL
global activation_rate_asym inactivation_rate_asym
global activation_rate_infected inactivation_rate_infected 
global p_link_const
global Nei_s not_Nei I1_s I2_s J1 J2 d P

%___________________________________________________________________
% Parameters
N = 200;
d = 4;
PL = [0.1 0.2 0.8];
% % PL = [0 0 0];
% 
beta_c = 0.15;
beta_I = 0.2;
delta = 0.2;
eta_prime = 0.8;
eta = 0.56;
act_rate_I = 0.2;
inact_rate_I = 0.8;
inact_rate_asym = 0.2;


% beta_c = 0.01;
% beta_I = 0.03;
% delta = 0.5;
% eta_prime = 0.8;
% eta = 0.56;
% act_rate_I = 0;
% inact_rate_I = 0.8;
% inact_rate_asym = 0.2;


% beta_I = 0.03;
% beta_c = 0.7*beta_I;
% 
% delta = 0.18;
% eta_prime = 0.14;
% eta = 0.098;
% act_rate_I = 0;
% inact_rate_I = 0.8;
% inact_rate_asym = 0.1;
% 


%___________________________________________________________________
% Graph and neighborhood

A = RandomRegularGraph(N, d);
G1 = graph(A);
%ER
% p = 0.2;
% B1 = rand(N,N)<p;
% B = triu(B1,1)+triu(B1,1).';
% B = B-A;
% B = max(B,0);
% G2 = graph(B);
%BA
B = BAgraph(N,20,10);
B = B-A;
B = max(B,0);
G2 = graph(B);

%______________________________________________________________________
% Real dataset

% Filep='Permanent_Layer.txt';
% Filet='Temporal_Layer.txt';
% Lp=load(Filep); %L1=L(:,1); L2=L(:,2); L3=L(:,3);
% Lt=load(Filet);
% EdgeList_p = unique(sort(Lp(:,1:2),2),'rows');
% EdgeList_t = unique(sort(Lt(:,1:2),2),'rows');
% H1 = graph(EdgeList_p(:,1),EdgeList_p(:,2));
% H2 = graph(EdgeList_t(:,1),EdgeList_t(:,2));
% bins = conncomp(H1,'OutputForm','cell');
% G1 = subgraph(H1,bins{4});
% G2 = subgraph(H2,bins{4});
% A = adjacency(G1);
% B = adjacency(G2);
% B = B-A;
% B = max(B,0);
% G2 = graph(B);
% N = length(bins{4});

%_______________________________________________________________________

% % A = RandomRegularGraph(N, d);
% % G1 = graph(A);
% % % 
% for i = 1:N
%     set_nei = neighbors(G2,i);
%     no_nei = length(set_nei);
%     pot_nei = setdiff(i:N,set_nei);
%     m = rand(1,length(pot_nei))<0.05;
%     B(i,pot_nei) = m;
%     B(pot_nei,i) = m;
% end
% G2 = graph(B);
% 
%  p = 0.2;
% B1 = rand(N,N)<p;
% B = triu(B1,1)+triu(B1,1).';
% B = B-A;
% B = max(B,0);
% G2 = graph(B);

%___________________________________________________________________
% Classes of activity
N1 = floor(N/6);
N2 = floor(2*N/3);
N3 = N-N1-N2;
P1 = PL(1)*ones(floor(N/3),1);
P2 = PL(2)*ones(floor(N/3),1);
P3 = PL(3)*ones(floor(N/3)+N-3*floor(N/3),1);


P1 = PL(1)*ones(N1,1);

P2 = PL(2)*ones(N2,1);

P3 = PL(3)*ones(N3,1);


P = [P1;P2;P3];
p_link_const = P;

%for real dataset 
% N1 = floor(N/2);
% N2 = N-N1;
% PL(1) = 1;
% PL(2) = 1;
% p_link_const = [PL(1)*ones(1,N1) PL(2)*ones(1,N2) ];
% P = p_link_const;
%___________________________________________________________________

E1 = G1.Edges.EndNodes;
E2 = G2.Edges.EndNodes;
E = [E1;E2];
% W1 = ones(size(E1,1),1);
W1 = G1.Edges.Weight;
W2 = p_link_const(E2(:,1))'.*p_link_const(E2(:,2))';
W = [W1;W2'];
G = graph(E(:,1),E(:,2),W);

% G.Edges.Weight = W;


%_________________________________________________________
inactivation_rate_asym = inact_rate_asym*ones(1,N);
activation_rate_infected = act_rate_I*ones(1,N);
inactivation_rate_infected = inact_rate_I*ones(1,N);


%___________________________________________________________________
% Cost
gamma1_max = 0.3*ones(N,1); 
gamma1_min = 0.08*ones(N,1);
cost_max = sum(1./gamma1_min);
cost_min = sum(1./gamma1_max);
step = (cost_max-cost_min)/10;
cost_vec = cost_min:step:cost_max;
cost_vec = cost_vec(4);

% Cost for real dataset

% gamma1_max = 3*ones(N,1); 
% gamma1_min = 0.01*ones(N,1);
% cost_max = sum(1./gamma1_min);
% cost_min = sum(1./gamma1_max);
% step = (cost_max-cost_min)/10;
% cost_vec = cost_min:step:cost_max;
% cost_vec = cost_vec(3);



%Neighborhood
J1 = zeros(1,N);
J2 = zeros(1,N);
Nei_s =[];
Nei_d = [];
cnt = 1;
cnt1 = 1;
d2 = zeros(1,N);
d1 = zeros(1,N);
for i =1:N
    set_neigh = neighbors(G1,i);
    Nei_s = [Nei_s; neighbors(G1,i) ];
    e = neighbors(G2,i);
    Nei_d = [Nei_d; e];
    I1(i) = cnt;
    I2(i) = cnt+length(set_neigh)-1;
    d1(i) = length(set_neigh);
    if ~isempty(e)
    J1(i) = cnt1;
    J2(i) = cnt1+length(e)-1;
    d2(i) = length(e);
    end
    cnt = cnt+length(set_neigh);
    cnt1 = cnt1+length(e);
end
%_____________________________________________________________________
% Main

gamma1_vec = zeros(N,length(cost_vec));
cnt = 0;
lambda1_vec = [];
lambda2_vec = [];
lambda_eigen_vec = [];

for cost = cost_vec
    
%     activation_rate_asym = 0.01*ones(1,N);
%     lambda_first_layer = R0(A,B);
%     
    
    cnt = cnt+1
    [lambda1,gamma1] = SGPA4_main(A,B,cost,gamma1_max,gamma1_min);
    gamma1_vec(:,cnt) = gamma1;
%     
    degree_vec = d1'+sum(diag(p_link_const)*B*diag(p_link_const),2);
    [deg,deg_ind] = sort(degree_vec);
    deg = fliplr(deg);
    deg_ind = fliplr(deg_ind);
    activation_rate_asym = ...
        high_degree_assignment_general(cost,deg,deg_ind,gamma1_min,gamma1_max);
    lambda2 = R0(A,B);
    
    
    eigen_vec = centrality(G,'closeness','Cost',G.Edges.Weight);
    [eigen,eigen_ind] = sort(eigen_vec);
    eigen = fliplr(eigen);
    eigen_ind = fliplr(eigen_ind);
    activation_rate_asym = ...
        high_degree_assignment_general(cost,eigen,eigen_ind,gamma1_min,gamma1_max);
    lambda_eigen = R0(A,B);
    
    
    
    lambda1_vec = [lambda1_vec lambda1];
    lambda2_vec = [lambda2_vec lambda2];
    lambda_eigen_vec = [lambda_eigen_vec lambda_eigen];
end
% 0.0191    0.5902    1.1613    1.7324
% 
%   Columns 5 through 8
% 
%     2.3035    2.8746    3.4456    4.0167
% 
%   Columns 9 through 11
% 
%     4.5878    5.1589    5.7300

plot(eigen_vec(1:N1),gamma1(1:N1),'bx')
hold on
plot(eigen_vec(N1+1:N1+N2),gamma1(N1+1:N1+N2),'ro')
hold on
plot(eigen_vec(N1+N2+1:N),gamma1(N1+N2+1:N),'g>')

plot(cost_vec,gamma1_vec(1,:),'-*')
hold on
plot(cost_vec,gamma1_vec(2,:),'-o')
hold on
plot(cost_vec,gamma1_vec(3,:),'-x')


% cost_vec =  (1.0e+03) *[0.3333    0.4250    0.5167    0.6083    0.7000 ...
% 0.7917    0.8833    0.9750    1.0667    1.1583  1.2500];


% lambda1_vec = [1.1166    0.9331    0.8010    0.6987    0.6137 ...
%    0.5417  0.4800  0.4268  0.3836    0.3477    0.3220];
% 
% lambda2_vec = [ 1.1165    0.9366   0.8271   0.7550  0.6481...
%   0.5590    0.5020    0.4459    0.3844    0.3469  0.3217];
% 
% plot(cost_vec,lambda1_vec,'*--',cost_vec,lambda2_vec,'--o')
plot(degree_vec(1:N1),gamma1(1:N1),'bx')
hold on
plot(degree_vec(N1+1:N1+N2),gamma1(N1+1:N1+N2),'ro')
hold on
plot(degree_vec(N1+N2+1:N),gamma1(N1+N2+1:N),'g>')

