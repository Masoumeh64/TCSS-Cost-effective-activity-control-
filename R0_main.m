function lambda = R0(A,B)
global activation_rate_infected inactivation_rate_infected
global activation_rate_asym inactivation_rate_asym
global beta_c beta_I delta eta_prime eta
global N PL P

activation_rate_I =  activation_rate_infected;
inactivation_rate_I =  inactivation_rate_infected;

t = 4*N;


% no_class1 = floor(N/3);
% no_class2 = floor(N/3);
% no_class3 = floor(N/3)+N-3*floor(N/3);
% 
% P1 = PL(1)*ones(floor(N/3),1);
% P2 = PL(2)*ones(floor(N/3),1);
% P3 = PL(3)*ones(floor(N/3)+N-3*floor(N/3),1);
% P = [P1;P2;P3];
% F = diag(diag(P));
% P = P-F;

R1_max = max(activation_rate_asym)*ones(N,1);

psi = max([2*R1_max'+eta_prime  inactivation_rate_asym+eta_prime ...
    activation_rate_I+delta  inactivation_rate_I+delta]);

psi = psi+0.01;
cnt = 0;
tic

cvx_begin gp
variables u(t) lambda
cvx_solver mosek

R1 = activation_rate_asym;
R_hat1 = psi-eta_prime-R1;
Y = R1+inactivation_rate_asym;

A_1 = diag(inactivation_rate_asym./Y)*A;
A_2 = diag(R1./Y)*A;
P_prime = diag(R1./Y)*((diag(P)*ones(N)*diag(P)).*B);

R_1 = diag(R1);
R_2 = diag(inactivation_rate_asym);
R_1_I = diag(activation_rate_I);
R_2_I = diag(inactivation_rate_I);

W = diag(R_hat1);

M = [beta_c*A_1+W        beta_c*A_1+R_2     beta_I*A_1     beta_I*A_1;...
    beta_c*A_2+R_1  beta_c*(A_2+P_prime)   beta_I*A_2     beta_I*(A_2+P_prime);...
    eta*eye(N)          zeros(N)               zeros(N)     R_2_I;...
    zeros(N)            eta*eye(N)             R_1_I    zeros(N)];
Dt = diag([zeros(1,N)     psi-inactivation_rate_asym-eta_prime ...
    psi-activation_rate_I-delta    psi-inactivation_rate_I-delta]);
M = M+Dt;


minimize( lambda )
subject to

M*u <= lambda*u;
u(1) <= 2;
1 <= u(1) ;

cvx_end
lambda = lambda-psi;
toc
