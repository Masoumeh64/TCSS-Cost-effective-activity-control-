function [lambda,R1] = SGPA4_main(A,B,cost,gamma1_max,gamma1_min)
global activation_rate_infected inactivation_rate_infected
global inactivation_rate_asym
global beta_c beta_I delta eta_prime eta
global N P

activation_rate_I =  activation_rate_infected;
inactivation_rate_I =  inactivation_rate_infected;

t = 4*N;

R1_max = gamma1_max;  % maximum activation rate of asymptomatic users
R1_min = gamma1_min;
%____________________________________________________
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
%__________________________________________________
% for real dataset
% no_class1 = floor(N/2);
% no_class2 = N-no_class1;
% 
% P1 = PL(1)*ones(floor(N/3),1);
% P2 = PL(2)*ones(floor(N/3),1);
% % P3 = PL(3)*ones(floor(N/3)+N-3*floor(N/3),1);
% P = [P1;P2;P3];
%_____________________________________________________

psi = max([2*R1_max'+eta_prime   inactivation_rate_asym+eta_prime ...
         activation_rate_I+delta   inactivation_rate_I+delta]);

psi= psi+0.01;     
     
R1_int = R1_max-0.01;
R_hat1_int = psi-eta_prime-R1_int;

R1 = zeros(N,1); % this is in fact the vec of gamma^1_i 
% 
R_hat1 = zeros(N,1); % this is the vector of \hat{gamma}^1_i
% q1 = R1_int./(R1_int+inactivation_rate_asym);
% q2 = inactivation_rate_asym./(R1_int+inactivation_rate_asym);
% %
% r1 = R1_int./(R1_int+R_hat1_int);
% r2 = R_hat1_int./(R1_int+R_hat1_int);



cnt = 0;
tic
while sum(abs([R1_int R_hat1_int]-[R1 R_hat1]))/sum([R1_int R_hat1_int]) > 0.001
    cnt = cnt+1;
    %     norm(x0-x)
    if cnt >1
        R1_int = R1;
        R_hat1_int = R_hat1;
        
    end
    
    % R1 is the vector of activity Rs of users
    % R = R1+ inact_R_asym
    % R_hat1 = psi-eta_prime-R1
    
    cvx_begin gp
    variables R1(N) R_hat1(N)  u(t) lambda
    cvx_solver mosek
%     R1 = [X1(1)*ones(no_class1,1);X1(2)*ones(no_class2,1);X1(3)*ones(no_class3,1)];
%     R_hat1 = [X_hat1(1)*ones(no_class1,1);X_hat1(2)*ones(no_class2,1);X_hat1(3)*ones(no_class3,1)];
    
    % We are going to approximate R1+inact_R_asym at initial point
    q1 = R1_int./(R1_int+inactivation_rate_asym');
    q2 = inactivation_rate_asym'./(R1_int+inactivation_rate_asym');
    Y = (R1./q1).^q1 .* (inactivation_rate_asym'./q2).^q2;
    
    A_1 = diag(inactivation_rate_asym'./Y)*A;
    A_2 = diag(R1./Y)*A;
%     P_prime = diag(R1./Y)*diag(P)*B;
    P_prime = diag(R1./Y)*((diag(P)*ones(N)*diag(P)).*B);
    
    
    R_1 = diag(R1);
    R_2 = diag(inactivation_rate_asym);
    R_1_I = diag(activation_rate_I);
    R_2_I = diag(inactivation_rate_I);
    
    %     P_prime = diag(0.5*ones(1,N))*P;
    
    W = diag(R_hat1);
    
    M = [beta_c*A_1+W        beta_c*A_1+R_2     beta_I*A_1     beta_I*A_1;...
        beta_c*A_2+R_1  beta_c*(A_2+P_prime)   beta_I*A_2     beta_I*(A_2+P_prime);...
        eta*eye(N)          zeros(N)               zeros(N)     R_2_I;...
        zeros(N)            eta*eye(N)             R_1_I    zeros(N)];
    %
    Dt = diag([zeros(1,N)     psi-inactivation_rate_asym-eta_prime ...
         psi-activation_rate_I-delta    psi-inactivation_rate_I-delta]);
    %     Qt = M+Dt;
    M = M+Dt;
    
    
    %      W = Qt*u;
    %     qt = [R1' rand(1,3*N)];
    %     Qt= [qt;rand(4*N-1,4*N)];
    %     Qt= rand(4*N,4*N);
    
    %
    %     % We are going to approximate R1+R_hat1 at initial point
    r1 = R1_int./(R1_int+R_hat1_int);
    r2 = R_hat1_int./(R1_int+R_hat1_int);
    %     Qt = M;
    minimize( lambda )
    subject to
   
    M*u <= lambda*u;
    u(1) <= 2;
    1 <= u(1) ;
    
    sum(1./R1) <= cost;
    R1_min <= R1 ;
    R1 <= R1_max;
    
    (psi-eta_prime)*ones(N,1) == (R1./r1).^r1 .* (R_hat1./r2).^r2;
    [R1;R_hat1] <= 1.1*[R1_int;R_hat1_int]
    (1/1.1)*[R1_int;R_hat1_int] <= [R1;R_hat1]
    
    cvx_end
    
%     gamma1 =[R1(1) R1(no_class1+1) R1(no_class1+no_class2+1)]
%      [R_hat1(1) R_hat1(no_class1+1) R_hat1(no_class1+no_class2+1)]
%     sum(abs(R1_int-R1))/sum(R1_int) > 0.01
    sum(1./R1)
    sum(abs([R1_int R_hat1_int]-[R1 R_hat1]))/sum([R1_int R_hat1_int])
    [lambda-psi]
    cost
   

    
    toc
    
end% optimization variables
lambda = lambda-psi;