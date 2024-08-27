function [R,R0_1,R0_2,discrim,root1,root2] = reproduction_no
global N  activation_rate_asym  inactivation_rate_asym activation_rate_infected inactivation_rate_infected
global p_link_const beta_c beta_I delta eta_prime eta
global act_rate_asym inact_rate_asym  act_rate_I  inact_rate_I p_link 
d1 = 15;
qs1 = inact_rate_asym/(inact_rate_asym+act_rate_asym);
qs2 = act_rate_asym/(inact_rate_asym+act_rate_asym);
% qs1 = 1;
% qs2 = 0;

F1 =[d1*qs1  d1*qs1;d1*qs2  (N*p_link+d1)*qs2];
F = [beta_c*F1  beta_I*F1; zeros(2,2) zeros(2,2)];

V = [act_rate_asym+eta_prime  -inact_rate_asym 0 0;...
     -act_rate_asym           inact_rate_asym+eta_prime  0  0;...
     -eta    0             act_rate_I+delta  -inact_rate_I ;...
     0      -eta            -act_rate_I      inact_rate_I+delta  ];
R = max(abs(eig(F*V^(-1))));



V1 =  [act_rate_asym+eta_prime  -inact_rate_asym ;...
     -act_rate_asym           inact_rate_asym+eta_prime ];
 
 V2 = [  act_rate_I+delta  -inact_rate_I ;...
      -act_rate_I      inact_rate_I+delta  ];
  
H = F1*(beta_c*eye(2)+eta*beta_I*V2^(-1))*V1^(-1);  
trace(H) 
det(H)
R1 = max(abs(eig(H)));


x = eta_prime+act_rate_asym+inact_rate_asym;
y = delta+act_rate_I+inact_rate_I;
u = beta_c/eta+beta_I/delta;
w = beta_c/eta+beta_I/y;
h1 = N*p_link*eta*inact_rate_asym*w/eta_prime;
h2 = N*p_link*d1*eta^2*inact_rate_asym*w*u/eta_prime;

E = N*p_link*eta*(act_rate_asym+eta_prime)*w/eta_prime/x+...
    N*p_link*eta*beta_I*act_rate_I/eta_prime/y/delta;
    
    
trace_L = d1*eta*u/eta_prime+qs2*E;
det_L = N*p_link*d1*eta^2*qs1*qs2*u*w/eta_prime/x;




R0_1 = d1*eta*(beta_c/eta+beta_I/delta)/eta_prime;
R0_2 = p_link*N*eta*(beta_c/eta+beta_I*(delta+act_rate_I)/delta/y)/eta_prime; 
h3 = R0_1+R0_2-1;
trace_L = -h1/x*qs2+R0_1+R0_2*qs2;
% trace_L = -h1/x+R0_1+R0_2;
det_L = h2*(x-eta_prime-inact_rate_asym)/x/(x-eta_prime)^2;
% det_L = h2/x/(x-eta_prime);

[R0_2-h1/x R0_1*qs1*N*p_link*eta*w/x]
[eta_prime+act_rate_asym R0_1*eta_prime*qs1]
R0_1+qs2*(R0_2-h1/x-R0_1*qs1*N*p_link*eta*w/x)
if x-eta_prime == 0
    det_L = N*p_link*d1*eta^2*qs1*w*u/eta_prime/x;
end
Delta = trace_L^2-4*det_L;
discrim = (eta_prime*h3+h1)^2-4*h3*(h1*eta_prime-h2);
root1 = ((eta_prime*h3+h1)+sqrt(discrim))/2/h3;
root2 = ((eta_prime*h3+h1)-sqrt(discrim))/2/h3;
R2 = (trace_L+sqrt(Delta))/2;

%  R2 = (trace_L+sqrt(Delta))/2;
% R2 = N*p_link*eta*w/x;



