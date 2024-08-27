function [rate_init,link_active_init] = initialization_SCIR_v2(is_sus,is_infected,is_carrier,is_active,is_inactive)
global N  activation_rate_asym  inactivation_rate_asym activation_rate_infected inactivation_rate_infected
global p_link_const beta_c beta_I delta eta_prime eta
global Nei_s Nei_d I1 I2 J1 J2 d

%Initializing the active neighbors
link_active_init = zeros(1,length(Nei_d));
for node = 1:N
    if is_active(node) == 1
        u = Nei_d(J1(node):J2(node));
        link_active_init(J1(node):J2(node)) = link_active_init(J1(node):J2(node))+...
            (rand(1,length(u))<= is_active(u)*p_link_const(node));
        ind = find(link_active_init(J1(node):J2(node))==2);
        f = J1(node):J2(node);
        link_active_init(f(ind)) = 1;
        ind = find(link_active_init(J1(node):J2(node))==1);
        temp_neigh = u(ind);
        for i = 1:length(ind)
            nd = temp_neigh(i);
            f = J1(nd):J2(nd);
            m = find(Nei_d(f)==node);
            link_active_init(f(m)) = 1;
        end
        
    end
end
% Initializing the rates
rate_init = zeros(2,N);
for i = 1:N
    if is_sus(i)==1
        u = Nei_d(J1(i):J2(i));
        ind = find(link_active_init(J1(i):J2(i))==1);
        temp_neigh = u(ind);
        set_neigh = [Nei_s(I1(i):I2(i)); temp_neigh];
        rate_init(:,i) = [sum(is_carrier(set_neigh))*beta_c+sum(is_infected(set_neigh))*beta_I;...
            is_active(i)*inactivation_rate_asym(i)+is_inactive(i)*activation_rate_asym(i)]; 
    elseif is_carrier(i)==1
            rate_init(:,i) = [eta_prime;...
            is_active(i)*inactivation_rate_asym(i)+is_inactive(i)*activation_rate_asym(i)]; 
    elseif is_infected(i) == 1
            rate_init(:,i) = [delta;...
            is_active(i)*inactivation_rate_infected(i)+is_inactive(i)*activation_rate_infected(i)];
    end
end