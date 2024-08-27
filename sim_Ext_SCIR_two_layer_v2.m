function [prevalence,z_mean] = sim_Ext_SCIR_two_layer_v2(no_run,y0)
global N  activation_rate_asym  inactivation_rate_asym activation_rate_infected inactivation_rate_infected
global p_link_const beta_c beta_I delta eta_prime eta kappa
global Nei_s Nei_d J1 J2 I1 I2 link t_final
%Pre-processing
%__________________________________________________________________________
% Constructing Neighbourhood vectors (input:Graph))
avg = 0;
%________________________________________________________________________
% Initial state and rate vectors
y0 = reshape(y0,8,N);
is_sus = (y0(1,:)+y0(2,:)==1);
is_carrier  = (y0(3,:)+y0(4,:)==1);
is_infected = (y0(5,:)+y0(6,:)==1);
is_recovered = (y0(7,:)+y0(8,:)==1);

is_active = (y0(2,:)+y0(4,:)+y0(6,:)+y0(8,:)==1);
is_inactive = (y0(1,:)+y0(3,:)+y0(5,:)+y0(7,:)==1);
rate_init = zeros(2,N);
link_active_init = zeros(2,N);
no_sus = sum(is_sus);
no_infected = sum(is_infected);
%___________________________________________________________________


%_____________________________________________________________________
% Main program
z1_mean = [];z2_mean = [];z3_mean = [];z4_mean = [];
ts_mean = [];
a1 = 0;a2 = 0; a3 = 0;
e1 = 0; e2 = 0; e3 = 0;
for count = 1:no_run
    
    if mod(count, 50) == 0
        count;
        [avg/count/N ];
    end
    is_sus = (y0(1,:)+y0(2,:)==1);
    is_carrier  = (y0(3,:)+y0(4,:)==1);
    is_infected = (y0(5,:)+y0(6,:)==1);
    is_recovered = (y0(7,:)+y0(8,:)==1);
    
    is_active = (y0(2,:)+y0(4,:)+y0(6,:)+y0(8,:)==1);
    is_inactive = (y0(1,:)+y0(3,:)+y0(5,:)+y0(7,:)==1);
    [rate_init,link_active_init] = initialization_SCIR_extended(is_sus,is_infected,is_carrier,is_active,is_inactive);
    rate = rate_init;
    link_active = link_active_init;
    no_sus = sum(is_sus);
    no_infected = sum(is_infected);
    no_carrier = sum(is_carrier);
    no_recovered = sum(is_recovered);
    cnt = 0;
    t = 0;
    z1 = [];z2 = [];z3 = [];z4 = [];
    ts = [];
    no_temp = zeros(1,N);
    % %     tic
    while (1)
        cnt = cnt+1;
        if mod(cnt,100) == 0
            [count t no_sus no_infected];
%             [e2/e1 e3/e1]
        end
        %____________________________________________________________
        R = sum(sum(rate,1));
        h= -log(rand)/R;
        ts(cnt)=t+h;
        t = t+h;
                z1(cnt) = no_sus;
                z2(cnt) = no_carrier;
                z3(cnt) = no_infected;
                z4(cnt) = no_recovered;
        
        %______________________________________________________________
        % Setting the new activity states and disease states
        
       
        node = rnd_draw(sum(rate,1));
        w1 = rnd_draw(rate(:,node))-1;

        is_active_new = w1*is_inactive(node)+(1-w1)*is_active(node);
        is_inactive_new = w1*is_active(node)+(1-w1)*is_inactive(node);
        is_active(node) = is_active_new;
        is_inactive(node) = is_inactive_new;
        
        ss = [eta eta_prime-eta];
        rr = [kappa 1-kappa];
        w2 = rnd_draw(ss)-1;
        w3 = rnd_draw(rr)-1;
        
       
        is_carrier_new = (1-w3)*is_sus(node);
        is_infected_new = w3*is_sus(node)+(1-w2)*is_carrier(node);
        is_recovered_new = w2*is_carrier(node)+is_infected(node);
        is_sus_new  = 0;
        is_sus_pre = is_sus(node);
        is_carrier_pre = is_carrier(node);
        
        is_sus(node) = (1-w1)*is_sus_new+w1*is_sus(node);
        is_carrier(node) = (1-w1)*is_carrier_new+w1*is_carrier(node);
        is_infected(node) = (1-w1)*is_infected_new+w1*is_infected(node);
        is_recovered(node) = (1-w1)*is_recovered_new+w1*is_recovered(node);
        [is_sus(node) is_carrier(node) is_infected(node) is_recovered(node)];
        
        
                        
        %________________________________________________________________
        % Setting the new neighbourhood vector according to activity states
        
        
        set_nodes = [];
        temp_neigh = [];
        %_________________________________________________________________
        % Acitivation
        if w1*is_active(node) == 1
           
            
            f = J1(node):J2(node);
            u = Nei_d(f);
            link_active(f) = (rand(1,length(u))<= is_active(u)*p_link_const(node));
            ind = find(link_active(f)==1);
            link_active(link(f(ind))) = 1;
                        
            temp_neigh = u(ind);
            per_neigh = Nei_s(I1(node):I2(node));
            
            rate(2,node) = (is_sus(node)+is_carrier(node))*inactivation_rate_asym(node)+...
                             is_infected(node)* inactivation_rate_infected(node);
            
            rn = sum(is_carrier(temp_neigh))*beta_c+sum(is_infected(temp_neigh))*beta_I;
            rate(1,node) = rate(1,node)+is_sus(node)*rn;
            
            rate(1,temp_neigh) = rate(1,temp_neigh)+is_sus(temp_neigh)*...
                                 (beta_c*is_carrier(node)+beta_I*is_infected(node));
            
                             

            %_________________________________________________________________
            % Inacitivation
            
        elseif w1*is_inactive(node) == 1
            f = J1(node):J2(node);
            u = Nei_d(J1(node):J2(node));
            ind = find(link_active(J1(node):J2(node))==1);
            link_active(link(f(ind))) = 0;
            link_active(J1(node):J2(node)) = 0;

            temp_neigh = u(ind);
           
            
            rate(2,node) = (is_sus(node)+is_carrier(node))*activation_rate_asym(node)+...
                             is_infected(node)* activation_rate_infected(node);
            
            rn = sum(is_carrier(temp_neigh))*beta_c+sum(is_infected(temp_neigh))*beta_I;            
            rate(1,node) = rate(1,node)-is_sus(node)*rn;
             
            rate(1,temp_neigh) = rate(1,temp_neigh)-is_sus(temp_neigh)*...
                                 (beta_c*is_carrier(node)+beta_I*is_infected(node));  
            

            %_________________________________________________________________
            % Change in disease state
        else
            u = Nei_d(J1(node):J2(node));
            ind = find(link_active(J1(node):J2(node))==1);
            temp_neigh = u(ind);
            per_neigh = Nei_s(I1(node):I2(node));
            
            rate(1,node) = is_carrier(node)*eta_prime+is_infected(node)*delta;
            
            rI = (1-is_carrier_pre)*beta_I+is_carrier_pre*(beta_I-beta_c);
            rR = (1-is_carrier_pre)*beta_I+is_carrier_pre*(beta_c);
            
            
            rate(1,temp_neigh) = rate(1,temp_neigh)+is_sus(temp_neigh')*...
                                 (is_carrier(node)*beta_c+is_infected(node)*rI-is_recovered(node)*rR);
                             
            rate(1,per_neigh) = rate(1,per_neigh)+is_sus(per_neigh')*...
                                 (is_carrier(node)*beta_c+is_infected(node)*rI-is_recovered(node)*rR);           
            

        end
       
        
        
        if ~isempty(find(rate<0))
            1;
        end
        
        
        
        
        %___________________________________________________________________
        no_sus = sum(is_sus);
        no_infected = sum(is_infected);
        no_carrier = sum(is_carrier);
        no_recovered = sum(is_recovered);
        no_active = sum(is_active);
        no_inactive = sum(is_inactive);
        if no_infected+no_carrier ==0
            break;
        end
       
    end
        for i = 1:t_final
            e = find(ts<=i);
            if isempty(e)
                z1_mean(i) = ((count-1)*z1_mean(i)+N-1)/count;
                z2_mean(i) = ((count-1)*z2_mean(i)+1)/count;
                z3_mean(i) = ((count-1)*z3_mean(i))/count;
                z4_mean(i) = ((count-1)*z4_mean(i))/count;
            else
                k = e(end);
                if count >1
                    z1_mean(i) = ((count-1)*z1_mean(i)+z1(k))/count;
                    z2_mean(i) = ((count-1)*z2_mean(i)+z2(k))/count;
                    z3_mean(i) = ((count-1)*z3_mean(i)+z3(k))/count;
                    z4_mean(i) = ((count-1)*z4_mean(i)+z4(k))/count;
                else
                    z1_mean(i) = z1(k);
                    z2_mean(i) = z2(k);
                    z3_mean(i) = z3(k);
                    z4_mean(i) = z4(k);
                end
            end
        end
%         z = [z1' z2' z3' z4'];
        z_mean = [z1_mean' z2_mean' z3_mean' z4_mean'];
%         [ts_mean,z1_mean,z2_mean,z3_mean,z4_mean] = time_average(count,ts_mean,ts,z_mean,z);
    avg = avg+no_recovered;
    %     toc
end
prevalence = avg/no_run/N;














