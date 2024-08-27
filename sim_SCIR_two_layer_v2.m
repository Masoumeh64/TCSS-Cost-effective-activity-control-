function [prevalence,z1_mean,z2_mean,z3_mean,z4_mean] = sim_SCIR_two_layer_v2(no_run,y0)
global N  activation_rate_asym  inactivation_rate_asym activation_rate_infected inactivation_rate_infected
global p_link_const beta_c beta_I delta eta_prime eta
global Nei_s Nei_d I1 I2 J1 J2 d
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
for count = 1:no_run
    
        if mod(count, 100) == 0
    %[a1/count a2/count a3/count]
    count
    [avg/count/N ]
        end
    is_sus = (y0(1,:)+y0(2,:)==1);
    is_carrier  = (y0(3,:)+y0(4,:)==1);
    is_infected = (y0(5,:)+y0(6,:)==1);
    is_recovered = (y0(7,:)+y0(8,:)==1);
    
    is_active = (y0(2,:)+y0(4,:)+y0(6,:)+y0(8,:)==1);
    is_inactive = (y0(1,:)+y0(3,:)+y0(5,:)+y0(7,:)==1);
    [rate_init,link_active_init] = initialization_SCIR_v2(is_sus,is_infected,is_carrier,is_active,is_inactive);
    rate = rate_init;
%     for h = 1:N
%         u = Nei_d(J1(h):J2(h));
%         ind1 = find(link_active_init(J1(h):J2(h))==1);
%         temp_neigh1 = u(ind1);
%         per_neigh1 = Nei_s(I1(h):I2(h));
%         neigh1 = [per_neigh1;temp_neigh1];
%         if is_sus(h)==1
%             
%             if rate(1,h) ~= sum(is_carrier(neigh1))*beta_c+sum(is_infected(neigh1))*beta_I
%                 h;
%                 1;
%             end;
%         elseif is_carrier(h)==1
%             if rate(1,h) ~= eta_prime
%                 h;
%                 1;
%             end
%         end
%     end
%     if ~isempty(find(rate<0))
%         1;
%     end
    link_active = link_active_init;
    no_sus = sum(is_sus);
    no_infected = sum(is_infected);
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
        end
        %____________________________________________________________
        R = sum(sum(rate,1));
        h= -log(rand)/R;
        ts(cnt)=t+h;
        t = t+h;
        %         z1(cnt) = no_sus;
        %         z2(cnt) = no_carrier;
        %         z3(cnt) = no_infected;
        %         z4(cnt) = no_recovered;
        
        %______________________________________________________________
        % Setting the new activity states and disease states
        
        node = rnd_draw(sum(rate,1));
        if size(rate(:,node),2)~=1
            1;
        end
        w1 = rnd_draw(rate(:,node))-1;
%         if cnt == 1
%             if node == 1
%                 a1 = a1+1;
%             elseif ismember(node, Nei_s(I1(1):I2(1)))
%                 a2 = a2+1;
%             else
%                 a3 = a3+1;
%             end
%         end
        
        is_active_new = w1*is_inactive(node)+(1-w1)*is_active(node);
        is_inactive_new = w1*is_active(node)+(1-w1)*is_inactive(node);
        is_active(node) = is_active_new;
        is_inactive(node) = is_inactive_new;
        
        ss = [eta eta_prime-eta];
%         if size(ss,1)~=1
%             1;
%         end
        w2 = rnd_draw(ss)-1;
        
        
        is_carrier_new = is_sus(node);
        is_infected_new = (1-w2)*is_carrier(node);
        is_recovered_new = w2*is_carrier(node)+is_infected(node);
        is_sus_new  = 0;
        is_carrier_pre = is_carrier(node);
        
        is_sus(node) = (1-w1)*is_sus_new+w1*is_sus(node);
        is_carrier(node) = (1-w1)*is_carrier_new+w1*is_carrier(node);
        is_infected(node) = (1-w1)*is_infected_new+w1*is_infected(node);
        is_recovered(node) = (1-w1)*is_recovered_new+w1*is_recovered(node);
        [is_sus(node) is_carrier(node) is_infected(node) is_recovered(node)];
%         if w1 == 0 && is_recovered(node) == 1
%             u = Nei_d(J1(node):J2(node));
%             ind = find(link_active(J1(node):J2(node))==1);
%             temp_neigh = u(ind);
%             per_neigh = Nei_s(I1(node):I2(node));
%         end
        %________________________________________________________________
        % Setting the new neighbourhood vector according to activity states
        
        
        set_nodes = [];
        temp_neigh = [];
        %_________________________________________________________________
        % Acitivation
        if w1*is_active(node) == 1
            if sum(link_active(J1(node):J2(node)))~=0
                1;
            end
            u = Nei_d(J1(node):J2(node));
            link_active(J1(node):J2(node)) = (rand(1,length(u))<= is_active(u)*p_link_const(node));
            ind = find(link_active(J1(node):J2(node))==1);
            temp_neigh = u(ind);
            for i = 1:length(ind)
                nd = temp_neigh(i);
                f = J1(nd):J2(nd);
                m = find(Nei_d(f)==node);
                link_active(f(m)) = 1;
            end
            
            per_neigh = Nei_s(I1(node):I2(node));
            
            if is_sus(node) == 1
%                 if rate(1,node) ~= sum(is_carrier(per_neigh))*beta_c+sum(is_infected(per_neigh))*beta_I
%                     1;
%                 end
                rate(2,node) = inactivation_rate_asym(node);
                rate(1,node) = rate(1,node)+sum(is_carrier(temp_neigh))*beta_c+sum(is_infected(temp_neigh))*beta_I;
            elseif is_carrier(node) == 1
                rate(2,node) = inactivation_rate_asym(node);
                rate(1,temp_neigh) = rate(1,temp_neigh)+is_sus(temp_neigh)*beta_c;
            elseif is_infected(node) == 1
                rate(2,node) = inactivation_rate_infected(node);
                rate(1,temp_neigh) = rate(1,temp_neigh)+is_sus(temp_neigh)*beta_I;
            end
%             if ~isempty(find(rate<0))
%                 1;
%             end
            
            %_________________________________________________________________
            % Inacitivation
            
        elseif w1*is_inactive(node) == 1
            u = Nei_d(J1(node):J2(node));
            ind = find(link_active(J1(node):J2(node))==1);
            temp_neigh = u(ind);
            for i = 1:length(ind)
                nd = temp_neigh(i);
                f = J1(nd):J2(nd);
                m = find(Nei_d(f)==node);
                link_active(f(m)) = 0;
            end
            link_active(J1(node):J2(node)) = 0;
            
            if is_sus(node) == 1
%                 if rate(1,node) < sum(is_carrier(temp_neigh))*beta_c-sum(is_infected(temp_neigh))*beta_I
%                     1;
%                 end
                rate(2,node) = activation_rate_asym(node);
                rate(1,node) = rate(1,node)-sum(is_carrier(temp_neigh))*beta_c-sum(is_infected(temp_neigh))*beta_I;
            elseif is_carrier(node) == 1
                rate(2,node) = activation_rate_asym(node);
                rate(1,temp_neigh) = rate(1,temp_neigh)-is_sus(temp_neigh)*beta_c;
            elseif is_infected(node) == 1
                rate(2,node) = activation_rate_infected(node);
                rate(1,temp_neigh) = rate(1,temp_neigh)-is_sus(temp_neigh)*beta_I;
            end
            if ~isempty(find(rate<0))
                1;
            end
            %_________________________________________________________________
            % Change in disease state
        else
            u = Nei_d(J1(node):J2(node));
            ind = find(link_active(J1(node):J2(node))==1);
            temp_neigh = u(ind);
            per_neigh = Nei_s(I1(node):I2(node));
            if is_carrier(node) == 1
                rate(1,temp_neigh) = rate(1,temp_neigh)+is_sus(temp_neigh')*beta_c;
                rate(1,per_neigh) = rate(1,per_neigh)+is_sus(per_neigh')*beta_c;
                rate(1,node) = eta_prime;
            elseif is_infected(node) == 1
                rate(1,temp_neigh) = rate(1,temp_neigh)+is_sus(temp_neigh')*(beta_I-beta_c);
                rate(1,per_neigh) = rate(1,per_neigh)+is_sus(per_neigh')*(beta_I-beta_c);
                rate(1,node) = delta;
            elseif is_recovered(node) == 1 && is_carrier_pre == 0
                rate(1,temp_neigh) = rate(1,temp_neigh)-is_sus(temp_neigh')*beta_I;
                rate(1,per_neigh) = rate(1,per_neigh)-is_sus(per_neigh')*beta_I;
                rate(1,node) = 0;
            elseif is_recovered(node) == 1 && is_carrier_pre == 1
                rate(1,temp_neigh) = rate(1,temp_neigh)-is_sus(temp_neigh')*beta_c;
                rate(1,per_neigh) = rate(1,per_neigh)-is_sus(per_neigh')*beta_c;
                rate(1,node) = 0;
            else
                1;
            end
        end
%         for h = 1:N
%             u = Nei_d(J1(h):J2(h));
%             ind1 = find(link_active_init(J1(h):J2(h))==1);
%             temp_neigh1 = u(ind1);
%             per_neigh1 = Nei_s(I1(h):I2(h));
%             neigh1 = [per_neigh1;temp_neigh1];
%             if is_sus(h)==1
%                 
%                 if rate(1,h) ~= sum(is_carrier(neigh1))*beta_c+sum(is_infected(neigh1))*beta_I
%                     h;
%                     1;
%                 end
%             elseif is_carrier(h)==1
%                 if rate(1,h) ~= eta_prime
%                     h;
%                     1;
%                 end
%             end
%         end
        
        
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
        if no_active+no_inactive ~=N
            1;
        end
        %         if no_carrier+no_infected == 0
        %             break
        %         end
        %
    end
    %     for i = 1:t_final
    %         e = find(ts<=i);
    %         if isempty(e)
    %             z1_mean(i) = ((count-1)*z1_mean(i)+N-1)/count;
    %             z2_mean(i) = ((count-1)*z2_mean(i)+1)/count;
    %             z3_mean(i) = ((count-1)*z3_mean(i))/count;
    %             z4_mean(i) = ((count-1)*z4_mean(i))/count;
    %         else
    %             k = e(end);
    %             if count >1
    %                 z1_mean(i) = ((count-1)*z1_mean(i)+z1(k))/count;
    %                 z2_mean(i) = ((count-1)*z2_mean(i)+z2(k))/count;
    %                 z3_mean(i) = ((count-1)*z3_mean(i)+z3(k))/count;
    %                 z4_mean(i) = ((count-1)*z4_mean(i)+z4(k))/count;
    %             else
    %                 z1_mean(i) = z1(k);
    %                 z2_mean(i) = z2(k);
    %                 z3_mean(i) = z3(k);
    %                 z4_mean(i) = z4(k);
    %             end
    %         end
    %     end
    %     z = [z1' z2' z3' z4'];
    %     z_mean = [z1_mean' z2_mean' z3_mean' z4_mean'];
    %     [ts_mean,z1_mean,z2_mean,z3_mean,z4_mean] = time_average(count,ts_mean,ts,z_mean,z);
    avg = avg+no_recovered;
%     toc
end
prevalence = avg/no_run/N;














