clear clc
% Part 1:subsystem sequence
% This part is to generate the subsystem switching sequence
deta_t = 6; %the length of each interval 
L = 10;   % the total number of continuous intervals 
sys_mode = []; 
t = [];
mode_index = []; %subsystem index for each time
for i = 1:L
     t0 = (i-1)*deta_t:0.1:i*deta_t; 
     L0 = length(t0);
     sys_mode(i) = mod(i-1,3)+1;
     mode_0 = sys_mode(i)*ones(1,L0);
     t = [t t0];
     mode_index = [mode_index mode_0];
end

figure(1)
plot(t,mode_index,'.',t,mode_index)
axis([0 deta_t*L,0,4])
xlabel('t')
ylabel('\sigma(t)')


% part2: system evolution without impulses
% This part is to show the dynamical evolution of the considered stochastic
% system without considering any impulse effect
x_0 = 0.0; % initial condition
t = 0;
X_temp_input = x_0;
X_sys_input = x_0; % system states 

N = 1;
randn('state',100)
T = deta_t; N_t = deta_t*2^8; dt = T/N_t;
R = 4; Dt = R*dt; L_t = N_t/R; %  EM steps of size Dt = R*dt
for i = 1:L
    dW = sqrt(dt)*randn(N,N_t); % Brownian increments
    for j = 1:L_t
        t = t+Dt;
        Winc = sum((dW(:,R*(j-1)+1:R*j))');    
          
        f_fun = [];   
        g_fun = [];
        k = L_t*(i-1)+j;  
        f_1 =2*cos(k)-2;  %subsystem activation function 
        f_2 =2*sin(k)+1;   %subsystem activation function 
        f_3 =2*cos(k)-1;   %subsystem activation function 
        g = 0.1*sin(k);    % noise strengthen function  

        u_c = cos(k);   % input function 
        
        if sys_mode(i) == 1
            X_temp_input = X_temp_input + Dt*(f_1*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        elseif sys_mode(i) == 2
            X_temp_input = X_temp_input + Dt*(f_2*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        else
            X_temp_input = X_temp_input + Dt*(f_3*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        end
        
        X_sys_input = [X_sys_input X_temp_input]; 
    end
end

figure(2)
t_plot = [0:Dt:L*T];
plot(t_plot,X_sys_input(1,:),'b');
legend('System states')
xlabel('t','FontSize',12)
ylabel('x(t)','FontSize',12)



%Part 3: frequent impulsive sequence
lamda = [1.2 1.4 1.6];  % parameters for Possion renewal process
lamda = 1./lamda;
uniform_end = [1.1 1.2 1.3]; % the endpoints for uniform distribution  
Tmax = deta_t;
impulive_t = [];
impulsive_str = [];

for i = 1:L
    j = 1;
    T = [];
    S = [];
    T(1) = Tmax*(i-1)+exprnd(lamda(sys_mode(i))); %random densities
    S(1)= unifrnd(1,uniform_end(sys_mode(i))); % random intensities 

    while(T(j)<Tmax*i)
        T(j+1) = T(j)+exprnd(lamda(sys_mode(i)));
        S(j+1) = unifrnd(1,uniform_end(sys_mode(i)));
        j=j+1;
    end     
    T_length = length(T);
    if T(T_length) > Tmax*i
        impulive_t = [impulive_t T(1:T_length-1)];
        impulsive_str = [impulsive_str S(1:T_length-1)];
    else
        impulive_t = [impulive_t T];
        impulsive_str = [impulsive_str S];
    end       
    
end

num_impulsive = length(impulive_t);

figure(3)
% show the random-produced mode-dependent impulse effects 
p1 = [];
p2 = [];
p3 = [];
for i = 1:num_impulsive
    if sys_mode(floor(impulive_t(i)/Tmax)+1)==1
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','r', 'LineWidth', 2);
        if isempty(p1)
            p1 = p;
        end
        hold on
    elseif sys_mode(floor(impulive_t(i)/Tmax)+1)==2
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','b', 'LineWidth', 2);
        if isempty(p2)
            p2 = p;
        end
        hold on
    else 
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','k', 'LineWidth', 2);
        if isempty(p3)
            p3 = p;
        end
        hold on
    end
end

axis([0 deta_t*L,1,1.5])
legend([p1,p2,p3],'impluses under subsystem mode 1','impluses under subsystem mode 2','impluses under subsystem mode 3')
xlabel('t')
ylabel('impulsive strength')


% Part4: system evolution with frequent impulsive sequence
% This part is similar to part 2, but impulse effects are taken into
% account
t = 0;
X_temp_input = x_0;
X_sys_input = x_0;

N = 1;
randn('state',100)
T = deta_t; N_t = deta_t*2^8; dt = T/N_t;
R = 4; Dt = R*dt; L_t = N_t/R; % L EM steps of size Dt = R*dt
for i = 1:L
    dW = sqrt(dt)*randn(N,N_t); % Brownian increments
    for j = 1:L_t
        t = t+Dt;
        Winc = sum((dW(:,R*(j-1)+1:R*j))');    
          
        f_fun = [];   
        g_fun = [];
        k = L_t*(i-1)+j;
        f_1 =2*cos(k)-2;
        f_2 =2*sin(k)+1;
        f_3 =2*cos(k)-1;
        g = 0.1*sin(k);

        u_c = cos(k);
        
        if sys_mode(i) == 1
            X_temp_input = X_temp_input + Dt*(f_1*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        elseif sys_mode(i) == 2
            X_temp_input = X_temp_input + Dt*(f_2*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        else
            X_temp_input = X_temp_input + Dt*(f_3*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        end
        
        index = find(abs((i-1)*T+Dt*j - impulive_t) <= Dt);
        if isempty(index) == 0 
            X_temp_input = prod(impulsive_str(index))*X_temp_input;
        end 
                  
        X_sys_input = [X_sys_input X_temp_input]; 
    end
end


figure(4)
t_plot = [0:Dt:L*T];
ISS_bound = [];
ISS_0 = [];
gamma = 0.16;
d_0 = 12;

plot(t_plot,X_sys_input(1,:),'b');
legend('System states')
xlabel('t','FontSize',12)
ylabel('x(t)','FontSize',12)

%Part 5: infrequent impulsive sequence
lamda = [0.20 0.25 0.30];
lamda = 1./lamda;
uniform_end = [1.02 1.04 1.06];
Tmax = deta_t;
impulive_t = [];
impulsive_str = [];

for i = 1:L
    j = 1;
    T = [];
    S = [];
    T(1) = Tmax*(i-1)+exprnd(lamda(sys_mode(i)));
    S(1)= unifrnd(1,uniform_end(sys_mode(i)));

    while(T(j)<Tmax*i)
        T(j+1) = T(j)+exprnd(lamda(sys_mode(i)));
        S(j+1) = unifrnd(1,uniform_end(sys_mode(i)));
        j=j+1;
    end     
    T_length = length(T);
    if T(T_length) > Tmax*i
        impulive_t = [impulive_t T(1:T_length-1)];
        impulsive_str = [impulsive_str S(1:T_length-1)];
    else
        impulive_t = [impulive_t T];
        impulsive_str = [impulsive_str S];
    end       
    
end

num_impulsive = length(impulive_t);
figure(5)
p1 = [];
p2 = [];
p3 = [];

for i = 1:num_impulsive
    if sys_mode(floor(impulive_t(i)/Tmax)+1)==1
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','r', 'LineWidth', 2);
        if isempty(p1)
            p1 = p;
        end
        hold on
    elseif sys_mode(floor(impulive_t(i)/Tmax)+1)==2
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','b', 'LineWidth', 2);
        if isempty(p2)
            p2 = p;
        end
        hold on
    else 
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','k', 'LineWidth', 2);
        if isempty(p3)
            p3 = p;
        end
        hold on
    end
end
axis([0 deta_t*L,1,1.07])
legend([p1,p2,p3],'impluses under subsystem mode 1','impluses under subsystem mode 2','impluses under subsystem mode 3')
xlabel('t')
ylabel('impulsive strength')


% Part6: system evolution with infrequent impulsive sequence
%This part is similar to part 4
t = 0;
X_temp_input = x_0;
X_sys_input = x_0;

N = 1;
T = deta_t; N_t = deta_t*2^8; dt = T/N_t;
R = 4; Dt = R*dt; L_t = N_t/R; % L EM steps of size Dt = R*dt
for i = 1:L
    dW = sqrt(dt)*randn(N,N_t); % Brownian increments
    for j = 1:L_t
        t = t+Dt;
        Winc = sum((dW(:,R*(j-1)+1:R*j))');    
          
        f_fun = [];   
        g_fun = [];
        k = L_t*(i-1)+j;
        f_1 =2*cos(k)-2;
        f_2 =2*sin(k)+1;
        f_3 =2*cos(k)-1;
        g = 0.1*sin(k);

        u_c = cos(k);
        
        if sys_mode(i) == 1
            X_temp_input = X_temp_input + Dt*(f_1*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        elseif sys_mode(i) == 2
            X_temp_input = X_temp_input + Dt*(f_2*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        else
            X_temp_input = X_temp_input + Dt*(f_3*X_sys_input(:,j)+u_c)+g*Winc(1:N)';
        end
        
        index = find(abs((i-1)*T+Dt*j - impulive_t) <= Dt);
        if isempty(index) == 0 
            X_temp_input = prod(impulsive_str(index))*X_temp_input;
        end 
                  
        X_sys_input = [X_sys_input X_temp_input]; 
    end
end


figure(6)
t_plot = [0:Dt:L*T];
ISS_bound = [];
ISS_0 = [];

for i = 1:length(t_plot)
    ISS_0 = [ISS_0 cos(t_plot(i))^2];
    ISS_bound = [ISS_bound max(ISS_0)];
end

plot(t_plot,X_sys_input(1,:),'b');
legend('System states','Reference ISS bound')
axis([0 60,-1,1.20])
xlabel('t','FontSize',12)
ylabel('x(t)','FontSize',12)
