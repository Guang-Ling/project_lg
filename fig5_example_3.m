%ISS stability of switched genetic regulatory networks
clear,clc
% Part 1:按指定概率生成子系统代码
deta_t = 4; 
L = 25;
sys_mode = randsample([1 2],L,true,[0.6 0.4]);

t = [];
mode_index = [];
for i = 1:L
    t0 = (i-1)*deta_t:0.1:i*deta_t;
    L0 = length(t0);
    mode_0 = sys_mode(i)*ones(1,L0);
    t = [t t0];
    mode_index = [mode_index mode_0];
end


figure(1)
plot(t,mode_index,'.',t,mode_index)
axis([0 deta_t*L,0,3])
xlabel('t','FontSize',14)
ylabel('\sigma(t)','FontSize',14)


%Part 2: 生成脉冲序列
lamda = [0.3 0.5];
uniform_end = [1.1 1.2];
Tmax = deta_t;
impulive_t = [];
impulsive_str = [];

for i = 1:L
    j = 1;
    T = [];
    S = [];
    T(1) = Tmax*(i-1)+exprnd(1/lamda(sys_mode(i)));
    S(1)= unifrnd(1,uniform_end(sys_mode(i)));

    while(T(j)<Tmax*i)
        T(j+1) = T(j)+exprnd(1/lamda(sys_mode(i)));
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
figure(2)
p = [];
p1 = [];
p2 = [];
for i = 1:num_impulsive
    if sys_mode(floor(impulive_t(i)/Tmax)+1)==1
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','r', 'LineWidth', 2);
        if isempty(p1)
            p1 = p;
        end
        hold on
    else
        p = plot([impulive_t(i),impulive_t(i)],[0,impulsive_str(i)],'linestyle','-', 'Color','b', 'LineWidth', 2);
        if isempty(p2)
            p2 = p;
        end
        hold on
    end
end
axis([0 30,1,1.25])
legend([p1,p2],'impluses under subsystem mode 1','impluses under subsystem mode 2')
xlabel('t','FontSize',14)
ylabel('impulsive strength','FontSize',14)

% Part3: Generation of random topologies
N_network = 1000;
P_network = 0.8;
adj_G1 = rand(N_network,N_network)*0.2-0.1;
adj_G2 = rand(N_network,N_network)*0.2-0.1;

for m=1:N_network
    for k = 1:N_network
        if(rand(1,1)<P_network)
            adj_G1(m,k)=0;
%            adj_G1(k,m)=1;
        end
        if(rand(1,1)<P_network)
            adj_G2(m,k)=0;
%            adj_G1(k,m)=1;
        end
    end
end

%Part4: system evolution
% parameters value
D1 = [-2.2,0;1,-2.3];
D2 = [-2.5,0;2.2,-2.9];
Gamma_sys = [0,1;0,0];
b = repmat({D1},N_network,1);
D_sys_1 = blkdiag(b{:});
b = repmat({D2},N_network,1);
D_sys_2 = blkdiag(b{:});

G_sys_1 = kron(adj_G1,Gamma_sys);
G_sys_2 = kron(adj_G2,Gamma_sys);
N = N_network*2;

% simulate Brownian movement
% randn('state',100)
% T = 30; N_t = 3000; dt = T/N_t;
% dW = sqrt(dt)*randn(N+1,N_t); % Brownian increments
% R = 4; Dt = R*dt; L = N_t/R; % L EM steps of size Dt = R*dt

x_0 = rand(N_network*2,1);
t = 0;
X_temp_input = x_0;
X_sys_input = x_0;
X_temp_noput = x_0;
X_sys_noput = x_0;

for i = 1:L
    randn('state',100)
    T = deta_t; N_t = deta_t*2^8; dt = T/N_t;
    dW = sqrt(dt)*randn(N,N_t); % Brownian increments
    R = 4; Dt = R*dt; L_t = N_t/R; %  EM steps of size Dt = R*dt
    for j = 1:L_t
        t = t+Dt;
        delay_t = (sin(t))^2;  %time-varying delays 
        delay_t = floor(delay_t/Dt);
        if j-delay_t <= 0
            delay_t =0;
        end
        Winc = sum((dW(:,R*(j-1)+1:R*j))');   
          
        k = L_t*(i-1)+j;
        
        f_sys = (X_sys_input(:,k-delay_t)).^2./(1+(X_sys_input(:,k-delay_t)).^2); %system activation function 
        f_sys(1:2:N) = 0;
        g_sys = 0.1*X_sys_input(:,k-delay_t); % noise strengthen function 
        
        b = [(cos(k)).^2;(sin(k)).^2];
        u_c = repmat(b,N_network,1); %input function
        
        if sys_mode(i) == 1
            X_temp_input = X_temp_input + Dt*(D_sys_1*X_sys_input(:,k)+G_sys_1*f_sys+u_c)+g_sys.*Winc(1:N)';
        else
            X_temp_input = X_temp_input + Dt*(D_sys_2*X_sys_input(:,k)+G_sys_2*f_sys+u_c)+g_sys.*Winc(1:N)';
        end
        
        index = find(abs((i-1)*T+Dt*j - impulive_t) <= Dt);
        if isempty(index) == 0 
            X_temp_input = prod(impulsive_str(index))*X_temp_input;
        end 
                        
        X_sys_input = [X_sys_input X_temp_input]; 
    end
end

figure(3)
t_plot = [0:Dt:L*T];
plot(t_plot,X_sys_input(1,:),t_plot,X_sys_input(100,:),t_plot,X_sys_input(500,:),t_plot,X_sys_input(1000,:));
legend('m_1(t)','m_{100}(t)','m_{500}(t)','m_{1000}(t)')
xlabel('t','FontSize',14)
ylabel('m(t)','FontSize',14)




