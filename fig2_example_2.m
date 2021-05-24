clear clc
% Part 1:subsystem sequence
% This part is to generate the subsystem switching sequence
deta_t = 2;  %the length of each interval 
L = 50;  % the total number of continuous intervals 
sys_mode = randsample([1 2],L,true,[0.2 0.8]);
%sys_mode=randsample([1 2],L,true,[0.6 0.4]);
% To generate a subsystem index sequence with given probability

t = [];
mode_index = []; %subsystem index for each time
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
xlabel('t')
ylabel('\sigma(t)')


%Part 2: Éú³ÉÂö³åÐòÁÐ
lamda = [0.2 0.3];  % % parameters for Possion renewal process
uniform_end = [1.02 1.04]; % the endpoints for uniform distribution
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
axis([0 30,1,1.10])
legend([p1,p2],'impluses under subsystem mode 1','impluses under subsystem mode 2')
xlabel('t')
ylabel('impulsive strength')


% Part3: system evolution
x_0 = [0.1;0.1;0.1];
t = 0;
X_temp_input = x_0;
X_sys_input = x_0;
X_temp_noput = x_0;
X_sys_noput = x_0;

N = 3;
for i = 1:L
    randn('state',100)
    T = deta_t; N_t = deta_t*2^8; dt = T/N_t;
    dW = sqrt(dt)*randn(N,N_t); % Brownian increments
    R = 4; Dt = R*dt; L_t = N_t/R; %  EM steps of size Dt = R*dt
    for j = 1:L_t
        t = t+Dt;
        delay_t = 1+sin(t);  %time-varying delays 
        delay_t = floor(delay_t/Dt);
        if j-delay_t <= 0
            delay_t =0;
        end
        Winc = sum((dW(:,R*(j-1)+1:R*j))');   
          
        f_fun = [];   
        g_fun = [];
        k = L_t*(i-1)+j;
        A_1 =[0.5 0.4 0.3;-0.4 0.3 0.4;0.2 -0.3 -0.2];  
        A_2 =[-1.4 0 0.3;0.5 -1.2 0.3;0.2 0 -2.8];
        f_1 = [tanh(X_sys_input(1,k-delay_t));tanh(X_sys_input(2,k-delay_t));tanh(X_sys_input(3,k-delay_t))].^2; %subsystem activation function 
        f_2 = [tanh(X_sys_input(1,k-delay_t));tanh(X_sys_input(2,k-delay_t));tanh(X_sys_input(3,k-delay_t))].^2; %subsystem activation function 
        g_1 = [0.1*X_sys_input(1,k-delay_t);0.2*X_sys_input(2,k-delay_t);0];   % noise strengthen function 
        g_2 = [0.2*X_sys_input(1,k-delay_t);0;0.1*X_sys_input(3,k-delay_t)];    % noise strengthen function 

        u_c = [(cos(k)+1);(cos(k));(cos(k)-1)]; %input function
        
        if sys_mode(i) == 1
            X_temp_input = X_temp_input + Dt*(A_1*tanh(X_sys_input(:,k))+0.2*A_1*f_1+u_c)+g_1.*Winc(1:N)';
        else
            X_temp_input = X_temp_input + Dt*(A_2*tanh(X_sys_input(:,k))+0.2*A_2*f_2+u_c)+g_2.*Winc(1:N)';
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
plot(t_plot,X_sys_input(1,:),t_plot,X_sys_input(2,:),t_plot,X_sys_input(3,:));
legend('x_1(t)','x_2(t)','x_3(t)')
xlabel('t','FontSize',12)
ylabel('x(t)','FontSize',12)

