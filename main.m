close all;
clear ;
%% x_dot = Ax, 
A = [-1   0    0   0   0   0;
     -1   1    1   0   0   0;
     1   -2   -1   -1  1   1;
     0    0    0   -1  0   0;
     -8   1   -1   -1  -2  0;
     4   -0.5 0.5   0  0   -4];

%% y = Hx
H1 = [  1 0 0 2 0 0;
        2 0 0 1 0 0];
    
H2 = [2 0 1 0 0 1];

H3 = [0 0 0 2 0 0];

H4 = [1 0 2 0 0 0;
      2 0 4 0 0 0];

H = [H1 ; H2 ; H3 ; H4];

%% Laplacian of this graph
L =[ 2  -1   0   -1
     0   1   -1   0
    -1  -1   2    0
    -1   0   0    1 ];
%% Check Observability

% Paper said (Hi , A) is observable
% obsv(H1,A)


% Paper said (H , A) is observable pair 

% Determine system observability
if rank(obsv(H,A)) == size(A, 1)
    disp('The system is observable.');
else
    disp('The system is not observable.');
end

%% The local observer gain matrices are computed as

gamma = 219.7065;

L1 = [0.1661 0.1661;
        0       0;
        0       0;
        0.4982 -0.4982;
        0       0;
        0       0 ];

L2 = [3.5770;
    -14.3468;
    0.5188;
    -5.4729;
    -14.8217;
    -5.9758];


L3 =   [0;
        0;
        0;
        0;
        0;
        0];

L4 =    [-2.2551    0.3970;
        -8.9847     -2.6908;
        1.6216      0.1859;
        6.2723      0.9189;
        -3.0056     -1.7857;
        15.9378     -1.8280];
    
    
M1 = [0.0037    0       0   0       0   0;
        0       1       0   0       0   0;
        0       0       1   0       0   0;
        0       0       0   0.0037  0   0;
        0       0       0   0       1   0;
        0       0       0   0       0   1];

M2 =    [0.0512     0.0879      -0.1315     -0.1416     -0.0429     0.0297
        0.0879      0.1796      -0.2330     -0.2519     -0.0576     0.0528
        -0.1315     -0.2330     0.3450      0.3644      0.1031      -0.0822
        -0.1416     -0.2519     0.3644      0.3981      0.1173      -0.0816
        -0.0429     -0.0576     0.1031      0.1173      0.0549      -0.0216
        0.0297      0.0528      -0.0822     -0.0816     -0.0216     0.0259];
    
M3 =    [1  0   0   0       0   0;
        0   1   0   0       0   0;
        0   0   1   0       0   0;
        0   0   0   0.0130  0   0;
        0   0   0   0       1   0;
        0   0   0   0       0   1];
    
M4 =    [0.0390     0.0293      -0.0195     -0.1186     -0.1410     0.0409;
        0.0293      0.0315      -0.0161     -0.0932     -0.1002     0.0303;
        -0.0195     -0.0161     0.0103      0.0598      0.0696      -0.0210;
        -0.1186     -0.0932     0.0598      0.3643      0.4299      -0.1273;
        -0.1410     -0.1002     0.0696      0.4299      0.5260      -0.1603;
        0.0409      0.0303      -0.0210     -0.1273     -0.1603     0.0621];
%%

x_old = [1, 3, -2, -3, -1, 2]';
x(:,1) = x_old ; 

%each local observer the initial state is taken to be zero
x_hat_1_old = [0, 0, 0, 0, 0, 0]';
x_hat_1(:,1) = x_hat_1_old ;

x_hat_2_old = x_hat_1_old ; 
x_hat_2(:,1) = x_hat_2_old;

x_hat_3_old = x_hat_1_old ;
x_hat_3(:,1) = x_hat_3_old;

x_hat_4_old = x_hat_1_old ;
x_hat_4(:,1) = x_hat_4_old;




% After Fig.2 in the paper , we have this matrices

L = [2  -1  0   -1;
     0   1  -1   0;
    -1  -1  2    0;
    -1  0   0    1];
    
% D = [];
aij = [0 1 0 1;
       0 0 1 0;
       1 1 0 0;
       1 0 0 0];
   
D = L+ aij;
   
%% Time 
%%
dt = 0.01; % 
time_sim = 10;   % time de simulation 

tot_time(:,1)=0;
t0 = 0;


for time_index = 1:time_sim/dt
    %%-----  System

    [x_dot,y1,y2,y3,y4] = system_fcn(A,H1,H2,H3,H4,x_old); 

    [t,x_new] = ode45(@(t,x)ode_system_fcn(x_dot),[t0 t0+dt],x_old);
    len_t = length(t);
    x_old = x_new(len_t,:)';
    
    %%----- Observer 

    %Formulation
    % x_hat = Ax_hat + L(y - H*x_hat) + gamma*M*Sum(aij (ˆxj - xˆi))
    
    x_hat_ij = [x_hat_1_old , x_hat_2_old , x_hat_3_old , x_hat_4_old ]; 

    x_hat_1_dot = Observer_fcn(1 , A, gamma, aij ,H1 ,M1 ,L1 ,x_hat_1(:,time_index) ,x_hat_ij , y1 );
    [t_1,x_hat_1_new ] = ode45(@(t,x)ode_system_fcn(x_hat_1_dot),[t0 t0+dt],x_hat_1_old);
    len_t = length(t_1);
    x_hat_1_old = x_hat_1_new(len_t,:)';

    x_hat_2_dot = Observer_fcn(2 , A, gamma, aij ,H2 ,M2 ,L2 ,x_hat_2(:,time_index) ,x_hat_ij , y2 );
    [t_2, x_hat_2_new] = ode45(@(t,x)ode_system_fcn(x_hat_2_dot),[t0 t0+dt],x_hat_2_old);
    len_t = length(t_2);
    x_hat_2_old = x_hat_2_new(len_t,:)';

    x_hat_3_dot = Observer_fcn(3 , A, gamma, aij ,H3 ,M3 ,L3 ,x_hat_3(:,time_index) ,x_hat_ij , y3 );
    [t_3,x_hat_3_new] = ode45(@(t,x)ode_system_fcn(x_hat_3_dot),[t0 t0+dt],x_hat_3_old);
    len_t = length(t_3);
    x_hat_3_old = x_hat_3_new(len_t,:)';

    x_hat_4_dot = Observer_fcn(4 , A, gamma, aij ,H4 ,M4 ,L4 ,x_hat_4(:,time_index) ,x_hat_ij , y4 );
    [t_4,x_hat_4_new] = ode45(@(t,x)ode_system_fcn(x_hat_4_dot),[t0 t0+dt],x_hat_4_old);
    len_t = length(t_4);
    x_hat_4_old = x_hat_4_new(len_t,:)';

    %% Save for plot
    % System State
    
    x(:,time_index+1) = x_old ;
    
    % Observer State
    x_hat_1(:,time_index+1) = x_hat_1_old;
    
    x_hat_2(:,time_index+1) = x_hat_2_old;
    
    x_hat_3(:,time_index+1) = x_hat_3_old;
    
    x_hat_4(:,time_index+1) = x_hat_4_old;
    
    t0 = time_index*dt;
    tot_time(:,time_index+1) = t0;
end

%% Plot 
figure(1)

plot(tot_time , x_hat_1(1,:) - x(1,:) )
hold on
plot(tot_time , x_hat_1(2,:) - x(2,:) )
hold on
plot(tot_time , x_hat_1(3,:) - x(3,:) )
hold on
plot(tot_time , x_hat_1(4,:) - x(4,:) )


