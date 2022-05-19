
%% system parameters
% idx_count: number of valid data points
% N = size(mpc.bus,1); % number of buses
% L = size(mpc.branch,1); % number of lines
% n_gen = size(mpc2.gen,1); % number of generators
% real and reactive load power
% load_file_real = zeros(N,idx_count);
% load_file_reac = zeros(N,idx_count);
%

%% line power
% check double lines
M = mpc.branch(:,1:2); % bus name, need to change to bus index
n_line = size(M,1);
M_bus = zeros(n_line,2);
for j = 1 : 2
    for i = 1 : n_line
        M_bus(i,j) = find(mpc.bus(:,1) == M(i,j));
    end
end
M = M_bus;

n0 = size(M,1);

double_line_idx = [];
for i = 1 : n0 - 1
    if M(i,1) == M(i+1,1) && M(i,2) == M(i+1,2)
        double_line_idx(end+1,1:2) = [i,i+1];
    end
end
n_double_line = size(double_line_idx,1);
fprintf('There are %d pairs of double lines. \n',n_double_line);

% line flow limit(real, apparent power)
p_line_max = mpc.branch(:,6);


% apparent power
% here we use 'from', 'to' is 15,16
f = zeros(L,idx_count);
for i = 1:idx_count
    p_real = result{1,index(i)}.branch(:,14);
    p_reac = result{1,index(i)}.branch(:,15);
    f(:,i) = sqrt(p_real.^2 + p_reac.^2);
end

% % concatenation for double lines
% M0 = M;
% f0 = f;
% p_line_max0 = p_line_max;
% for i = 1 : n_double_line
%     % connection
%     M0(double_line_idx(i,1):end-1,:) = M(double_line_idx(i,2):end,:);
%     % flow
%     f0(double_line_idx(i,1),:) = f(double_line_idx(i,1),:) + f(double_line_idx(i,2),:);
%     f0(double_line_idx(i,2):end-1,:) = f(double_line_idx(i,2)+1:end,:);
%     % flow limit
%     p_line_max0(double_line_idx(i,1),:) = p_line_max(double_line_idx(i,1),:) + p_line_max(double_line_idx(i,2),:);
%     p_line_max0(double_line_idx(i,2):end-1,:) = p_line_max(double_line_idx(i,2)+1:end,:);
% end
% M = M0(1:L-n_double_line,:);
% f = f0(1:L-n_double_line,:);
% p_line_max = p_line_max0(1:L-n_double_line,:);

% concatenation for double lines
M0 = M;
f0 = f;
p_line_max0 = p_line_max;
for i = 1 : n_double_line
    % incidence
    M0(double_line_idx(i,2),:) = [0,0];
    % flow
    f0(double_line_idx(i,1),:) = f(double_line_idx(i,1),:) + f(double_line_idx(i,2),:);
    f0(double_line_idx(i,2),:) = zeros(size(f0(double_line_idx(i,2),:)));
    % flow limit
    p_line_max0(double_line_idx(i,1),:) = p_line_max(double_line_idx(i,1),:) + p_line_max(double_line_idx(i,2),:);
    p_line_max0(double_line_idx(i,2),:) = zeros(size(p_line_max0(double_line_idx(i,2),:)));
end
for i = 1 : n0
    if M0(i,1) == 0
        M0(i:end-1,:) = M0(i+1:end,:);
        f0(i:end-1,:) = f0(i+1:end,:);
        p_line_max0(i:end-1,:) = p_line_max0(i+1:end,:);
    end
end
M = M0(1:L-n_double_line,:);
f = f0(1:L-n_double_line,:);
p_line_max = p_line_max0(1:L-n_double_line,:);
n1 = size(M,1);


%% generator concatenation
% every bus has only single generator (unlike the 24bus system)
gen_idx = mpc2.gen(:,1); % bus number of generators
gen_real_lim = mpc2.gen(:,9:10);
gen_reac_lim = mpc2.gen(:,4:5);

% generator data
% AC case
gen_data0_real = zeros(n_gen,idx_count);
gen_data0_reac = zeros(n_gen,idx_count);
% already concatenated
% gen_data_real = zeros(n_bus,idx_count);
% gen_data_reac = zeros(n_bus,idx_count);
for i = 1:idx_count
    gen_data0_real(:,i) = result{1,index(i)}.gen(:,2);
    gen_data0_reac(:,i) = result{1,index(i)}.gen(:,3);
end
% for i = 1 : n_gen
%     gen_data_real(gen_idx(i),:) = gen_data(gen_idx(i),:) + gen_data0(i,:);
%     gen_data_reac(gen_idx(i),:) = gen_data(gen_idx(i),:) + gen_data0(i,:);
% end

%% generator Pmax&cost datasets
% concatenate generators to nodes
gen_cost_node = zeros(N,idx_count);
gen_cost_node_quad = zeros(N,idx_count);
% gen_pmax_node = zeros(N,idx_count);
gen_real_lim_node = zeros(N,idx_count);
gen_reac_lim_node = zeros(N,idx_count);
% power limit
gen_data0_real_node = zeros(N,idx_count);
gen_data0_reac_node = zeros(N,idx_count);

gen_idx0 = mpc2.gen(:,1);
gen_idx1 = zeros(length(gen_idx0),1);
for i = 1 : length(gen_idx0)
    gen_idx1(i) = find(mpc.bus(:,1) == gen_idx0(i));
end

for i = 1 : N_gen
    gen_idx = gen_idx1(i);
%     gen_idx = mpc.gen(i,1);
    gen_cost_node(gen_idx,:) = gen_cost_node(gen_idx,:)+ gen_cost(i,:);  
    gen_cost_node_quad(gen_idx,:) = gen_cost_node_quad(gen_idx,:)+ gen_cost_quad(i,:); 
%     gen_pmax_node(gen_idx,:) = gen_pmax_node(gen_idx,:)+ gen_pmax(i,:);
    gen_real_lim_node(gen_idx,:) = gen_real_lim_node(gen_idx,:) + gen_pmax(i,:);
    gen_reac_lim_node(gen_idx,:) = gen_reac_lim_node(gen_idx,:) + gen_qmax(i,:);
    % power limit
    gen_data0_real_node(gen_idx,:) = gen_data0_real_node(gen_idx,:) + gen_data0_real(i,:);
    gen_data0_reac_node(gen_idx,:) = gen_data0_reac_node(gen_idx,:) + gen_data0_reac(i,:);
end


%% Generate dataset
n = idx_count;
% system parameters, line-bus index, AC load, AC generator limit, 
% AC generation
data = zeros(L-n_double_line,10*n+5);

a_sys = zeros(L-n_double_line,1);
a_sys(1:4) = [N,L-n_double_line,n_load,n];
a_sys(5) = n_gen;% number of generators

% 1: index
data(:,1) = a_sys;
% 2-3: connection
data(:,2:3) = M;
% 4-(2n+3): real and reactive load
data(1:N,4:2*n+3) = [load_file_real,load_file_reac];
% (2n+4)-(2n+7): real and reactive generator power limit
% data(1:n_gen,2*n+4:2*n+7) = [gen_real_lim,gen_reac_lim];
data(1:N,2*n+4:4*n+3) = [gen_real_lim_node,gen_reac_lim_node];
% (2n+8)-(4n+7): real and reactive generator power limit
% data(1:n_gen,2*n+8:4*n+7) = [gen_data0_real,gen_data0_reac];
data(1:N,4*n+4:6*n+3) = [gen_data0_real_node,gen_data0_reac_node];
% (4n+8)-(5n+8): line power limit and data
data(:,6*n+4:7*n+4) = [p_line_max,f];
data(1:n_gen,7*n+5) = gen_idx;
data(1:N,7*n+6:8*n+5) = lambda_opt;% dataset for lmp

data(1:N,8*n+6:9*n+5) = gen_cost_node; % generation cost lin
data(1:N,9*n+6:10*n+5) = gen_cost_node_quad; % generation cost quad


%% Write the csv/txt file

filename1 = 'ieee1354pegase_ac_9974_cases_v_04022022.txt';
csvwrite(filename1,v);

filename2 = 'ieee1354pegase_ac_9974_cases_theta_04022022.txt';
csvwrite(filename2,theta);

filename = 'ieee1354pegase_ac_9974_cases_04022022.txt';

csvwrite(filename,data);














