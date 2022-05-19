clc; clear all;
% mpc = loadcase('case14');
% mpc = loadcase('case24_ieee_rts');
% mpc = loadcase('pglib_opf_case24_ieee_rts');
% mpc = loadcase('pglib_opf_case118_ieee');
mpc = loadcase('pglib_opf_case1354_pegase');
% mpc = loadcase('pglib_opf_case2383wp_k');

mpopt = mpoption('out.all', 0);
mpopt.opf.flow_lim = 'P'; % default "S" is the apparent 

tic

n = 10000; % number of data sets 

N = size(mpc.bus,1); % number of buses
N_gen = size(mpc.gen,1); % number of generator
L = size(mpc.branch,1); % number of lines

index = []; % index of nonzero load
for i = 1 : N
    if (mpc.bus(i,3) ~= 0)
        index = [index;i];
    end
end
n_load = size(index,1); % number of load

load_file0_real = zeros(N,n); % load p
load_file0_reac = zeros(N,n); % load q
gen_cost0 = zeros(N_gen,n);
gen_cost0_quad = zeros(N_gen,n);
gen_pmax0 = zeros(N_gen,n);
gen_qmax0 = zeros(N_gen,n);

%% Set the cost function to be quadtatic (118 case, 24 is already quadratic)
rng(5);
quad_scale = rand(N_gen,1);

mpc.gencost(:,5) = 0.05 .* quad_scale .* mpc.gencost(:,6); % Quad term



%% Set the line limit manually
% x = mpc.branch(:,4); % reactance
% p_max = 2*pi ./ x; % approximate the max line flow
% mpc.branch(:,6) = p_max; % line flow limits on long term real power

% 24-bus system has line limit
mpc.branch(:,6) = mpc.branch(:,6) .* 1; % reduce the flow limit
p_max = mpc.branch(:,6);

% %% adjust the load for 14-bus system
% % load 5 and 3 too large
% % mpc.bus(3,3) = mpc.bus(3,3).* .5;
% % mpc.bus(5,3) = mpc.bus(5,3).* .1;
% mpc.bus(:,3) = mpc.bus(:,3).* 1;

mpc2 = mpc; % store the reference system file

%% Generate random perturbation and then run optimal power flow
% count number of effective data sets
index = [];
idx_count = 0;
index0 = zeros(n,1);
% run
for k = 1:n  % data sets

    disp('case:');
    disp(k);
    r0 = unifrnd(-0.08,0.1,N,1) .*2; % load p
    r1 = unifrnd(-0.08,0.1,N,1) .*1.5; % load q
    r2 = unifrnd(-0.1,0.1,N_gen,1) .*2; % generation cost lin
    r2_1 = unifrnd(-0.1,0.1,N_gen,1) .*2; % generation cost quad
    r3 = unifrnd(-0.01,0.2,N_gen,1) .*2; % generation limit p
    r3_1 = unifrnd(-0.01,0.2,N_gen,1) .*2; % generation limit q
    mpc.bus(:,3) = mpc.bus(:,3).*(1 + r0); % load p
    mpc.bus(:,4) = mpc.bus(:,4).*(1 + 0 + r1 .* 0.2); % load q
    mpc.gencost(:,6) = mpc.gencost(:,6).*(1 + r2); % gen cost for lin
    mpc.gencost(:,5) = mpc.gencost(:,5).*(1 + r2_1); % gen cost for quad
    mpc.gen(:,9) = mpc.gen(:,9).*(1 + r3); % gen P_max
    mpc.gen(:,4) = mpc.gen(:,4).*(1 + r3_1); % gen Q_max
    mpc.gen(:,5) = mpc.gen(:,5).*(1 - r3_1); % gen Q_min = -Q_max
    load_file0_real(:,k) = mpc.bus(:,3);
    load_file0_reac(:,k) = mpc.bus(:,4);
    gen_cost0(:,k) = mpc.gencost(:,6); % lin
    gen_cost0_quad(:,k) = mpc.gencost(:,5); % quad
    gen_pmax0(:,k) = mpc.gen(:,9);
    gen_qmax0(:,k) = mpc.gen(:,4);
    
%     result{k} = runopf(mpc);
    if min(load_file0_real(:,k)) < -500
        result{k}.success = 0;
    else
%         result{k} = rundcopf(mpc,mpopt); % DC-OPF
        result{k} = runopf(mpc,mpopt); % AC_OPF
        
        % count number of effective data sets
        if result{1,k}.success == 1
            index(end+1) = k;
            index0(k) = 1;
            idx_count = idx_count + 1;
        else
            result{k} = 0;
        end
    end
    mpc = mpc2; % reset the mpc after it is modified
end




% initializaion of voltage, bus angle and line flow
v = zeros(N,idx_count);
theta = zeros(N,idx_count);
f = zeros(L,idx_count);
load_file_real = zeros(N,idx_count);
load_file_reac = zeros(N,idx_count);
gen_cost = zeros(N_gen,idx_count);
gen_cost_quad = zeros(N_gen,idx_count);
gen_pmax = zeros(N_gen,idx_count);
gen_qmax = zeros(N_gen,idx_count);


% initialization of dual variables
lambda_opt = zeros(N,idx_count);
mu_opt_l = zeros(L,idx_count);
mu_opt_u = zeros(L,idx_count);

for i = 1:idx_count
%     v{i} = result{1,i}.bus(:,8); % voltage magnitude
%     theta{i} = result{1,i}.bus(:,9); % phase angles
%     f{i} = result{1,i}.branch(:,14); % line flow, transmission limit is unknown in the original file
    v(:,i) = result{1,index(i)}.bus(:,8);
    theta(:,i) = result{1,index(i)}.bus(:,9);
    f(:,i) = result{1,index(i)}.branch(:,14);
    load_file_real(:,i) = load_file0_real(:,index(i));
    load_file_reac(:,i) = load_file0_reac(:,index(i));
    gen_cost(:,i) = gen_cost0(:,index(i));
    gen_cost_quad(:,i) = gen_cost0_quad(:,index(i));
    gen_pmax(:,i) = gen_pmax0(:,index(i));
    gen_qmax(:,i) = gen_pmax0(:,index(i));
    % dual variables for LMP
    lambda_opt(:,i) = result{1,index(i)}.bus(:,14);
    mu_opt_l(:,i) = result{1,index(i)}.branch(:,19);
    mu_opt_u(:,i) = result{1,index(i)}.branch(:,18);
end

% line constraints
line_index_lo = zeros(L,idx_count);
line_index_up = zeros(L,idx_count);
for i = 1 : idx_count
    for j = 1 : L
        if p_max(j) - f(j,i) < 1e-3 && f(j,i) > 0
            line_index_up(j,i) = 1;
        end
        if p_max(j) + f(j,i) < 1e-3 && f(j,i) < 0
            line_index_lo(j,i) = 1;
        end
    end
end

%% bus constraints
n_gen = size(mpc2.gen,1);
n_bus = size(mpc2.bus,1);
gen_idx0 = mpc2.gen(:,1);
gen_idx = zeros(length(gen_idx0),1);
for i = 1 : length(gen_idx0)
    gen_idx(i) = find(mpc.bus(:,1) == gen_idx0(i));
end
gen_lim = mpc2.gen(:,9:10);

% generator limit
g_min = zeros(n_bus,1);
g_max = zeros(n_bus,1);
for i = 1 : n_gen
    g_min(gen_idx(i)) = g_min(gen_idx(i)) + gen_lim(i,2);
    g_max(gen_idx(i)) = g_max(gen_idx(i)) + gen_lim(i,1);
end

% get data from results
gen_data0 = zeros(n_gen,idx_count);
gen_data = zeros(n_bus,idx_count);
for i = 1:idx_count
    gen_data0(:,i) = result{1,index(i)}.gen(:,2);
end
for i = 1 : n_gen
    gen_data(gen_idx(i),:) = gen_data(gen_idx(i),:) + gen_data0(i,:);
end
% gen constraints
gen_index_lo = zeros(n_bus,idx_count);
gen_index_up = zeros(n_bus,idx_count);
for i = 1 : idx_count
    for j = 1 : n_bus
        if g_max(j) - gen_data(j,i) < 1e-3 && g_max(j) > 0
            gen_index_up(j,i) = 1;
        end
        if gen_data(j,i) -  g_min(j)< 1e-3 && g_min(j) > 0
            gen_index_lo(j,i) = 1;
        end
        % remove the idle bus w/out generators
    end
end

fprintf('Effective samples %d. \n',idx_count);
line=line_index_lo+line_index_up;
lsum=sum(line,1);
disp('Max active line:');
disp(max(lsum));
disp('Active ratio:');
disp(sum(sum(line))/169/idx_count);

toc