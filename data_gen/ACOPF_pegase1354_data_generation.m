clc; clear all;
% mpc = loadcase('case14');
% mpc = loadcase('case24_ieee_rts');
mpc = loadcase('pglib_opf_case1354_pegase');
% mpc = loadcase('pglib_opf_case118_ieee');

mpopt = mpoption('out.all', 0);
mpopt.opf.flow_lim = 'S'; % or P

n = 10; % number of data sets 

N = size(mpc.bus,1); % number of buses
L = size(mpc.branch,1); % number of lines

index = []; % index of nonzero load
for i = 1 : N
    if (mpc.bus(i,3) ~= 0)
        index = [index;i];
    end
end
n_load = size(index,1); % number of load

load_file0_real = zeros(N,n); % did not use n_load due to the storage in package
load_file0_reac = zeros(N,n);

%% Set the line limit manually for 14-bus
% x = mpc.branch(:,4); % reactance
% p_max = 2*pi ./ x; % approximate the max line flow
% mpc.branch(:,6) = p_max; % line flow limits on long term real power

% 24-bus system has line limit
mpc.branch(:,6) = mpc.branch(:,6) .* 1; % reduce the flow limit
p_max = mpc.branch(:,6);

%% adjust the load for 14-bus system
% load 5 and 3 too large
% mpc.bus(3,3) = mpc.bus(3,3).* .5;
% mpc.bus(5,3) = mpc.bus(5,3).* .1;
mpc.bus(:,3) = mpc.bus(:,3).* 1;

mpc2 = mpc; % store the reference system file

%% Generate random load and then run optimal power flow
for k = 1:n  % data sets
%     for i = 1:size(index,1)
%         r = (rand - 0.5)*0.2; % random number from -0.1 to 0.1
%         mpc.bus(i,3) = mpc.bus(i,3).*(1 + r);
%     end
    r0 = unifrnd(-0.08,0.1,N,1) .*2; % load p
    r1 = unifrnd(-0.08,0.1,N,1) .*2; % load q
    % active power
    mpc.bus(:,3) = mpc.bus(:,3).*(1 + 0 + r0);
    % reactive power
    mpc.bus(:,4) = mpc.bus(:,4).*(1 + 0 + r1 .* 0.2);
    load_file0_real(:,k) = mpc.bus(:,3);
    load_file0_reac(:,k) = mpc.bus(:,4);
    
%     result{k} = runopf(mpc);
    if min(load_file0_real(:,k)) < -300
        result{k}.success = 0;
        disp('load too samll')
    else
%         result{k} = rundcopf(mpc,mpopt);
        result{k} = runopf(mpc,mpopt);
    end
    mpc = mpc2; % reset the mpc after it is modified
end

% count number of effective data sets
index = [];
idx_count = 0;
index0 = zeros(n,1);
for i = 1 : n
    if result{1,i}.success == 1
        index(end+1) = i;
        index0(i) = 1;
        idx_count = idx_count + 1;
    end
end

% initializaion of voltage, bus angle and line flow
v = zeros(N,idx_count);
theta = zeros(N,idx_count);
f = zeros(L,idx_count);
load_file_real = zeros(N,idx_count);
load_file_reac = zeros(N,idx_count);

for i = 1:idx_count
%     v{i} = result{1,i}.bus(:,8); % voltage magnitude
%     theta{i} = result{1,i}.bus(:,9); % phase angles
%     f{i} = result{1,i}.branch(:,14); % line flow, transmission limit is unknown in the original file
    v(:,i) = result{1,index(i)}.bus(:,8);
    theta(:,i) = result{1,index(i)}.bus(:,9);
    f(:,i) = result{1,index(i)}.branch(:,14);
    load_file_real(:,i) = load_file0_real(:,index(i));
    load_file_reac(:,i) = load_file0_reac(:,index(i));
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

% % generator limit
% g_min = zeros(n_bus,1);
% g_max = zeros(n_bus,1);
% for i = 1 : n_gen
%     g_min(gen_idx(i)) = g_min(gen_idx(i)) + gen_lim(i,2);
%     g_max(gen_idx(i)) = g_max(gen_idx(i)) + gen_lim(i,1);
% end

% get data from results
gen_data0 = zeros(n_gen,idx_count);
gen_data = zeros(n_bus,idx_count);
for i = 1:idx_count
    gen_data0(:,i) = result{1,index(i)}.gen(:,2);
end
for i = 1 : n_gen
    gen_data(gen_idx(i),:) = gen_data(gen_idx(i),:) + gen_data0(i,:);
end
% % gen constraints
% gen_index_lo = zeros(n_bus,idx_count);
% gen_index_up = zeros(n_bus,idx_count);
% for i = 1 : idx_count
%     for j = 1 : n_bus
%         if g_max(j) - gen_data(j,i) < 1e-3 && g_max(j) > 0
%             gen_index_up(j,i) = 1;
%         end
%         if gen_data(j,i) -  g_min(j)< 1e-3 && g_min(j) > 0
%             gen_index_lo(j,i) = 1;
%         end
%         % remove the idle bus w/out generators
%     end
% end

fprintf('Effective samples %d. \n',idx_count);
line=line_index_lo+line_index_up;
lsum=sum(line,1);
disp('Max active line:');
disp(max(lsum));
disp('Active ratio:');
disp(sum(sum(line))/34/idx_count);