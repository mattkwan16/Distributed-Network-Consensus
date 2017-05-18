%% EE130: Class Presentation Simulation

% Matt Kwan
% Due: 3/16/16

%% Setting the parameters
clear
close all

TIME_ITER = .001;
MAX_TIMESTEPS = 25;
GRAPH = [ 1 1 0 0;
          1 1 1 1;
          0 1 1 1;
          0 1 1 1; ];
NUM_NODES = 4;
      
% GRAPH = [ ...
% 	0	0	1	1	0	0	0	0	0	0;
% 	0	0	1	1	0	0	0	0	0	0;
% 	1	1	0	0	1	1	0	0	0	0;
% 	1	1	0	0	1	1	0	0	0	0;
% 	0	0	1	1	0	0	1	1	0	0;
% 	0	0	1	1	0	0	1	1	0	0;
% 	0	0	0	0	1	1	0	0	1	1;
% 	0	0	0	0	1	1	0	0	1	1;
% 	0	0	0	0	0	0	1	1	0	0;
% 	0	0	0	0	0	0	1	1	0	0; ];
% NUM_NODES = 10;
          
INIT_OFFSET = [2 3 8 1 12 1 3 3 9 10] * TIME_ITER;  
INIT_OFFSET = [2 3 8 1 ] * TIME_ITER;
rates = [1.1 .9 .75 1.3] * TIME_ITER; % values are per TIME_ITER
p = .6; % Tuning -- Set by authors; see page 2293

%% Applying the algorithm
skew_vir = zeros(MAX_TIMESTEPS, NUM_NODES);
skew_rel = cell(MAX_TIMESTEPS, 1);
time_loc = zeros(MAX_TIMESTEPS, NUM_NODES); 
time_vir = zeros(MAX_TIMESTEPS, NUM_NODES); 
offset_vir = zeros(MAX_TIMESTEPS, NUM_NODES);

skew_vir(1, :) = ones(1, NUM_NODES);
skew_rel(:) = {GRAPH};
time_loc(1, :) = INIT_OFFSET;
time_vir(1, :) = time_loc(1,:);

% Assumes that there is TX/RX between all nodes at every time step
for t = 2:1:MAX_TIMESTEPS
    time_loc(t, :) = time_loc(t-1, :) + rates;
    % Go through the graph for links
    for i=1:1:NUM_NODES
        for j=1:1:NUM_NODES
            if GRAPH(i, j) ~= 0 % link is found where i RXs from j
                % A. RELATIVE SKEW ESTIMATION
                skew_rel{t}(i,j) = p*skew_rel{t-1}(i,j) + (1-p) ...
                                   *(time_loc(t,j)-time_loc(t-1,j)) ...
                                   /(time_loc(t,i)-time_loc(t-1,i));
                % B. SKEW COMPENSATION
                skew_vir(t,i) = p*skew_vir(t-1,i) ...
                                + (1-p)*skew_rel{t-1}(i,j)*skew_vir(t-1,j);
                % C. OFFSET COMPENSATION
                offset_vir(t,i) = offset_vir(t-1,i) + (1-p) ...
                               * (skew_vir(t-1,j)*time_loc(t-1,j) ...
                               + offset_vir(t-1, j) ...
                               - skew_vir(t-1,i)*time_loc(t-1,i) ...
                               - offset_vir(t-1, i));
            end
        end
    end
    time_vir(t, :) = skew_vir(t, :).*time_loc(t, :) + offset_vir(t,:);
end

%% Plotting the data
figure;
for i=1:1:NUM_NODES
    plot(1:1:MAX_TIMESTEPS, time_loc(:, i),'color', rand(1,3)), hold on
end
title('Local Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;

figure;
for i=1:1:NUM_NODES
    plot(1:1:MAX_TIMESTEPS, time_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Time within Nodes');
xlabel('Iterations');
ylabel('Time (s)');
hold off;

figure;
temp = zeros(1,MAX_TIMESTEPS);
for i=1:1:NUM_NODES
    for j=1:1:MAX_TIMESTEPS
        temp(j) = skew_rel{j}(2,i);
    end
    plot(1:1:MAX_TIMESTEPS, temp,'color', rand(1,3)), hold on
end
% Note: a line at 0 indicates no link. A line at 1 indicates self.
% Expectation: ([Other Node Speed] / [This Node Speed]) convergence
title('Relative Skew Estimation as Seen by Node 2')
xlabel('Iterations');
ylabel('Skew (rate/rate)');
hold off;

figure;
for i=1:1:NUM_NODES
    plot(1:1:MAX_TIMESTEPS, skew_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Skew Estimation within Nodes');
xlabel('Iterations');
ylabel('Skew (rate/rate)');
hold off;

figure;
for i=1:1:NUM_NODES
    plot(1:1:MAX_TIMESTEPS, offset_vir(:, i),'color', rand(1,3)), hold on
end
title('Virtual Offset Estimation within Nodes');
xlabel('Iterations');
ylabel('Offset (sec)');
hold off;

figure;
error = zeros(MAX_TIMESTEPS, NUM_NODES);
for j=1:1:MAX_TIMESTEPS
    for i=1:1:NUM_NODES
        error(j,i) = 100*(time_vir(j,i) - mean(time_vir(j,:)))/mean(time_vir(j,:));
    end
end
for i=1:1:NUM_NODES
    plot(1:1:MAX_TIMESTEPS, error(:, i),'color', rand(1,3)), hold on
end
title('Error Percentage from Instantaneous Mean of Local Times');
xlabel('Iterations');
ylabel('Error Percentage (%)');
hold off;