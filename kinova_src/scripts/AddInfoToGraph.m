%% Add Info to Already Made Graph

%% Load Info

config_struct = load('PlannerGraphResult2.mat');
q_valid_list = config_struct.q_valid_list;

%% Load Graph

graph_struct = load('QGraph700kPartial.mat');
Q_Graph = graph_struct.Q_Graph;

%% Add Info to Graph

Q_Graph.Nodes.Configuration = q_valid_list';

%% Save New Graph

save('QGraph700kPartial_AllInfo.mat','Q_Graph','q_valid_list')
