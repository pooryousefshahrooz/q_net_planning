# This is the configuration file for the experiments. 

class Config():
  """ this class includes all the experiments setup"""
  path_searching_cut_off = 8
  memory_size_max = 60
  memory_size_min = 1
  network_topology = "ESnet2.gml"# Dumbbell,Random,"ESnet2.gml",Repeater_chain,SurfnetCore.gml,"Random
  left_ESnet_subgraph = False
  number_of_nodes_in_random = [1]# set one if you use SURFnet or ESnet
  number_of_events = [1] # for simulating dynamic phase of the network
  mu_values = [1] 
    
  schemes = ["linear_link_based"]#["link_based","exhaustive","linear_link_based"] # different ways of solving the optimization problem
  link_based_solver = "CPLEX"#"exhaustive"#"Bonmin","Baron
  repeating_times = 1 # repeatting the experiment 
  D_values= [10] # set of values for repeaters memory budget
  repeater_min_decoherence = 100 # minimum repeater coherence time in microseconds
  repeater_max_decoherence= 1400#micro seconds
  repeater_decoherence_incresing_step = 50 # increaenting the coherence time of repeaters step size (used when experimenting with different values) 
  end_node_min_decoherence = 3000 # end nodes minimum coherence time in micoseconds 
  end_node_max_decoherence= 5000# micro seconds
  end_node_decoherence_incresing_step = 50
  checking_repeater_memory_life_time_flag = False # set to True of we want to restrict paths 
  checking_end_node_memory_life_time_flag = False 
  
  R_values = [10,4,6,8] # set of values for number of repetaers budget
  
  
  edge_F_values = [1.0]# set of values to experiment for link level fidelity. we assume all links have this fidelity
  include_fidelity_in_utility = True 
  relaxing_QCAST_formulation = True
  utility_type = "NGTV" # Nagativity utility fucntion
  dynamic_system_flag = False
  weighted_sum_rate_objective = False
  include_log = True
  planning_with_wrong_weights = False
  save_QCAST_equation_results = False
  toplogy_each_path_W_QCAST_equation_result = "results/toplogy_each_path_QCAST_result_final.csv"
  weight_upper_bounds = [1]
  q_values = [0.5]# Different values of q (swap success probability) that we want to evaluate in the experiment
  set_of_number_of_user_pairs = [3] # number of user pairs in the network
  lower_bound_distance_between_pairs = 200 # in km
  upper_bound_distance_between_pairs = 250 # in km

  square_dim_in_km_values = [100] # in km
  number_of_paths_values = [6000] # number of paths between each pair of user pairs
  K_values = [1] # set of number of paths allowed to be used for each user pair. One means each user pair can use at most one path.
  end_user_memory_set = [10] # the memory budget of end nodes
  W_step_granularities = [1] # experiemnting on the width of paths for brute force search
  lengths_of_middle_link = [100] # in km for dumbell shape topology
  max_dist =50 # max distance between auxiliary nodes in ESnet topology
  results_file_name = "results/utility_maximization_results.csv"

def get_config(FLAGS):
  config = Config
  for k, v in FLAGS.__flags.items():
    if hasattr(config, k):
      setattr(config, k, v.value)

  return config
