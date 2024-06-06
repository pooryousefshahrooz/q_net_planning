# This is the configuration file for the experiments. 

class Config():
  """ this class includes all the experiments setup"""
  edge_fidelity_ranges = [1.0]
  q_values = [0.5]# Different values of q (swap success probability) that we want to evaluate in the experiment
  path_searching_cut_off = 8
  memory_size_max = 60
  memory_size_min = 1
  network_topology = "ESnet2.gml"# Dumbbell,Random,"ESnet2.gml",Repeater_chain,SurfnetCore.gml,"Random
  left_ESnet_subgraph = False
  number_of_nodes_in_random = [1]
  number_of_events = [1]
  mu_values = [1]
    
  schemes = ["linear_link_based"]#["link_based","exhaustive","linear_link_based"]
  link_based_solver = "CPLEX"#"exhaustive"#"Bonmin","Baron
  repeating_times = 1
  M_values = [100]
  D_values= [10]
  repeater_min_decoherence = 100
  repeater_max_decoherence= 1400#micro seconds
  repeater_decoherence_incresing_step = 50
  end_node_min_decoherence = 3000
  end_node_max_decoherence= 5000#micro seconds
  end_node_decoherence_incresing_step = 50
  checking_repeater_memory_life_time_flag = False
  checking_end_node_memory_life_time_flag = False
  
  R_values = [10,4,6,8]
  R_values = [2,8,12,16,18,5,10,15,20]
  
  edge_F_values = [1.0]# we assume all links have this fidelity
  include_fidelity_in_utility = True
  relaxing_QCAST_formulation = True
  utility_type = "NGTV"
  dynamic_system_flag = False
  weighted_sum_rate_objective = False
  include_log = True
  planning_with_wrong_weights = False
  save_QCAST_equation_results = False
  toplogy_each_path_W_QCAST_equation_result = "results/toplogy_each_path_resultv_final.csv"
  weight_upper_bounds = [1]
  q_values = [0.5]
  set_of_number_of_user_pairs = [3]
  lower_bound_distance_between_pairs = 200
  upper_bound_distance_between_pairs = 250

  square_dim_in_km_values = [100]
  number_of_paths_values = [6000]
  K_values = [1]
  end_user_memory_set = [10]
  W_step_granularities = [1]
  lengths_of_middle_link = [100]
  max_dist =50
  results_file_name = "results/utility_maximization_results.csv"

def get_config(FLAGS):
  config = Config
  for k, v in FLAGS.__flags.items():
    if hasattr(config, k):
      setattr(config, k, v.value)

  return config
