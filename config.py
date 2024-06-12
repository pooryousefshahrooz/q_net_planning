# This is the configuration file for the experiments. 

class Config():
  """ this class includes all the experiments setup """

  """ general hyperparameters """
  repeating_times = 1 # repeating the experiment 
  include_fidelity_in_utility = True 
  include_log = True
  q_values = [0.5]# Different values of q (swap success probability) that we want to evaluate in the experiment
  utility_type = "NGTV" # Nagativity utility fucntion
  number_of_events = [1] # for simulating dynamic phase of the network
  
  """ solver setup """
  schemes = ["exhaustive"]#["link_based","exhaustive","linear_link_based"] # different ways of solving the optimization problem
  link_based_solver = "CPLEX"#"exhaustive"#"Bonmin","Baron
  relaxing_QCAST_formulation = True
  W_step_granularities = [1] # experiemnting on the width of paths for brute force search
    
  """ saving files """
  save_QCAST_equation_results = False
  results_file_name = "results/utility_maximization_results.csv"


    
  """ network topology """
  toplogy_each_path_W_QCAST_equation_result = "results/toplogy_each_path_QCAST_result_final.csv"
  network_topology = "Dumbbell"# Dumbbell,Random,"ESnet2.gml",Repeater_chain,SurfnetCore.gml
  left_ESnet_subgraph = False
  max_dist =50 # max distance between auxiliary nodes in ESnet topology
  number_of_nodes_in_random = [1]# set one if you use SURFnet or ESnet or Dumbbell
  edge_F_values = [1.0]# set of values to experiment for link level fidelity. we assume all links have this fidelity
  set_of_number_of_user_pairs = [3] # number of user pairs in the network
  lower_bound_distance_between_pairs = 200 # in km
  upper_bound_distance_between_pairs = 250 # in km
  square_dim_in_km_values = [100] # in km
  number_of_paths_values = [600] # number of paths between each pair of user pairs
  K_values = [1] # set of number of paths allowed to be used for each user pair. One means each user pair can use at most one path.
  lengths_of_middle_link = [100] # in km for dumbell shape topology
    
    
  
  """ network planning assumptions experiment hyperparameters """  
  dynamic_system_flag = False
  weighted_sum_rate_objective = False
  planning_with_wrong_weights = False
  weight_upper_bounds = [1]
  mu_values = [1] 


    
  """ repeaters and end nodes hyperparameters """
  D_values= [100] # set of values for repeaters memory budget
  repeater_min_decoherence = 100 # minimum repeater coherence time in microseconds
  repeater_max_decoherence= 1400#micro seconds
  repeater_decoherence_incresing_step = 50 # increaenting the coherence time of repeaters step size (used when experimenting with different values) 
  end_node_min_decoherence = 3000 # end nodes minimum coherence time in micoseconds 
  end_node_max_decoherence= 5000# micro seconds
  end_node_decoherence_incresing_step = 50
  checking_repeater_memory_life_time_flag = False # set to True of we want to restrict paths 
  checking_end_node_memory_life_time_flag = False 
  R_values = [10] # set of values for number of repetaers budget
  end_user_memory_set = [100] # the memory budget of end nodes
  
  
  
  

def get_config(FLAGS):
  config = Config
  for k, v in FLAGS.__flags.items():
    if hasattr(config, k):
      setattr(config, k, v.value)

  return config
