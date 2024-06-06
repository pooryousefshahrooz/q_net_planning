#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import print_function

import numpy as np
import random
import multiprocessing as mp
from absl import app
from absl import flags
import ast
from network import Network
import networkx as nx
from solver import Solver
from config import get_config
import time
import os
FLAGS = flags.FLAGS


# In[ ]:


utility_default_value = -200
config = get_config(FLAGS) or FLAGS
network = Network()
solver = Solver()
if config.include_fidelity_in_utility:
    network.include_fidelity_in_utility = True
else:
    network.include_fidelity_in_utility = False
network.relaxing_QCAST_formulation = config.relaxing_QCAST_formulation
global results_file_name
network.results_file_path = config.results_file_name
network.utility_type = config.utility_type
network.link_based_solver = config.link_based_solver
network.network_topology = config.network_topology
network.checking_repeater_memory_decoherence_flag = config.checking_repeater_memory_life_time_flag
network.checking_end_node_memory_decoherence_flag = config.checking_end_node_memory_life_time_flag
network.dynamic_system_flag = config.dynamic_system_flag
network.weighted_sum_rate_objective = config.weighted_sum_rate_objective
network.lower_bound_distance_between_pairs = config.lower_bound_distance_between_pairs
network.upper_bound_distance_between_pairs = config.upper_bound_distance_between_pairs
network.include_log = config.include_log
network.max_dist = config.max_dist
network.left_ESnet_subgraph = config.left_ESnet_subgraph
network.save_QCAST_equation_results_flag = config.save_QCAST_equation_results
network.toplogy_each_path_W_QCAST_equation_result = config.toplogy_each_path_W_QCAST_equation_result
if config.network_topology!="Random":
    graphs_with_nodes  = [1]
else:
    graphs_with_nodes = config.number_of_nodes_in_random

for round_number in range(config.repeating_times):
    for number_of_nodes_in_random in graphs_with_nodes:
        print("round ",round_number)
        for edge_F in config.edge_F_values:
            network.edge_F = edge_F# we assume all links have this fidelity
            for q_value in config.q_values:
                network.q  = q_value
                for number_of_user_pairs in config.set_of_number_of_user_pairs:
                    print("for number of user pairs %s from %s "%(number_of_user_pairs,config.set_of_number_of_user_pairs))

                    for end_user_memory in config.end_user_memory_set:
                        selected_user_pairs = []
                        network.end_user_memory = end_user_memory
                        for R in config.R_values:
                            network.R = R
                            for distance in config.lengths_of_middle_link:
                                network.reset_all_variables()
                                network.number_of_user_pairs = number_of_user_pairs
                                #for square_dim_in_km in square_dim_in_km_values:
                                network.metro_area_link_length = distance
                                network.L = distance # metro_area_link_length in km
                                network.N = 0
                                connected_graph_flag = True
                                if "Repeater_chain" in network.network_topology or "repeater_chain" in network.network_topology:
                                    network.generate_repeater_chain(distance,R)
                                elif "Dumbbell" in network.network_topology or "dumbbell" in network.network_topology:
                                    network.generate_dumbbell_shape_graph_with_n_repeater(distance,R)
                                elif "SURFnet" in network.network_topology or "Surfnet" in network.network_topology or "surfnet" in network.network_topology:
                                    network.parse_SURFnet(selected_user_pairs)
                                elif "ESnet" in network.network_topology or "esnet" in network.network_topology or "surfnet" in network.network_topology:
                                    selected_user_pairs = [("ORNL","ANL"),
                                                           ("SLAC","GA"),
                                                        ("NERSC", "INL"),
                                                        ("FNAL","Y12"),
                                                        ("PSFC","JLAB"),
                                                        ("KCNSC", "AMES"),
                                                        ]
                                    selected_user_pairs = [
                                                           ("SLAC","GA")

                                                        ]
                                    if network.left_ESnet_subgraph:
                                        selected_user_pairs = [
#                                             ("NETLPGH","PSFC"),
#                                               ("NETLMGN","PPPL"),
#                                               ("BNL","JLAB"),
                                            
                                            
                                            
                                              ("SRS","ORAU"),
                                              ("Y12","FNAL"),
                                              ("ORNL","ANL")
                                        ]
                                    else:
                                        selected_user_pairs = [
                                            ("NETLPGH","PSFC"),
                                              ("NETLMGN","PPPL"),
                                             ("BNL","JLAB"),
                                            
                                            
                                            
                                            
#                                               ("SRS","ORAU"),
#                                               ("Y12","FNAL"),
#                                               ("ORNL","ANL")
                                                              ]
                                        
                                    #network.parse_ESnet(selected_user_pairs)
                                    network.parse_ESnet_subgraph(selected_user_pairs)
                                    import pdb
                                    #pdb.set_trace()
                                elif "Random" in network.network_topology or "random" in network.network_topology:
                                    connected_graph_flag = network.generate_random_graph(number_of_nodes_in_random)
                                if connected_graph_flag:
    #                                 if round_number ==0:
    #                                     network.user_pairs = each_number_of_users_pairs[number_of_user_pairs]
                                    selected_user_pairs=network.user_pairs
                                    for pair in network.user_pairs:
                                        if pair not in selected_user_pairs:
                                            selected_user_pairs.append(pair)
                                    if network.dynamic_system_flag:
                                        if network.weighted_sum_rate_objective:
                                            random_weight = round(random.uniform(0.1,1),3)
    #                                         random_weight = 1
                                            random_user_poir = selected_user_pairs[random.randint(0,len(selected_user_pairs)-1)]
                                            network.each_user_pair_weight  ={}
                                            network.each_user_pair_weight[random_user_poir] = random_weight
                                            selected_weights = [random_weight]
                                            assigned_user_pairs = [random_user_poir]
                                            while(len(assigned_user_pairs)<len(selected_user_pairs)):
                                                random_weight = round(random.uniform(0.1,1),3)
                                                if random_weight not in selected_weights:
                                                    selected_weights.append(random_weight)
                                                    random_user_poir = selected_user_pairs[random.randint(0,len(selected_user_pairs)-1)]
                                                    if random_user_poir not in assigned_user_pairs:
                                                        network.each_user_pair_weight[random_user_poir] = random_weight
                                                        assigned_user_pairs.append(random_user_poir)

                                            selected_weights.sort(reverse=True)
                                            ordered_user_pairs = []
                                            for weight in selected_weights:
                                                for user_pair in assigned_user_pairs:
                                                    if network.each_user_pair_weight[user_pair]==weight and user_pair not in ordered_user_pairs:
                                                        ordered_user_pairs.append(user_pair)




                                    print("network.user_pairs ",network.user_pairs)

                                    for pair in network.user_pairs:
                                        if pair not in selected_user_pairs:
                                            selected_user_pairs.append(pair)
                                    network.square_dim_in_km = 0
                                    network.covered_W_edges_as_a_path=[]
                                    network.each_W_edges_E_F={}
                                    network.each_W_edges_E_F={}
                                    network.covered_W_path_edges_lengths_as_a_path = {}
                                    for cut_off in config.number_of_paths_values:
                                        network.reset_path_variables()
                                        network.number_of_paths = cut_off
                                        start_time = time.time()
                                        for scheme in config.schemes:
                                            network.scheme = scheme
                                            for D in config.D_values:
                                                network.D  = D
                                                if scheme =="Exhaustive_search" or "exhaustive" in scheme:
                                                    network.compute_paths_between_pairs()
                                                    """we enumerate all paths and compute their versions based on the value of W"""
                                                    network.M = 0
                                                    print("for scheme %s D %s R %s user pairs %s "%(scheme,
                                                                                           network.D,network.R,network.user_pairs))
                                                    for step_size in config.W_step_granularities:
                                                        this_round_starting_time = time.time()
                                                        network.set_of_paths = {}
                                                        network.each_user_pair_paths = {}
                                                        network.each_path_e_value  ={}
                                                        network.each_path_e2e_F_value = {}
                                                        start_time2 = time.time()
                                                        path_id = 0
                                                        for user_pair,path_ids in network.each_user_pair_dafualt_paths.items():    #                                         print("we have %s paths for user pair %s "%(len(path_ids),user_pair))
                                                            print("user pair is done ",user_pair)
                                                            set_of_path_ids = []
                                                            path_counter = 0
                                                            start_time2 = time.time()
                                                            if network.dynamic_system_flag and network.weighted_sum_rate_objective:

                                                                user_pair_weight = network.each_user_pair_weight[user_pair] 
                                                            else:
                                                                user_pair_weight = 1
                                                            for old_path_id  in path_ids:
                                                                current_time2 = time.time()
                                                                procssing_time2 = int(current_time2 -start_time2)
                                                                path_edges = network.set_of_paths_default[old_path_id]
                                                                start_time2 = time.time()

                                                                for edge in path_edges:
                                                                    try:
                                                                        p_value = network.transmission[edge]
                                                                    except:
                                                                        distance = nx.shortest_path_length(network.G, source=edge[0], target=edge[1],weight = "weight")
                                                                        p_value =  10**(-0.2*distance/10)
                                                                        network.transmission[edge] = p_value
                                                                    try:
                                                                        length = network.weights[edge]
                                                                    except:
                                                                        distance = nx.shortest_path_length(network.G, source=edge[0], target=edge[1],weight = "weight")
                                                                        network.weights[edge] = distance

                                                                retrived_flag,retrived_W,retrived_path_id,retrived_set_of_path_ids = network.extract_QCAST_equation_results(user_pair,path_edges)
                                                                for r_p_id in retrived_set_of_path_ids:
                                                                    set_of_path_ids.append(r_p_id)
                                                                path_id = path_id+retrived_path_id 
                                                                for W in range(max(1,retrived_W),min(D,network.end_user_memory)+1,step_size):

                                                                    start_time3 = time.time()
                                                                    E_p =network.compute_e2e_rate(path_id,path_edges,W)
                                                                    e2e_F = network.compute_e2e_fidleity(path_edges)
                                                                    current_time3 = time.time()
                                                                    procssing_time3 = int(current_time3 -start_time3)

                                                                    network.each_path_e_value[path_id] =E_p
                                                                    network.each_path_e2e_F_value[path_id] = e2e_F
                                                                    network.each_path_width[path_id] = W
                                                                    network.set_of_paths[path_id] = path_edges
                                                                    set_of_path_ids.append(path_id)
                                                                    network.each_path_weight[path_id] = user_pair_weight
                                                                    if network.save_QCAST_equation_results_flag:
                                                                        network.save_QCAST_equation_results(user_pair,path_edges,path_id,E_p,e2e_F,W,user_pair_weight)
                                                                    path_id+=1
                                                                path_counter+=1

                                                            network.each_user_pair_paths[user_pair] = set_of_path_ids
                                                        if not network.checking_end_node_memory_decoherence_flag:
                                                            config.end_node_decoherence_incresing_step = config.end_node_max_decoherence
                                                        for end_node_decoherence in range(config.end_node_min_decoherence,config.end_node_max_decoherence,config.end_node_decoherence_incresing_step):
                                                            if not network.checking_end_node_memory_decoherence_flag:
                                                                network.end_node_memory_decoherence_time = 10000000000
                                                            else:
                                                                network.end_node_memory_decoherence_time = end_node_decoherence

                                                            if not network.checking_repeater_memory_decoherence_flag:
                                                                config.repeater_decoherence_incresing_step = config.repeater_max_decoherence
                                                            for repeater_memory_decoherence_time in range(config.repeater_min_decoherence,config.repeater_max_decoherence,config.repeater_decoherence_incresing_step):
                                                                if not network.checking_repeater_memory_decoherence_flag:
                                                                    network.repeater_memory_decoherence_time = 100000000
                                                                else:
                                                                    network.repeater_memory_decoherence_time= repeater_memory_decoherence_time

                                                                network.set_paths_required_decoherence_time()
                                                                current_time_for_preparation_tiume = time.time()
                                                                preparation_tiume = current_time_for_preparation_tiume -start_time2 
                                                                for num_of_allowed_paths in config.K_values:
                                                                    network.K = num_of_allowed_paths
                                                                    try:
                                                                        static_utility,used_repeaters= solver.path_based_utility_maximization(False,network,[])
                                                                    except:
                                                                        static_utility = network.utility_default_value
                                                                    current_time = time.time()
                                                                    processing_time_in_seconds = int(current_time -this_round_starting_time)
                                                                    network.save(round_number,distance,
                                                                     static_utility,
                                                                    step_size,processing_time_in_seconds,
                                                                                 preparation_tiume,scheme,
                                                                                 number_of_user_pairs,"Static",-1,
                                                                                 number_of_nodes_in_random
                                                                                )

                                                                    if network.dynamic_system_flag:
                                                                        if config.planning_with_wrong_weights:
                                                                            if network.weighted_sum_rate_objective:
                                                                                random_weight = round(random.uniform(0.1,1),3)
                                        #                                         random_weight = 1
                                                                                random_user_poir = selected_user_pairs[random.randint(0,len(selected_user_pairs)-1)]
                                                                                network.each_user_pair_weight  ={}
                                                                                network.each_user_pair_weight[random_user_poir] = random_weight
                                                                                selected_weights = [random_weight]
                                                                                assigned_user_pairs = [random_user_poir]
                                                                                while(len(assigned_user_pairs)<len(selected_user_pairs)):
                                                                                    random_weight = round(random.uniform(0.1,1),3)
                                                                                    if random_weight not in selected_weights:
                                                                                        selected_weights.append(random_weight)
                                                                                        random_user_poir = selected_user_pairs[random.randint(0,len(selected_user_pairs)-1)]
                                                                                        if random_user_poir not in assigned_user_pairs:
                                                                                            network.each_user_pair_weight[random_user_poir] = random_weight
                                                                                            assigned_user_pairs.append(random_user_poir)

                                                                                selected_weights.sort(reverse=True)
                                                                                ordered_user_pairs = []
                                                                                for weight in selected_weights:
                                                                                    for user_pair in assigned_user_pairs:
                                                                                        if network.each_user_pair_weight[user_pair]==weight and user_pair not in ordered_user_pairs:
                                                                                            ordered_user_pairs.append(user_pair)




                                                                        each_list_optimal_utility = {}
                                                                        each_list_dynamic_utility = {}
                                                                        each_list_dynamic_utility_solution = {}
                                                                        each_list_optimal_utility_solution = {}
                                                                        back_up_to_store_all_paths_of_user_pairs = {}
                                                                        for user_pair,path_ids in network.each_user_pair_paths.items():
                                                                            back_up_to_store_all_paths_of_user_pairs[user_pair] = path_ids
                                                                        back_up_data_structure = {}
                                                                        for user_pair_p,each_p_edges in network.solution_each_user_pair_used_paths.items():
                                                                            back_up_data_structure[user_pair_p] = each_p_edges
                                                                        back_up_set_of_paths = {}
                                                                        for path_id,path_edges in network.set_of_paths.items():
                                                                            back_up_set_of_paths[path_id] = path_edges
                                                                        backup_user_pairs = []
                                                                        for user_pair in network.user_pairs:
                                                                            backup_user_pairs.append(user_pair)


                                                                        for mu in config.mu_values:
                                                                            events = network.generate_user_pairs_arrving_time(mu,300)
                                                                            for time_unit in range(1,int(max(events))): 

                                                                                utility = network.utility_default_value
                                                                                counter=0
                                                                                for event in events:
                                                                                    if event <= time_unit and event >time_unit-1 :
                                                                                        counter+=1
                                                                                counter = min(number_of_user_pairs,counter)
                                                                                arrived_user_pairs = []
                                                                                while(len(arrived_user_pairs)<counter):
                                                                                        random_user_pair = selected_user_pairs[random.randint(0,len(selected_user_pairs)-1)]
                                                                                        if random_user_pair not in arrived_user_pairs:
                                                                                            arrived_user_pairs.append(random_user_pair)

                                                                                network.user_pairs = arrived_user_pairs
                                                                                if tuple(arrived_user_pairs) in each_list_dynamic_utility:
                                                                                    dynamic_utility = each_list_dynamic_utility[tuple(arrived_user_pairs)]
                                                                                    network.solution_each_user_pair_used_paths = each_list_dynamic_utility_solution[tuple(arrived_user_pairs)]
                                                                                else:
                                                                                    network.each_user_pair_paths = {}
                                                                                    for user_pair_p,each_p_edges in back_up_data_structure.items():
                                                                                        user_pair = user_pair_p[0]
                                                                                        if user_pair in arrived_user_pairs:
                                                                                            path_id = user_pair_p[1]
                                                                                            path_edges = each_p_edges[0]

                                                                                            set_of_path_ids = []
                                                                                            path_counter = 0
                                                                                            start_time2 = time.time()

                                                                                            current_time2 = time.time()
                                                                                            procssing_time2 = int(current_time2 -start_time2)

                                                                                            start_time2 = time.time()


                                                                                            for W in range(1,min(D,network.end_user_memory)+1,step_size):
                                                                                                start_time3 = time.time()
                                                                                                E_p =network.compute_e2e_rate(path_id,path_edges,W)
                                                                                                e2e_F = network.compute_e2e_fidleity(path_edges)
                                                                                                current_time3 = time.time()
                                                                                                procssing_time3 = int(current_time3 -start_time3)

                                                                                                network.each_path_e_value[path_id] =E_p
                                                                                                network.each_path_e2e_F_value[path_id] = e2e_F
                                                                                                network.each_path_width[path_id] = W
                                                                                                network.set_of_paths[path_id] = path_edges
                                                                                                set_of_path_ids.append(path_id)
                                                                                                path_id+=1
                                                                                            path_counter+=1

                                                                                            network.each_user_pair_paths[user_pair] = set_of_path_ids
                                                                                    network.solution_each_user_pair_used_paths = {}
                                                                                    print("***********dynamic *********** going to solve dynamic")
                                                                                    try:
                                                                                        dynamic_utility,_= solver.path_based_utility_maximization(True,network,used_repeaters)
                                                                                    except ValueError:
                                                                                        print("ValueError",ValueError)
                                                                                        dynamic_utility = network.utility_default_value

                                                                                    each_list_dynamic_utility[tuple(arrived_user_pairs)] = dynamic_utility
                                                                                    each_list_dynamic_utility_solution[tuple(arrived_user_pairs)] = network.solution_each_user_pair_used_paths
                                                                                network.save(round_number,distance,
                                                                                 dynamic_utility,
                                                                                step_size,processing_time_in_seconds,
                                                                                             preparation_tiume,scheme,mu,
                                                                                             "Dynamic",time_unit,
                                                                                             number_of_nodes_in_random
                                                                                            )

                                                                                network.each_user_pair_paths= {}
                                                                                arrived_user_pairs.sort()
                                                                                if tuple(arrived_user_pairs) in each_list_optimal_utility:
                                                                                    optimal_utility = each_list_optimal_utility[tuple(arrived_user_pairs)]
                                                                                    network.solution_each_user_pair_used_paths = each_list_optimal_utility_solution[tuple(arrived_user_pairs)]

                                                                                else:
                                                                                    for user_pair,path_ids in back_up_to_store_all_paths_of_user_pairs.items():
                                                                                        if user_pair in arrived_user_pairs:
                                                                                            network.each_user_pair_paths[user_pair] = path_ids
                                                                                    try:
                                                                                        optimal_utility,_= solver.path_based_utility_maximization(False,network,[])
                                                                                    except ValueError:
                                                                                        optimal_utility = network.utility_default_value

                                                                                    each_list_optimal_utility[tuple(arrived_user_pairs)] = optimal_utility
                                                                                    each_list_optimal_utility_solution[tuple(arrived_user_pairs)] = network.solution_each_user_pair_used_paths
                                                                                network.save(round_number,distance,
                                                                                 optimal_utility,
                                                                                step_size,processing_time_in_seconds,
                                                                                             preparation_tiume,scheme,mu,
                                                                                             "Optimal",time_unit,
                                                                                             number_of_nodes_in_random
                                                                                            )





                                                elif scheme=="linear_link_based":
                                                    network.reset_path_variables()
                                                    network.set_of_paths = {}
                                                    network.each_user_pair_paths = {}
                                                    network.each_path_e_value  ={}
                                                    network.each_path_e2e_F_value = {}
                                                    network.path_id = 0
                                                    network.solution_repeater_counter = 0
                                                    network.K = 1
                                                    start_time3 = time.time()

                                                    try:
                                                        utility_value = solver.linear_link_based_in_CPLEX(network)
                                                        current_time = time.time()
                                                        processing_time_in_seconds = int(current_time -start_time3)
                                                        print("solved in %s seconds"%(processing_time_in_seconds))
                                                        network.save(round_number,distance,
                                                                 utility_value,0,processing_time_in_seconds,
                                                                     0,scheme,number_of_user_pairs,
                                                                     "Static",-1,
                                                                                 number_of_nodes_in_random
                                                                    )
                                                    except ValueError:
                                                        print(ValueError)
                                                else:
                                                    network.reset_path_variables()
                                                    network.set_of_paths = {}
                                                    network.each_user_pair_paths = {}
                                                    network.each_path_e_value  ={}
                                                    network.each_path_e2e_F_value = {}
                                                    network.path_id = 0
                                                    network.K = 1
                                                    start_time3 = time.time()
                                                    try:
                                                        utility_value = solver.link_based_formulation(network)
                                                        current_time = time.time()
                                                        processing_time_in_seconds = int(current_time -start_time3)
                                                        network.save(round_number,distance,
                                                                 utility_value,0,processing_time_in_seconds,
                                                                     0,scheme,number_of_user_pairs,
                                                                     "Static",-1,
                                                                                 number_of_nodes_in_random
                                                                    )
                                                    except:
                                                        pass





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




