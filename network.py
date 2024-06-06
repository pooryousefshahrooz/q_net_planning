#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import random
import matplotlib.pyplot as plt
import networkx as nx
from itertools import islice
import itertools
import math
import time
import csv
from math import *
import csv
import os
import sys
from absl import app
from absl import flags
import ast

from config import get_config
from solver import Solver
from docplex.mp.progress import *
from docplex.mp.progress import SolutionRecorder
import docplex.mp.model as cpx
# from pulp import LpMinimize, LpMaximize, LpProblem, LpStatus, lpSum, LpVariable, value, GLPK
import networkx as nx
import time
import matplotlib.pyplot as plt
# import pyomo.environ as pe
FLAGS = flags.FLAGS
OBJ_EPSILON = 1e-12


# In[32]:





# In[ ]:





# # A class that includes the codes related to generating topologies, computing paths and selecting repeater node candidtes etc. 

# In[ ]:





# In[ ]:


class Network:
    def __init__(self):
        self.N = 10# number of nodes
        self.L0 = 100 # square dim in km
        self.L = 100 # square dim in km
        self.network_topology = "Dumbbell"
        self.scheme ="Exhaustive_search"
        self.user_pairs = []
        self.G = nx.random_geometric_graph(10, 35, dim=2, p=2)
        self.weights = dict()
        self.transmission = dict()
        self.each_path_width = {}
        self.each_repeater_memory = {}
        self.utility_default_value = -1000
        self.upper_bound_distance_between_pairs = 0
        self.lower_bound_distance_between_pairs = 0
        self.p={}
        self.q = 1.0
        self.repeater_places = []
        self.set_of_paths_default = {}
        self.utility_values = set([])
        self.each_user_pair_dafualt_paths ={}
        self.each_path_e_value ={}
        self.each_path_e2e_F_value={}
        self.each_path_width = {}
        self.path_counter = 0
        self.set_of_paths = {}
        self.multi_path_state = False
        self.save_QCAST_equation_results_flag = False
        self.set_of_edges = []
        self.num_of_allowed_paths = 1
        self.metro_area_link_length = 1
        self.include_fidelity_in_utility  = True
        self.relaxing_QCAST_formulation = False
        self.edge_F = 0.99
        self.forbidden_edges = []
        self.checking_memory_decoherence_flag = False
        self.repeater_memory_decoherence_time= 100000000
        self.end_node_memory_decoherence_time = 100000000
        self.each_path_repeater_required_decoherence_time = {}
        self.each_path_end_node_required_decoherence_time = {}
        self.solution_each_user_pair_used_paths = {}
        
        self.checking_repeater_memory_decoherence_flag = False
        self.checking_end_node_memory_decoherence_flag = False
    
        self.left_ESnet_subgraph = False
        self.toplogy_each_path_W_QCAST_equation_result = ""
        
        self.dynamic_system_flag = False
        self.weighted_sum_rate_objective = False
        self.include_log = True
        self.each_path_weight = {}
        
    def split(self,x, n):
        points = []

        if (x % n == 0):
            for i in range(n):
                points.append(x//n)
        else:

            zp = n - (x % n)
            pp = x//n
            for i in range(n):
                if(i>= zp):
                    points.append(pp + 1)
                else:
                    points.append(pp)
        return points    
    def check_a_valid_path(self,path_edges):
        """ checking if a path is not using any forbidden edges"""
        e2e_F = self.compute_e2e_fidleity(path_edges)
        for edge in path_edges:
            if edge in self.forbidden_edges:
                return False
        return True
    def compute_paths_between_pairs(self):
        """computes the shortest number_of_paths paths between user pairs"""
        path_counter = 0
        for pair in self.user_pairs:
            start_time = time.time()
            k_paths_user_pair1  = self.get_paths_between_user_pairs(pair,self.number_of_paths)
            k_paths_user_pair1.append([pair[0],pair[1]])
            user_pair_path_ids = []
            for path in k_paths_user_pair1:
                node_indx = 0
                path_edges = []
                for node_indx in range(len(path)-1):                
                    path_edges.append((path[node_indx],path[node_indx+1]))
                    node_indx+=1
                if self.check_a_valid_path(path_edges):
#                     print("for user pair %s we are adding path %s "%(pair,path_edges))
                    self.set_of_paths_default[path_counter] = path_edges
                    user_pair_path_ids.append(path_counter)
                    path_counter+=1
            self.each_user_pair_dafualt_paths[pair] = user_pair_path_ids
            end_time = time.time()
            duration = end_time-start_time
    def get_paths_between_user_pairs(self,user_pair,cut_off_for_path_searching):
        return self.k_shortest_paths(user_pair[0], user_pair[1], cut_off_for_path_searching,"weight")
    def k_shortest_paths(self,source, target, k, weight):
        return list(
           
                islice(nx.shortest_simple_paths(self.G, source, target, weight=weight), k)
           
        )

    
    def compute_e2e_fidleity(self,path_edges):
        """computing the end to end fidelity of a path"""
        if path_edges:
            F_product = (4*self.edge_F-1)/3 
            for edge in path_edges[1:]:
                F_product  = F_product*(((4*self.edge_F)-1)/3)
        else:
            print("Error")
            return 1.0
        N = len(path_edges)
        p1 = 1
        p2 = 1
        F_final = (1/4)+(3/4)*((p1*(((4*(p2**2))-1)/3))**(N-1))*(F_product)

        return round(F_final,3)
    

    def f(self,W,P,Q,i,k):
        """computing the QCAST recursive function in Shi, Shouqian, and Chen Qian SIGCOM202 paper"""
        if k==1:
            return Q[i][k]
        else:
            sum_Q = 0
            for l in range(i,W+1):
                sum_Q = sum_Q+Q[l][k]
            sum_P = 0
            for l in range(i+1,W+1):
                sum_P = sum_P+self.f(W,P,Q,l,k-1)
            return self.f(W,P,Q,i,k-1)* sum_Q+Q[i][k]* sum_P

    def compute_e2e_rate(self,path_id,path,W):
        """using QCAST or approximate function to compute end to end rate on a path"""
        if (self.relaxing_QCAST_formulation or W>=100):
            p_values = []
            for edge in path:
                try:
                    p_value = self.transmission[edge]
                except:
                    distance = nx.shortest_path_length(self.G, source=edge[0], target=edge[1],weight = "weight")
                    p_value =  10**(-0.2*distance/10)
                    self.transmission[edge] = p_value
                    self.weights[edge] = distance
                p_values.append(p_value)
            return W * min(p_values) * (self.q**(len(path)-1))
        else:
            scheme_key = str(W)+"-path"
            h =len(path)
            P={}
            Q={}
            q_value = self.q
            self.p={}
            each_scheme_each_point_value = {}
            edge_indx = 1

            for edge in path:
                try:
                    p_value = self.transmission[edge]
                except:
                    distance = nx.shortest_path_length(self.G, source=edge[0], target=edge[1],weight = "weight")
                    p_value =  10**(-0.2*distance/10)
                    self.transmission[edge] = p_value
                    self.weights[edge] = distance
                self.p[edge_indx]=p_value
                edge_indx+=1
            for i in range(1,W+1):
                for k in range(1,h+1):
                    try:
                        Q[i][k] =math.comb(W, i)*self.p[k]**i * (1-self.p[k])**(W-i)
                    except:
                        Q[i]={}
                        Q[i][k] =math.comb(W, i)*self.p[k]**i * (1-self.p[k])**(W-i)

                    if k==1:
                        try:
                            P[i][1] = Q[i][1]
                        except:
                            P[i]={}
                            P[i][1] = Q[i][1]

            for i in range(1,W+1):
                for k in range(1,h+1):
                    value = self.f(W,P,Q,i,k)
                    P[i][k] = value
            sum_i = 0
            for i in range(1,W+1):
                sum_i = sum_i+(i * P[i][k])
            EXT = q_value**(h-1) *sum_i
            return EXT

    def save_QCAST_equation_results(self,user_pair,edges,path_id,E_p,e2e_F,W,user_pair_weight):
        with open(self.toplogy_each_path_W_QCAST_equation_result, 'a') as newFile:                                
            newFileWriter = csv.writer(newFile)
            edges_string = ""
            for edge in edges:
                if edges_string:
                    edges_string = edges_string+","+str(edge[0])+":"+str(edge[1])
                else:
                    edges_string = str(edge[0])+":"+str(edge[1])
            
            user_pair_string = str(user_pair[0])+":"+str(user_pair[1])

            newFileWriter.writerow([self.network_topology,user_pair_string,
                                    path_id,edges_string,
                                    E_p,e2e_F,
                                    W,user_pair_weight
                                   ])
    def extract_QCAST_equation_results(self,user_pair,path_edges):
        set_of_path_ids = []
        exist_flag =False
        W = 1
        retrived_path_id = 0
        try:
            with open(self.toplogy_each_path_W_QCAST_equation_result, "r") as f:
                reader = csv.reader( (line.replace('\0','') for line in f) )
                for line in reader:
                    network = line[0]
                    saved_user_pair = (line[1]).split(":")
                    saved_user_pair = ((saved_user_pair[0]),(saved_user_pair[1]))
                    if saved_user_pair==user_pair and network==self.network_topology and saved_path_edges == path_edges:
                        exist_flag =True
                        retrived_path_id = int(line[2])
                        saved_edges = (line[3])
                        saved_edges = saved_edges.split(",")
                        saved_path_edges =[]
                        for edge in saved_edges:
                            edge = edge.split(":")
                            saved_path_edges.append(((edge[0]),(edge[1])))
                        E_p =float(line[4])
                        e2e_F = float(line[5])
                        W = int(line[6])
                        user_pair_weight = line[7]
                        self.each_path_e_value[saved_path_id] =E_p
                        self.each_path_e2e_F_value[saved_path_id] = e2e_F
                        self.each_path_width[saved_path_id] = W
                        self.set_of_paths[saved_path_id] = saved_path_edges
                        set_of_path_ids.append(saved_path_id)
                        self.each_path_weight[saved_path_id] = user_pair_weight
            network.each_user_pair_paths[saved_user_pair] = set_of_path_ids
        except:
            pass
        return exist_flag,W,retrived_path_id,set_of_path_ids
    def check_if_path_using_any_end_nodes(self,path_edges,given_user_pair):
        for user_pair in self.user_pairs:
            if user_pair != given_user_pair:
                for edge in path_edges:
                    if edge[0] in user_pair or edge[1] in user_pair:
                        return True
        return False
    def check_path_uses_repeater(self,u,p):
        for edge in self.set_of_paths[p]:
            if u in edge:
                return True
        return False


    def compute_equation(self,path_id,e2e_rate,e2e_F):
        if e2e_rate==0.0 or e2e_F==0:
            return self.utility_default_value
        else:
            if self.utility_type=="DE":
                if self.include_fidelity_in_utility:
                    return np.log(e2e_rate*(e2e_F-0.5))
                else:
                    return math.log(e2e_rate)
            elif self.utility_type=="SKF":
                if self.include_fidelity_in_utility:
                    return np.log(e2e_rate*(e2e_F-0.5))
                else:
                    return np.log(e2e_rate)
            else:
                if self.include_fidelity_in_utility:
                    if not self.include_log:
                        return e2e_rate*(e2e_F-0.5)
                    else:
                        return  np.log(e2e_rate*(e2e_F-0.5))
                else:
                        return np.log(e2e_rate)
    def compute_utility(self,path_id):
        e_value = self.each_path_e_value[path_id]
        e2e_F_value =self.each_path_e2e_F_value[path_id]
        if self.each_path_e_value[path_id]>0 and e2e_F_value>0.5:
            utility_on_this_path = self.compute_equation(path_id,e_value,e2e_F_value)
        else:
            return self.utility_default_value
        return utility_on_this_path

    
    
    def reset_all_variables(self):
        self.weights = dict()
        self.transmission = dict()
        self.each_path_width = {}
        self.each_repeater_memory = {}
        self.utility_default_value = -1000
        self.p={}
        self.repeater_places = []
        self.set_of_paths_default = {}
        self.each_user_pair_dafualt_paths ={}
        self.each_path_e_value ={}
        self.each_path_e2e_F_value={}
        self.each_path_width = {}
        self.path_counter = 0
        self.set_of_paths = {}
        self.multi_path_state = False
        self.set_of_edges = []
        self.num_of_allowed_paths = 1
    def reset_path_variables(self):
        self.set_of_paths_default = {}
        self.each_user_pair_dafualt_paths ={}
        self.each_path_e_value ={}
        self.each_path_e2e_F_value={}
        self.each_path_width = {}
        self.path_counter = 0
        self.set_of_paths = {}
        
    def save(self,round_number,l_of_middle_link,reported_max_utility,
            step_size,processing_time_in_seconds,preparation_time,scheme,mu,system_model,time_unit,
            number_of_nodes_in_graph):
        if self.include_fidelity_in_utility:
            include_fidelity_in_utility = "True"
        else:
            include_fidelity_in_utility = "False"
            
        if self.relaxing_QCAST_formulation:
            relaxing_QCAST_formulation_key = "True"
        else:
            relaxing_QCAST_formulation_key = "False"
        if scheme=="link_based" or scheme =="linear_link_based":
            computed_max_utility = 0
            for user_pair_p,each_p_edges in self.solution_each_user_pair_used_paths.items():
                user_pair = user_pair_p[0]
                used_path = user_pair_p[1]
                edges = each_p_edges
                W_p = self.each_path_width[used_path]
                E_p = self.each_path_e_value[used_path]
                e2e_F  = self.each_path_e2e_F_value[used_path]
                computed_max_utility = computed_max_utility + self.compute_utility(used_path)
        else:
            computed_max_utility = reported_max_utility
        with open(self.results_file_path, 'a') as newFile:                                
            newFileWriter = csv.writer(newFile)
            if self.solution_each_user_pair_used_paths:
                for user_pair_p,each_p_edges in self.solution_each_user_pair_used_paths.items():
                    user_pair = user_pair_p[0]
                    used_path = user_pair_p[1]
                    edges = each_p_edges[0]
                    W_p = self.each_path_width[used_path]
                    E_p = self.each_path_e_value[used_path]
                    e2e_F  = self.each_path_e2e_F_value[used_path]
                    p_values = []
                    for edge in edges:
                        try:
                            p_value = self.transmission[edge]
                        except:
                            distance = nx.shortest_path_length(self.G, source=edge[0], target=edge[1],weight = "weight")
                            p_value =  10**(-0.2*distance/10)
                            self.transmission[edge] = p_value
                        
                        p_values.append(p_value)

                    e2e_rate = W_p * min(p_values) * (self.q**(len(edges)-1))
                    if self.weighted_sum_rate_objective:
                        utility = e2e_rate * self.each_path_weight[used_path]
                    else:
                        try:
                            utility  = np.log(e2e_rate)
                        except:
                            utility = self.utility_default_value
                    computed_again = self.compute_utility(used_path)
                    

                    newFileWriter.writerow([round_number,self.N,
                                            self.number_of_user_pairs,
                                            l_of_middle_link,
                                            self.square_dim_in_km,
                                            self.number_of_paths,
                                            self.K,
                                            self.D,0,self.end_user_memory,
                                            self.R,reported_max_utility,
                                            user_pair,self.solution_repeater_counter,

                                            self.solution_used_nodes_by_paths,

                                            self.edge_F,self.solution_path_counter,
                                           used_path,W_p,E_p,e2e_F,edges,min(p_values),
                                           include_fidelity_in_utility,
                                           self.q,step_size,processing_time_in_seconds,
                                           relaxing_QCAST_formulation_key,
                                            preparation_time,scheme,#29
                                            self.network_topology,
                                           self.link_based_solver,
                                           self.repeater_memory_decoherence_time,
                                           self.end_node_memory_decoherence_time,
                                           computed_max_utility,
                                           mu,system_model,#36
                                           self.dynamic_system_flag,
                                            self.weighted_sum_rate_objective,
                                           time_unit,#39
                                           number_of_nodes_in_graph,
                                           self.lower_bound_distance_between_pairs,
                                           self.upper_bound_distance_between_pairs,
                                           self.max_dist])
            else:
                newFileWriter.writerow([round_number,self.N,
                                    self.number_of_user_pairs,
                                    l_of_middle_link,
                                    self.square_dim_in_km,
                                    self.number_of_paths,
                                    self.K,
                                    self.D,0,self.end_user_memory,
                                    self.R,reported_max_utility,
                                    (0,0),0,

                                    0,

                                    self.edge_F,0,
                                   0,0,0,0,"edges",0,
                                   include_fidelity_in_utility,
                                   self.q,step_size,processing_time_in_seconds,
                                   relaxing_QCAST_formulation_key,
                                    preparation_time,scheme,
                                    self.network_topology,
                                   self.link_based_solver,
                                   self.repeater_memory_decoherence_time,
                                   self.end_node_memory_decoherence_time,
                                   reported_max_utility,
                                       mu,system_model,
                                       self.dynamic_system_flag,
                                            self.weighted_sum_rate_objective,
                                       time_unit,
                                       number_of_nodes_in_graph,
                                       self.lower_bound_distance_between_pairs,
                                       self.upper_bound_distance_between_pairs,
                                       self.max_dist])
                
                
    def generate_repeater_chain(self,dist,number_of_repetaers):
        self.G = nx.Graph()
        self.user_pairs = []
        self.repeater_places=[]
        self.weights = {}
        self.set_of_edges=[]
        self.transmission = {}
        self.each_pair_id_nodes ={}
        self.each_edge_lenght={}
        self.pair_id = 0
        end_nodes = []
        for i in range(number_of_repetaers+1):
            source = i
            end = i+1
            e = (source,end)
            self.G.add_edge(source,end,weight=dist)
            p = 10**(-0.2*dist/10)
            self.weights[e] = dist
            self.weights[(e[1],e[0])] = dist
            self.transmission[e] = 10**(-0.2*dist/10)
        self.pos = nx.spring_layout(self.G)
        self.user_pairs.append((0,number_of_repetaers+1))
        user_pair_nodes = [0,number_of_repetaers+1]
        for node in self.G.nodes:
            if node not in user_pair_nodes:
                self.repeater_places.append(node)
        self.N =len(self.G.nodes)
        self.node_list = np.arange(self.N)
        length = nx.shortest_path_length(self.G, source=0, target=number_of_repetaers+1,weight = "weight")
        
    def generate_dumbbell_shape_graph_with_n_repeater(self,distance,R):
        print("generating a dumble shape topology with %s user pairs and middle link %s divided by R repetaer %s "%(self.number_of_user_pairs,
                                                                                 distance,
                                                                               R))
        def check_not_on_same_side(node1,node2,first_side_nodes,second_side_nodes):
            """this function checked that one one on one side of dumble
            does not connect a node on the same side of the dumble"""
            if node1 in first_side_nodes and node2 in first_side_nodes:
                return False
            elif node1 in second_side_nodes and node2 in second_side_nodes:
                return False
            return True
        
        self.G = nx.Graph()
        self.user_pairs = []
        self.repeater_places=[]
        self.weights = {}
        self.set_of_edges=[]
        self.transmission = {}
        self.each_pair_id_nodes ={}
        self.each_edge_lenght={}
        
        """these are the edges that connect a node on one side of dumble to the another node on the same side of dumble"""
        self.forbidden_edges = []
        
        self.pair_id = 0
        end_nodes = []
        self.metro_area_link_length = distance
        
            
        first_side_nodes = []
        second_side_nodes = []
            
        for n in range(self.number_of_user_pairs):
            first_side_nodes.append(n)
        for n in range(self.number_of_user_pairs):
            source = n
            end = self.number_of_user_pairs
            e = (n,self.number_of_user_pairs)
            
            dist = distance
            p = 10**(-0.2*dist/10)
            self.forbidden_edges.append((e[1],e[0]))
            self.each_pair_id_nodes[self.pair_id] = [source]
            end_nodes.append(source)
            self.pair_id+=1
            self.weights[e] = dist
            self.weights[(e[1],e[0])] = dist
            self.transmission[e] = 10**(-0.2*dist/10)
#             print("distance between %s and %s is %s with p %s"%(e[0],e[1],dist,self.transmission[e]))
            self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)
            self.G.add_edge(source,end,weight=dist)


        end = self.number_of_user_pairs
        new_repeaters_id = self.number_of_user_pairs+1
        edge = (self.number_of_user_pairs,self.number_of_user_pairs+1)
        middle_flag = False
        distance_between_repeaters = distance/(R-1)
        for i in range(R-1):
            source = end
            end = new_repeaters_id
            e = (source,end)
            dist = distance
            p = 10**(-0.2*dist/10)
            self.weights[e] = dist
            self.weights[(e[1],e[0])] = dist
            self.transmission[e] = 10**(-0.2*dist/10)
            self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)
            self.G.add_edge(source,end,weight=dist)
            new_repeaters_id+=1
        for n in range(new_repeaters_id,new_repeaters_id+self.number_of_user_pairs):
            second_side_nodes.append(n)
                        
        source = new_repeaters_id-1
        self.pair_id = 0
        for n in range(self.number_of_user_pairs):
            end = new_repeaters_id
            e = (source,end)
            dist = distance
            p = 10**(-0.2*dist/10)
            self.forbidden_edges.append((e[1],e[0]))
            self.each_pair_id_nodes[self.pair_id].append(end)
            end_nodes.append(end)
            self.pair_id +=1
            self.weights[e] = dist
            self.weights[(e[1],e[0])] = dist
            self.transmission[e] = p
            self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)
            self.G.add_edge(source,end,weight=dist)
            new_repeaters_id+=1
#         nx.set_edge_attributes(self.G, values = self.weights, name = 'weight')
#         nx.set_edge_attributes(self.G, values = self.transmission, name = 'trans')

        # Get the layout
        self.pos = nx.spring_layout(self.G)
        
        
        labels = nx.get_edge_attributes(self.G,'weight')
        nx.draw_networkx_edge_labels(self.G,self.pos,edge_labels=labels)
#         plt.figure(figsize=(4,3))
#         nx.draw(self.G,with_labels=True,edge_labels=labels)
        plt.show()
#         nx.draw(self.G)
#         labels = nx.get_edge_attributes(self.G,'weight')
#         nx.draw_networkx_edge_labels(self.G,self.pos,edge_labels=labels)
        plt.savefig("plotting/plots/Graph_for_N_"+str(self.number_of_user_pairs)+"_R_"+str(R)+"_L"+str(distance)+".png", format="PNG")
        

        user_pair_nodes = []
        for pair,nodes in self.each_pair_id_nodes.items():
            self.user_pairs.append((nodes[0],nodes[1]))
            user_pair_nodes.append(nodes[0])
            user_pair_nodes.append(nodes[1])
        for node in self.G.nodes:
            if node not in user_pair_nodes:
                self.repeater_places.append(node)
        self.N =len(self.G.nodes)
        self.node_list = np.arange(self.N)
        length = nx.shortest_path_length(self.G, source=self.user_pairs[0][0], target=self.user_pairs[0][1],weight = "weight")
    
    def read_graph_from_gml(self,file):
        G = nx.read_gml(file)
        pos = {}
        for node, nodedata in G.nodes.items():
            if "position" in nodedata:
                pos[node] = ast.literal_eval(nodedata["position"])
            elif "Longitude" in nodedata and "Latitude" in nodedata:
                pos[node] = [nodedata['Longitude'], nodedata['Latitude']]
            else:
                raise ValueError("Cannot determine node position.")
        nx.set_node_attributes(G, pos, name='pos')
        return G

    def compute_dist_lat_lon(self,edge,graph):
        """Compute the distance in km between two points based on their latitude and longitude.
        Assumes both are given in radians."""
        R = 6371  # Radius of the earth in km
        node1, node2 = edge
        lon1 = np.radians(graph.nodes[node1]['Longitude'])
        lon2 = np.radians(graph.nodes[node2]['Longitude'])
        lat1 = np.radians(graph.nodes[node1]['Latitude'])
        lat2 = np.radians(graph.nodes[node2]['Latitude'])
        delta_lat = lat2 - lat1
        delta_lon = lon2 - lon1
        a = np.sin(delta_lat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin(delta_lon / 2) ** 2)
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return np.round(R * c, 5)
    
    def parse_ESnet(self,selected_user_pairs):
        self.G = nx.Graph()
        self.user_pairs = []
        self.repeater_places=[]
        self.weights = {}
        self.set_of_edges=[]
        self.transmission = {}
        self.each_pair_id_nodes ={}
        self.each_edge_lenght={}
        """these are the edges that connect a node on one side of dumble to the another node on the same side of dumble"""
        self.forbidden_edges = []
        self.pair_id = 0
        end_nodes = []
        self.G = self.read_graph_from_gml("data/"+self.network_topology)
        pos_list = nx.get_node_attributes(self.G, 'pos')
        # nx.draw_networkx(G=G,  width=1, node_size = 20, with_labels=False, font_size=6, pos=pos_list)
        # states.boundary.plot()
        max_dist = self.max_dist

        orig_edges = np.copy(self.G.edges())
        c_edge = 0
        c_rep = 0
        self.weights = dict()
        for edge in orig_edges:
            dist = self.compute_dist_lat_lon(edge,self.G)
            if dist > max_dist:
                n1, n2 = edge
                self.G.remove_edge(n1,n2)
                c_edge += 1
                n_rep = int(dist/max_dist)
                self.G.add_nodes_from(np.arange(c_rep,c_rep+n_rep))
                self.G.add_edge(n1,c_rep)
                self.weights[(n1,c_rep)] = dist/(n_rep+1)                
                self.weights[(c_rep,n1)] = dist
                self.transmission[(n1,c_rep)] = 10**(-0.2*dist/10)
                self.transmission[(c_rep,n1)] = 10**(-0.2*dist/10)
                

                i_r = 0
                pos_list[(i_r+c_rep)] = list(np.array(pos_list[n1])+ (np.array(pos_list[n2])-np.array(pos_list[n1]))/(n_rep+1)*(i_r+1))
                self.G.add_edge(c_rep+n_rep-1,n2)
                self.weights[(c_rep+n_rep-1,n2)] = dist/(n_rep+1)
                self.weights[(n2,c_rep+n_rep-1)] = dist
                self.transmission[(c_rep+n_rep-1,n2)] = 10**(-0.2*dist/10)
                self.transmission[(n2,c_rep+n_rep-1)] = 10**(-0.2*dist/10)
                
                i_r = n_rep -1
                pos_list[(i_r+c_rep)] = list(np.array(pos_list[n1])+ (np.array(pos_list[n2])-np.array(pos_list[n1]))/(n_rep+1)*(i_r+1))
                for i_r in range(0,n_rep-1):
                    self.G.add_edge(i_r+c_rep,i_r+c_rep+1)
                    pos_list[(i_r+c_rep)] = list(np.array(pos_list[n1])+ (np.array(pos_list[n2])-np.array(pos_list[n1]))/(n_rep+1)*(i_r+1))
                    self.weights[(i_r+c_rep,i_r+c_rep+1)] = dist/(n_rep+1)
                    self.weights[(i_r+c_rep+1,i_r+c_rep)] = dist
                    self.transmission[(i_r+c_rep,i_r+c_rep+1)] = 10**(-0.2*dist/10)
                    self.transmission[(i_r+c_rep+1,i_r+c_rep)] = 10**(-0.2*dist/10)

                c_rep += n_rep
            else:
                n1, n2 = edge
                self.weights[(n1,n2)] = dist
                self.weights[(edge[1],edge[0])] = dist
                self.transmission[(edge[0],edge[1])] = 10**(-0.2*dist/10)
                self.transmission[(edge[1],edge[0])] = 10**(-0.2*dist/10)

        nx.set_edge_attributes(self.G, values = self.weights, name = 'weight')



        self.pos = nx.spring_layout(self.G)
        labels = nx.get_edge_attributes(self.G,'weight')
        nx.draw_networkx_edge_labels(self.G,self.pos,edge_labels=labels)
#         plt.figure(figsize=(4,3))
#         nx.draw(self.G,with_labels=True,edge_labels=labels)
#         plt.show()
        nx.draw(self.G)
        plt.savefig("plotting/plots/ESnet_Graph.png", format="PNG")
        node_counter = 0
        user_pairs = []
        for node1 in self.G.nodes:
            node_counter+=1
            if node1 not in end_nodes:
                if node1 not in self.repeater_places:
                    self.repeater_places.append(node1)
            for node2 in self.G.nodes:
                if node1!=node2:
                    e = (node1,node2)
                    if e not in self.G.edges:
                        shortest_path = self.get_paths_between_user_pairs(e,1)
                        path_edges = []
                        node_indx = 0
                        shortest_path = shortest_path[0]
                        for node_indx in range(len(shortest_path)-1):                
                            path_edges.append((shortest_path[node_indx],shortest_path[node_indx+1]))
                            node_indx+=1
                        dist = 0
                        for edge in path_edges:
                            dist = dist+self.weights[edge]

                        user_pairs.append((node1,node2))
                        self.weights[e] = dist
                        self.weights[(e[1],e[0])] = dist
                        self.G.add_edge(node1,node2,weight=dist)
                        self.transmission[e] = 10**(-0.2*dist/10)
                        p = 10**(-0.2*dist/10)
                        self.set_of_edges.append(e)
                        self.set_of_edges.append((e[1],e[0]))
                        self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)
        for pair in selected_user_pairs:
            self.user_pairs.append(pair)
        user_pairs_nodes = []
        self.repeater_places = []
        end_node_list = ["PNNL","LIGO","HLAN","NETLALB","INL","SNLL",\
                        "LBNL","NERSC","SLAC","LLNL","GA","NNSS","LANL-NNSS",\
                        "NREL","LANL","SNLA","KCNSC-NM","NGA-SW","KCNSC",\
                        "PANTEX","ORNL","DOE-SC-OSTI","FNAL","AMES","ANL",\
                        "ORCC","Y12","SRS","ORAU","LNS","PSFC","NREL-DC",\
                        "BNL","JLAB","PPPL","NETLPGH","NETLMGN"]

        color_map = []
        for node, nodedata in self.G.nodes.items():
            if type(node)==np.int64:
                nodedata['type'] = 'auxiliary'
                self.repeater_places.append(node)
            else:
                if node in end_node_list:
                    nodedata['type'] = 'site'
                else:
                    nodedata['type'] = 'router'
                    self.repeater_places.append(node)
        
        
        self.N =len(self.G.nodes)
        self.node_list = np.arange(self.N)
                                
    def parse_SURFnet(self,selected_user_pairs):
        self.G = nx.Graph()
        self.user_pairs = []
        self.repeater_places=[]
        self.weights = {}
        self.set_of_edges=[]
        self.transmission = {}
        self.each_pair_id_nodes ={}
        self.each_edge_lenght={}
        """these are the edges that connect a node on one side of dumble to the another node on the same side of dumble"""
        self.forbidden_edges = []
        self.pair_id = 0
        end_nodes = []
        print('[*] Loading topology...',self.network_topology)
        self.G = nx.read_gml("data/"+self.network_topology, label='id', destringizer=None)
        for s, d, dist in self.G.edges(data=True):
            dist = float(dist["length"])
            e = (int(s),int(d))
            source = int(s)
            end = int(d)
            self.weights[e] = dist
            self.weights[(e[1],e[0])] = dist
            self.transmission[e] = 10**(-0.2*dist/10)
            self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)
            self.G.add_edge(source,end,weight=dist)
        #f.close()
            
        # Get the layout
        self.pos = nx.spring_layout(self.G)
        
        
        labels = nx.get_edge_attributes(self.G,'weight')
        nx.draw_networkx_edge_labels(self.G,self.pos,edge_labels=labels)
#         plt.figure(figsize=(4,3))
#         nx.draw(self.G,with_labels=True,edge_labels=labels)
#         plt.show()
        nx.draw(self.G)

        user_pairs = []
        for node1 in self.G.nodes:
            if node1 not in end_nodes:
                if node1 not in self.repeater_places:
                    self.repeater_places.append(node1)
            for node2 in self.G.nodes:
                if node1!=node2:
                    e = (node1,node2)
                    if e not in self.G.edges:
                        shortest_path = self.get_paths_between_user_pairs(e,1)
                        path_edges = []
                        node_indx = 0
                        shortest_path = shortest_path[0]
                        for node_indx in range(len(shortest_path)-1):                
                            path_edges.append((shortest_path[node_indx],shortest_path[node_indx+1]))
                            node_indx+=1
                        dist = 0
                        for edge in path_edges:
                            dist = dist+self.weights[edge]
                        if dist<=self.upper_bound_distance_between_pairs and dist>=self.lower_bound_distance_between_pairs:
                            user_pairs.append((node1,node2))
                        self.weights[e] = dist
                        self.weights[(e[1],e[0])] = dist
                        self.G.add_edge(node1,node2,weight=dist)
                        self.transmission[e] = 10**(-0.2*dist/10)
                        p = 10**(-0.2*dist/10)
                        self.set_of_edges.append(e)
                        self.set_of_edges.append((e[1],e[0]))
                        self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)
        
        for pair in selected_user_pairs:
            self.user_pairs.append(pair)
        while(len(self.user_pairs)<self.number_of_user_pairs):
            random_user_pair = user_pairs[random.randint(0,len(user_pairs)-1)]
            if random_user_pair not in self.user_pairs and (random_user_pair[1],random_user_pair[0]) not in self.user_pairs:
                self.user_pairs.append((random_user_pair[0],random_user_pair[1]))
        
        user_pairs_nodes = []
        self.repeater_places = []
        for pair in self.user_pairs:
            user_pairs_nodes.append(pair[0])
            user_pairs_nodes.append(pair[1])
        for node1 in self.G.nodes:
            if node1 not in  self.repeater_places and node1 not in user_pairs_nodes:
                self.repeater_places.append(node1)
        self.N =len(self.G.nodes)
        self.node_list = np.arange(self.N)
                                
    def generate_user_pairs_arrving_time(self,mu,number_of_events):
        #generate number_of_events events with a mu events per unit time
        import random
        n=number_of_events
        time_span=n/mu
        events=[]
        for j in range(0,n):
            events.append(random.random())
        events.sort()

        for j in range(0,n)  :
            events[j]*=time_span

        return events
        
    def generate_random_graph(self,number_of_nodes):
        print("we are going to generate a random topology with %s nodes"%(number_of_nodes))
        self.G = nx.Graph()
        self.user_pairs = []
        self.repeater_places=[]
        self.weights = {}
        self.set_of_edges=[]
        self.transmission = {}
        self.each_pair_id_nodes ={}
        self.each_edge_lenght={}
        """these are the edges that connect a node on one side of dumble to the another node on the same side of dumble"""
        self.forbidden_edges = []
        self.pair_id = 0
        end_nodes = []
        # random network graph of N nodes
        N = 14# number of nodes
        N = number_of_nodes# number of nodes
        node_list = np.arange(N)
        L0 = self.L # square dim in km
        length = L0
        width = L0

        dmax = L0
        random.seed()
        pos = dict()
        for node in range(N):
            pos[node] = (random.random()*length,random.random()*width)

        self.G = nx.random_geometric_graph(N, dmax, dim=2, pos=pos, p=2)

        self.weights = dict()
        self.transmission = dict()
        for e in self.G.edges():
            dist = np.linalg.norm([pos[e[0]][0]-pos[e[1]][0],pos[e[0]][1]-pos[e[1]][1]])
            reverse_e = (e[1],e[0])
            self.weights[e] = dist
            self.transmission[e] = 10**(-0.2*dist/10)
            
            self.weights[reverse_e] = dist
            self.transmission[reverse_e] = 10**(-0.2*dist/10)

        nx.set_edge_attributes(self.G, values = self.weights, name = 'weight')
        nx.set_edge_attributes(self.G, values = self.transmission, name = 'trans')

        plt.figure(figsize=(4,3))
        nx.draw(self.G,pos,with_labels=True)
        plt.show()
        if not nx.is_connected(self.G):
            return False
        else:
            user_pairs = []
            for node1 in self.G.nodes:
                if node1 not in end_nodes:
                    if node1 not in self.repeater_places:
                        self.repeater_places.append(node1)
                for node2 in self.G.nodes:
                    if node1!=node2:
                        e = (node1,node2)
                        if e not in self.G.edges:
                            shortest_path = self.get_paths_between_user_pairs(e,1)
                            if (node1,node2) not in user_pairs and (node2,node1) not in user_pairs:
                                user_pairs.append((node1,node2))
                            path_edges = []
                            node_indx = 0
                            shortest_path = shortest_path[0]
                            for node_indx in range(len(shortest_path)-1):                
                                path_edges.append((shortest_path[node_indx],shortest_path[node_indx+1]))
                                node_indx+=1
                            dist = 0
                            for edge in path_edges:
                                dist = dist+self.weights[edge]
                            self.weights[e] = dist
                            self.weights[(e[1],e[0])] = dist
                            self.G.add_edge(node1,node2,weight=dist)
                            self.transmission[e] = 10**(-0.2*dist/10)
                            p = 10**(-0.2*dist/10)
                            self.set_of_edges.append(e)
                            self.set_of_edges.append((e[1],e[0]))
                            self.transmission[(e[1],e[0])] = 10**(-0.2*dist/10)

            if len(user_pairs)>=self.number_of_user_pairs:
                while(len(self.user_pairs)<self.number_of_user_pairs):
                    random_user_pair = user_pairs[random.randint(0,len(user_pairs)-1)]
                    if random_user_pair not in self.user_pairs and (random_user_pair[1],random_user_pair[0]) not in self.user_pairs:
                        self.user_pairs.append((random_user_pair[0],random_user_pair[1]))

                user_pairs_nodes = []
                self.repeater_places = []
                for pair in self.user_pairs:
                    user_pairs_nodes.append(pair[0])
                    user_pairs_nodes.append(pair[1])
                for node1 in self.G.nodes:
                    if node1 not in  self.repeater_places and node1 not in user_pairs_nodes:
                        self.repeater_places.append(node1)
                self.N =len(self.G.nodes)
                self.node_list = np.arange(self.N)
                return True
            else:
                return False
    def set_paths_required_decoherence_time_new_version(self):

        for path_id,path_edges in self.set_of_paths.items():
            sum_length = 0
            repeater_loc = []
            for edge in path_edges:
                d2 =self.weights[edge]
                sum_length = sum_length+d2
                repeater_loc.append(sum_length)
            repeater_loc = np.array(repeater_loc)
            Nrepeater = len(repeater_loc)-2 # number of repeater nodes
            events = np.zeros((len(repeater_loc),3)) # 3 events at each node
            if elementary_link_sync:
                events[1:,0] = repeater_loc[1:]-repeater_loc[:-1]
                events[1,2] = events[2,0] + (repeater_loc[2]-repeater_loc[1])
                for i in range(2,Nrepeater+1):
                    events[i,1]= np.max([events[i,0],events[i+1,0]+repeater_loc[i+1]-repeater_loc[i]])
                    events[i,2]= events[i-1,2] + repeater_loc[i]-repeater_loc[i-1]
                    events[i,2]= np.max(events[i,1:])
            else:
                events[:,0] = repeater_loc # receiving quantum signal
                events[1,2] = events[1,0] + 2* (repeater_loc[2]-repeater_loc[1])
                for i in range(2,Nrepeater+1):
                    # for node i
                    # receiving classical acknowledgement from next node and possible swap result from previous node
                    events[i,1]= repeater_loc[i]+ 2*np.max([repeater_loc[i]-repeater_loc[i-1],repeater_loc[i+1]-repeater_loc[i]])
                    events[i,2]= events[i-1,2] + repeater_loc[i]-repeater_loc[i-1] # receiving swap result
                    events[i,2]= np.max(events[i,1:]) # whichever event takes longer determines the swapping schedule at node i
            #     return events
            total_time = events[Nrepeater,2] + repeater_loc[Nrepeater]
            node_time = np.max(events[1:-1,2]-events[1:-1,0])
            
            self.each_path_repeater_required_decoherence_time[path_id] = node_time
            self.each_path_end_node_required_decoherence_time[path_id] = total_time
            
    def set_paths_required_decoherence_time(self):

        for path_id,path_edges in self.set_of_paths.items():
            
            end_memory_waiting_times = []
            for edge in path_edges:
                d2 =self.weights[edge]
                time_for_d2 = (1.44*d2)/ 299792
                end_memory_waiting_times.append(3*time_for_d2*1000000)# we get the waiting time per micro seconds
                
            first_hop = path_edges[0]
            first_hop_length =self.weights[first_hop]
            first_hop_waiting_time = (1.44*first_hop_length)/ 299792
            first_hop_waiting_time = 2*first_hop_waiting_time*1000000
            if len(path_edges)<=1:
                waiting_times = [0]
            else:
                waiting_times = []
                for edge in path_edges:
                #     print("edge ",edge)#2*(1.44*link_length)/ 299792
                    d2 =self.weights[edge]
                    time_for_d2 = (1.44*d2)/ 299792
                    waiting_times.append(2*time_for_d2*1000000)# we get the waiting time per micro seconds
            self.each_path_repeater_required_decoherence_time[path_id] = max(waiting_times)
            self.each_path_end_node_required_decoherence_time[path_id] = sum(end_memory_waiting_times)

       

    
    
    
    


# In[2]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




