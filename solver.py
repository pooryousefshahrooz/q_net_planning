#!/usr/bin/env python
# coding: utf-8

import numpy as np
import random
import csv
import os
import sys
from docplex.mp.progress import *
from docplex.mp.progress import SolutionRecorder
import docplex.mp.model as cpx
import networkx as nx
import time
import matplotlib.pyplot as plt
from config import get_config
from absl import flags
FLAGS = flags.FLAGS

class Solver:
    def __init__(self):
        pass
    
    def pyomo_linear_link_based_utility_maximization(self,network):
        ### minimum capacity
        Nmax = 2
        q = 0.5
        F = 0.95
        D0 = self.D
        D = D0*np.ones(network.N)
        # w = 4
        W_list = np.arange(1,D0+1) # maximum memory allowed for end users
        s_list = [0] # source
        t_list = [network.N-1] # destination

        self.s_list = [] # source
        self.t_list = [] # destination

        for user_pair in network.user_pairs:
            self.s_list.append(user_pair[0])
            self.t_list.append(user_pair[1])
        self.C = len(self.s_list)
        R_list = np.array(list(set(network.node_list)-set(np.concatenate((self.s_list,self.t_list)))))

        #Create a simple model
        model = pe.ConcreteModel()

        model.R = pe.Set(initialize=R_list)
        model.N = pe.Set(initialize=range(self.N)) 
        model.C = pe.Set(initialize=range(self.C))
        model.W = pe.Set(initialize=W_list)
        model.x = pe.Var(model.C,model.W,model.N,model.N, domain=pe.Binary)#,initialize=0)
        model.y = pe.Var(model.R, domain=pe.Binary)#,initialize=1)
        model.minCapacity = pe.Var(model.C,domain=pe.NonNegativeReals)
        model.whichmemory = pe.Var(model.C,model.W)#,domain=pe.NonNegativeReals)
        model.hops = pe.Var(model.C)#,domain=pe.NonNegativeReals)

        model.constraints = pe.ConstraintList()

        for n1 in R_list:
            n_list = np.array(list(set(R_list)-{n1}))
            sum1_const = 0 
            for user_pair in range(self.C):
                s = self.s_list[user_pair]
                t = self.t_list[user_pair]
                sum1_const += sum( i*model.x[user_pair,i,n1,t] for i in model.W) 
                for n2 in n_list:
                    sum1_const += sum( i*model.x[user_pair,i,n1,n2] for i in model.W) 
                    model.constraints.add( sum(model.x[user_pair,i,n1,n2] for i in model.W)<=1 )

            model.constraints.add( sum1_const <= D[n1]*model.y[n1] )

        path_length = 0
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]
            path_s = 0
            path_t = 0
            number = 0
            for memory in model.W:

                for n1 in range(self.N):
                    model.constraints.add(  model.x[user_pair,memory,n1,n1] == 0 ) 
                    model.constraints.add(  model.x[user_pair,memory,n1,s] == 0 ) 
                    model.constraints.add(  model.x[user_pair,memory,t,n1] == 0 ) 

                number +=  model.x[user_pair,memory,s,t]+ sum(model.x[user_pair,memory,s,i]+ model.x[user_pair,memory,i,t] for i in model.R)
                for n1 in R_list:
                    n_list = np.array(list(set(R_list)-{n1}))
                    number += sum(model.x[user_pair,memory,n1,i] for i in n_list)

                    if nx.has_path(self.G, source=s, target=n1):
                        l_uv = nx.shortest_path_length(self.G, source=s, target=n1, weight="weight")
                        path_length += l_uv*  model.x[user_pair,memory,s,n1]
                        model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,memory,s,n1] )

                    if nx.has_path(self.G, source=n1, target=t):
                        l_uv = nx.shortest_path_length(self.G, source=n1, target=t, weight="weight")
                        path_length += l_uv* model.x[user_pair,memory,n1,t]
                        model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,memory,n1,t] )

                    path = model.x[user_pair,memory,n1,t] - model.x[user_pair,memory,s,n1]
                    for n2 in R_list:
                        if n2 != n1:
                            path += model.x[user_pair,memory,n1,n2] - model.x[user_pair,memory,n2,n1]
                            if nx.has_path(self.G, source=n1, target=n2):
                                l_uv = nx.shortest_path_length(self.G, source=n1, target=n2, weight="weight")
                                path_length += l_uv* model.x[user_pair,memory,n1,n2]
                                model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,memory,n1,n2] )

                    model.constraints.add( path == 0 ) 
                l_uv = nx.shortest_path_length(self.G, source=s, target=t, weight="weight")
                model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,memory,s,t] )

                path_s += model.x[user_pair,memory,s,t] + sum(model.x[user_pair,memory,s,i] for i in model.R) 
                path_t += model.x[user_pair,memory,s,t] + sum(model.x[user_pair,memory,i,t] for i in model.R)
                model.constraints.add( model.whichmemory[user_pair,memory]== model.x[user_pair,memory,s,t] + sum(model.x[user_pair,memory,s,i] for i in model.R) )
        #         model.constraints.add( expr= model.x[user_pair,memory,s,t] + sum(model.x[user_pair,memory,i,t] for i in model.R) == 1)
        #     #     model.constraints.add( model.fidelity[user_pair] == 3*((4*F-1)/3)**number-1)
            model.constraints.add( model.hops[user_pair] == number )
            model.constraints.add( path_s== 1)
            model.constraints.add( path_t== 1)

        # model.constraints.add(expr= sum(model.y[i] for i in model.R) <= Nmax)
        # model.objective = pe.Objective(rule=ObjRule, sense=pe.maximize)
        α = 0.02*np.log(10)
        # model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+ np.log(q)*model.hops[i] +pe.log(model.fidelity[i]) for i in model.C), sense=pe.maximize)
        # model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+ pe.log(model.memory[i])+ np.log(q)*model.hops[i] +pe.log(3*((4*F-1)/3)**model.hops[i]-1) for i in model.C), sense=pe.maximize)
        model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+
                                                 sum(model.whichmemory[i,m]*np.log(m) for m in model.W) +
                                                 np.log(q)*model.hops[i] for i in model.C), sense=pe.maximize)
        # λ = 1e-2
        # model.objective = pe.Objective(expr= (sum(pe.log(model.minCapacity[i]*fidelity[i]) for i in model.C) - λ* path_length), sense=pe.maximize)
        # model.objective = pe.Objective(expr= model.minCapacity[0]*fidelity[0] , sense=pe.maximize)

        # opt = pe.SolverFactory("bonmin",tee=True)
#         opt = pe.SolverFactory("cbc",tee=True,executable="/work/spooryousefd_umass_edu/bonmin/bin/cbc")
        opt = pe.SolverFactory("cplex",tee=True,executable="/modules/apps/cplex/2210/cplex/bin/x86-64_linux/cplex")
        
        results = opt.solve(model)
        results.write()

        # model.objective.display()
        # model.display()
        # model.pprint()

        # x_opt = np.zeros((C,w,N,N))
        x_opt = np.zeros((self.C,len(W_list),self.N,self.N))
        y_opt = np.zeros(self.N)
        for i in R_list:
            y_opt[i] = model.y[i].value

        for i in range(self.N):
            for j in range(self.N):
                for user_pair in range(self.C):
                    for i_m, memory in enumerate(W_list):
                        x_opt[user_pair,i_m,i,j] =  model.x[user_pair,memory,i,j].value
                        
        each_user_pair_allocated_W = {}
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
#             plt.plot(pos[s][0],pos[s][1],"s", color = colors[user_pair])
#             plt.text(pos[s][0],pos[s][1],"%d" % s)
            t = self.t_list[user_pair]
#             plt.plot(pos[t][0],pos[t][1],"s", color = colors[user_pair])
#             plt.text(pos[t][0],pos[t][1],"%d" % t)
            for i_m, memory in enumerate(W_list):
                edges = np.argwhere(x_opt[user_pair,i_m,:,:]>0.5)
                if len(edges)>0:
                    print("optimal memory for ", (s,t)," :", memory)
                    each_user_pair_allocated_W[(s,t)] =memory
        self.path_extraction(x_opt,each_user_pair_allocated_W)
        utility_value = pe.value(model.objective)
        return utility_value
    
    
    
    def path_extraction(self,x_opt,each_user_pair_allocated_W):
        used_edges = []
        self.solution_path_counter = 0
        self.solution_used_nodes_by_paths = []
        self.solution_each_user_pair_used_paths = {}
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]
            edges = np.argwhere(x_opt[user_pair,:,:]>0.5)
            path_output = []
            for e in edges:
                if self.G.has_edge(e[0],e[1]):
                    path_output.append(list(e))
                else:
                    path = nx.shortest_path(self.G, source=e[0], target=e[1], weight="weight")
                    path_output.append(path)
                    for i in range(len(path)-1):
                        e1 = path[i]
                        e2 = path[i+1]

            user_pair_src_dst = (self.s_list[user_pair],self.t_list[user_pair])
            
            path_edges = []
            for edge in path_output:
                path_edges.append((edge[0],edge[1]))
            self.solution_each_user_pair_used_paths[user_pair_src_dst,self.path_counter] = [path_edges]
            
            
            W = int(each_user_pair_allocated_W[user_pair_src_dst])
            E_p =self.compute_e2e_rate(self.path_counter,path_edges,W)
            e2e_F = self.compute_e2e_fidleity(path_edges)
            self.each_path_width[self.path_counter] = W
            self.each_path_e2e_F_value[self.path_counter] = e2e_F
            self.each_path_e_value[self.path_counter] =E_p
            self.path_counter+=1
            self.solution_path_counter+=1
            
            for edge in path_output:
                if edge[0] not in self.solution_used_nodes_by_paths and edge[0] not in user_pair_src_dst:
                    self.solution_used_nodes_by_paths.append(edge[0])
                if edge[1] not in self.solution_used_nodes_by_paths and edge[1] not in user_pair_src_dst:
                    self.solution_used_nodes_by_paths.append(edge[1])
        self.solution_used_nodes_by_paths_string = ""
        for node in self.solution_used_nodes_by_paths:
            if self.solution_used_nodes_by_paths_string:
                self.solution_used_nodes_by_paths_string = self.solution_used_nodes_by_paths_string+":"+str(node)
            else:
                self.solution_used_nodes_by_paths_string = str(node)
        
    def link_based_non_linear_formulation(self):
        D = self.D*np.ones(self.N)


        self.s_list = [] # source
        self.t_list = [] # destination

        for user_pair in self.user_pairs:
            self.s_list.append(user_pair[0])
            self.t_list.append(user_pair[1])

        self.C = len(self.s_list)
        self.R_list = np.array(list(set(self.node_list)-set(np.concatenate((self.s_list,self.t_list)))))
        #Create a simple model
        model = pe.ConcreteModel()

        model.R = pe.Set(initialize=self.R_list)
        model.N = pe.Set(initialize=range(self.N))
        model.C = pe.Set(initialize=range(self.C))
        model.x = pe.Var(model.C,model.N,model.N, domain=pe.Binary)#,initialize=0)
        model.y = pe.Var(model.R, domain=pe.Binary)#,initialize=1)
        model.minCapacity = pe.Var(model.C,domain=pe.NonNegativeReals)
        model.memory = pe.Var(model.C, bounds=(1, self.end_user_memory),domain=pe.PositiveIntegers)
        # model.fidelity = pe.Var(model.C)#,domain=pe.NonNegativeReals)
        model.hops = pe.Var(model.C)#,domain=pe.NonNegativeReals)

        model.constraints = pe.ConstraintList()

        for n1 in self.R_list:
            n_list = np.array(list(set(self.R_list)-{n1}))
            sum1_const = 0 
            for user_pair in range(self.C):
                s = self.s_list[user_pair]
                t = self.t_list[user_pair]
                sum1_const += model.memory[user_pair]*model.x[user_pair,n1,t] 
                for n2 in n_list:
                    sum1_const += model.memory[user_pair]*model.x[user_pair,n1,n2]
            model.constraints.add( sum1_const <= D[n1]*model.y[n1] )

        path_length = 0
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]

        #     model.constraints.add( model.memory[user_pair] <= w )

            for n1 in range(self.N):
                model.constraints.add(  model.x[user_pair,n1,n1] == 0 ) 
                model.constraints.add(  model.x[user_pair,n1,s] == 0 ) 
                model.constraints.add(  model.x[user_pair,t,n1] == 0 ) 

            number =  model.x[user_pair,s,t]+ sum(model.x[user_pair,s,i]+
                                                  model.x[user_pair,i,t] for i in model.R)
            for n1 in self.R_list:
                n_list = np.array(list(set(self.R_list)-{n1}))
                number += sum(model.x[user_pair,n1,i] for i in n_list)

                if nx.has_path(self.G, source=s, target=n1):
                    l_uv = nx.shortest_path_length(self.G, source=s, target=n1, weight="weight")
                    path_length += l_uv*  model.x[user_pair,s,n1]
                    model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,s,n1] )

                if nx.has_path(self.G, source=n1, target=t):
                    l_uv = nx.shortest_path_length(self.G, source=n1, target=t, weight="weight")
                    path_length += l_uv* model.x[user_pair,n1,t]
                    model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,n1,t] )

                path = model.x[user_pair,n1,t] - model.x[user_pair,s,n1]
                for n2 in self.R_list:
                    if n2 != n1:
                        path += model.x[user_pair,n1,n2] - model.x[user_pair,n2,n1]
                        if nx.has_path(self.G, source=n1, target=n2):
                            l_uv = nx.shortest_path_length(self.G, source=n1, target=n2, weight="weight")
                            path_length += l_uv* model.x[user_pair,n1,n2]
                            model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,n1,n2] )

                model.constraints.add( path == 0 ) 
            l_uv = nx.shortest_path_length(self.G, source=s, target=t, weight="weight")
            model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,s,t] )

            model.constraints.add( expr= model.x[user_pair,s,t] + 
                                  sum(model.x[user_pair,s,i] for i in model.R) == 1)
            model.constraints.add( expr= model.x[user_pair,s,t] + 
                                  sum(model.x[user_pair,i,t] for i in model.R) == 1)
            model.constraints.add( model.hops[user_pair] == number )

        α = 0.02*np.log(10)

        model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+
                                                 pe.log(model.memory[i])+ 
                                                 np.log(self.q)*model.hops[i] 
                                                 for i in model.C), 
                                       sense=pe.maximize)

        # opt = pe.SolverFactory("bonmin",tee=True)
#         opt = pe.SolverFactory("bonmin", executable="/work/spooryousefd_umass_edu/Bonmin-1.8.8/bin/bonmin",
#                                tee=True)
        if self.link_based_solver =="Baron":
            opt = pe.SolverFactory("baron", executable="/path/to/baron-lin64/baron",
                                   tee=True)
        else:
            print("setting Bonmin as solver...")
            opt = pe.SolverFactory("bonmin", executable="/path/to/Bonmin-1.8.8/bin/bonmin",
                               tee=True)
        print("start solving ....")
        results = opt.solve(model)
        results.write()

        # model.objective.display()
        # model.display()
        # model.pprint()

        x_opt = np.zeros((self.C,self.N,self.N))
        y_opt = np.zeros(self.N)
        for i in self.R_list:
            y_opt[i] = model.y[i].value

        for i in range(self.N):
            for j in range(self.N):
                for user_pair in range(self.C):
                    x_opt[user_pair,i,j] =  model.x[user_pair,i,j].value


        each_user_pair_allocated_W = {}
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]
#             print("optimal memory for ", (s,t)," :", model.memory[user_pair].value)#, pe.value(rate[user_pair]) )
            each_user_pair_allocated_W[(s,t)] = model.memory[user_pair].value
        #model.objective.display()
        utility_value = pe.value(model.objective)
#         utility2 = model.objective.value()
        used_repeaters = set([])
        y_opt = np.zeros(self.N)
        for n1 in self.R_list:
            n_list = np.array(list(set(self.R_list)-{n1}))
            sum1_const = 0 
            for user_pair in range(self.C):
                s = self.s_list[user_pair]
                t = self.t_list[user_pair]
                sum1_const += x_opt[user_pair,n1,t] 
                for n2 in n_list:
                    sum1_const += x_opt[user_pair,n1,n2]
            if sum1_const>0.5:
                y_opt[n1] = 1
                used_repeaters.add(n1)
                
                
                
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]

        self.solution_repeater_counter = len(used_repeaters)
        
#         self.plot_output(x_opt,y_opt)       
        self.path_extraction(x_opt,each_user_pair_allocated_W)
#         self.plot_colored_graph(x_opt,y_opt,used_repeaters,used_edges)

        return utility_value    
        
    def path_based_utility_maximization(self,dynamic_system_flag,network,used_repeaters):
        """we solve the utility maximization using R repeaters and M memories"""

        opt_model = cpx.Model(name="maximizing_utility")
        x_p_vars  = {(p): opt_model.integer_var(lb=0, ub= 1,
                                  name="x_p_{0}".format(p))  
                   for user_pair in network.user_pairs
                for p in network.each_user_pair_paths[user_pair]}

        r_u_vars  = {(r): opt_model.integer_var(lb=0, ub= 1,
                                  name="r_{0}".format(r))  
                   for r in network.repeater_places}
        if dynamic_system_flag:
            for u in used_repeaters:
#                 print("used repeater ",u)
                opt_model.add_constraint(
                    r_u_vars[u]
                     >= 1.0, ctname="repeater_constraint1_{0}".format(u))

            for u in network.repeater_places:
                if u not in used_repeaters:
#                     print("not used repeater ",u)
                    opt_model.add_constraint(
                        r_u_vars[u]
                         <= 0.0, ctname="repeater_constraint2_{0}".format(u))
        # Repeater memory constraint
        for u in network.repeater_places:
            opt_model.add_constraint(
                opt_model.sum(x_p_vars[p]*network.each_path_width[p]
                for user_pair in network.user_pairs
                for p in network.each_user_pair_paths[user_pair]
                if network.check_path_uses_repeater(u,p)
                )
                 <= network.D *r_u_vars[u], ctname="repeater_memory_constraint_{0}".format(u))

        # number of memory and repeaters constraint
        opt_model.add_constraint(
                opt_model.sum(r_u_vars[node]
                for node in network.repeater_places
                )<= network.R, ctname="number_of_repeaters_constraint")

        

        for user_pair in network.user_pairs:
            opt_model.add_constraint(
                opt_model.sum(x_p_vars[p]* network.each_path_width[p]
                    for p in network.each_user_pair_paths[user_pair]
                )
                     <= network.end_user_memory, ctname="end_node_capacity_constraint_{0}".format(user_pair))

        # each user pair can use at most network.K number of paths
        for user_pair in network.user_pairs:
            opt_model.add_constraint(
                opt_model.sum(x_p_vars[p]
                    for p in network.each_user_pair_paths[user_pair]
                )
                     <= network.K, ctname="using_multiple_paths_constraint_{0}".format(user_pair))
            
        all_end_node_available_paths = set([])
        all_repeater_available_paths = set([])
        if network.checking_repeater_memory_decoherence_flag:
            for user_pair in network.user_pairs:
                for p in network.each_user_pair_paths[user_pair]:
                    all_repeater_available_paths.add(network.each_path_repeater_required_decoherence_time[p])
                    opt_model.add_constraint(
                        x_p_vars[p] * network.each_path_repeater_required_decoherence_time[p]
                             <= network.repeater_memory_decoherence_time, ctname="path_repeater_decoherence_time_constraint_{0}_{1}".format(user_pair,p))

        # using the paths that have the waiting time less than the end node decoherence time
        if network.checking_end_node_memory_decoherence_flag:
            for user_pair in network.user_pairs:
                for p in network.each_user_pair_paths[user_pair]:
                    all_end_node_available_paths.add(network.each_path_end_node_required_decoherence_time[p])
                    opt_model.add_constraint(
                        x_p_vars[p] * network.each_path_end_node_required_decoherence_time[p]
                             <= network.end_node_memory_decoherence_time, ctname="path_end_node_decoherence_time_constraint_{0}_{1}".format(user_pair,p))
                
                   
        for user_pair in network.user_pairs:
            opt_model.add_constraint(
                opt_model.sum(x_p_vars[p]
                    for p in network.each_user_pair_paths[user_pair]
                )
                     >= 1, ctname="using_at_least_one_path_constraint_{0}".format(user_pair))
 
        if network.weighted_sum_rate_objective:
            objective = opt_model.sum(network.compute_utility(p)*x_p_vars[p]*network.each_user_pair_weight[pair]
                              for pair in network.user_pairs for p in network.each_user_pair_paths[pair] 
                              )
            
        else:
            
            objective = opt_model.sum(network.compute_utility(p)*x_p_vars[p]
                              for pair in network.user_pairs for p in network.each_user_pair_paths[pair] 
                              )
        
        objective_value = network.utility_default_value
        # for maximization
        opt_model.maximize(objective)
        print("going to solve the optimization")
        try:
            opt_model.solve()
        except:
            objective_value = network.utility_default_value

        set_of_used_repeaters = []
        set_of_used_paths = []


        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
        except:
            objective_value=network.utility_default_value
#             print(ValueError)
        nodes_string = ""
        
        network.solution_repeater_counter=0 
        used_nodes = []
        used_repeaters_list = []
        network.solution_each_user_pair_used_paths = {}
        network.solution_each_user_pair_paths_edges_str  ={}
        if objective_value!=network.utility_default_value:
            
            for user_pair in network.user_pairs:
                path_edges_str = ""
                for p in network.each_user_pair_paths[user_pair]:
                    this_path_string = ""
                    if x_p_vars[p].solution_value>0.5:
                        try:
                            network.solution_each_user_pair_used_paths[user_pair,p].append(network.set_of_paths[p])
                        except:
                            network.solution_each_user_pair_used_paths[user_pair,p] = [network.set_of_paths[p]]
                        for edge in network.set_of_paths[p]:
                            if this_path_string:
                                this_path_string = this_path_string+"-"+str(edge[0])+":"+str(edge[1])
                            else:
                                this_path_string = str(edge[0])+":"+str(edge[1])
                                
                            
                    if this_path_string:
                        if path_edges_str:
                            path_edges_str = path_edges_str+"|"+this_path_string
                        else:
                            path_edges_str = this_path_string
                network.solution_each_user_pair_paths_edges_str[user_pair] = path_edges_str
            
            
            for r in network.repeater_places:
                if r_u_vars[r].solution_value>0.5:
                    network.solution_repeater_counter+=1
                    used_repeaters_list.append(r)


        network.solution_path_counter = 0
        network.solution_used_paths_by_user = {}
        network.solution_used_nodes_by_paths = []
        for user_pair in network.user_pairs:
            user_pair_path_counter = 0
            for p in network.each_user_pair_paths[user_pair]:
                if x_p_vars[p].solution_value>0.5:
                    user_pair_path_counter+=1
                    network.solution_path_counter+=1
                    network.solution_used_paths_by_user[p]=network.set_of_paths[p]
                    for edge in network.set_of_paths[p]:
                        if edge[0] not in network.solution_used_nodes_by_paths and edge[0] not in user_pair:
                            network.solution_used_nodes_by_paths.append(edge[0])
                        if edge[1] not in network.solution_used_nodes_by_paths and edge[1] not in user_pair:
                            network.solution_used_nodes_by_paths.append(edge[1])


        return objective_value,used_repeaters_list
    
    
    def linear_link_based_in_CPLEX(self,network):
        
        
        used_repeaters = set([])
        D = network.D*np.ones(len(network.G.nodes)+1)
        w0 = network.end_user_memory
        
        
        W_list = np.arange(1,w0+1) # maximum memory allowed for end users

        self.s_list = [] # source
        self.t_list = [] # destination
        print("W_list ",W_list)
        for user_pair in network.user_pairs:
            self.s_list.append(user_pair[0])
            self.t_list.append(user_pair[1])
        self.C = len(self.s_list)
        R_list = np.array(list(set(network.node_list)-set(np.concatenate((self.s_list,self.t_list)))))
        
        
        #Create a simple model
        opt_model = cpx.Model(name="linear_link_based_maximizing_utility")
        R_list = network.repeater_places
        list_R = R_list
#         list_N = [i for i in range(network.N)]
        list_N = []
        for node in network.G.nodes:
            list_N.append(node)
        list_C = [i for i in range(self.C)]
        list_W = W_list
        
        x  = {(i,j,k,w): opt_model.integer_var(lb=0, ub= 1,
                                  name="x_{0}_{1}_{2}_{3}".format(i,j,k,w))  
                   for i in list_C for j in list_W for k in list_N for w in list_N}
        y  = {(i): opt_model.integer_var(lb=0, ub= 1,
                                  name="y_{0}".format(i))  
                   for i in list_R}
        minCapacity = {(i): opt_model.continuous_var(lb=0,
                                  name="minCapacity_{0}".format(i))  
                   for i in list_C}
        whichmemory = {(i,w): opt_model.integer_var(lb=0, ub= min(network.D,network.end_user_memory),
                                  name="whichmemory_{0}_{1}".format(i,w))  
                   for i in list_C for w in list_W}
        
        hops = {(i): opt_model.integer_var(lb=0, ub= len(network.G.nodes),
                                  name="hops_{0}".format(i))  
                   for i in list_C}


        for n1 in R_list:
            n_list = np.array(list(set(R_list)-{n1}))

            sum1_const = 0 
            for user_pair in range(self.C):
                s = self.s_list[user_pair]
                t = self.t_list[user_pair]
                sum1_const += sum( i*x[user_pair,i,n1,t] for i in list_W) 
                for n2 in n_list:
                    sum1_const += sum( i*x[user_pair,i,n1,n2] for i in list_W)
                    opt_model.add_constraint(opt_model.sum(x[user_pair,i,n1,n2] for i in list_W)<= 1)

            opt_model.add_constraint(sum1_const<= D[n1]*y[n1])
            
            
        # number of memory and repeaters constraint
        opt_model.add_constraint(
                opt_model.sum(y[node]
                for node in network.repeater_places
                )<= network.R, ctname="number_of_repeaters_constraint")   
            
        path_length = 0
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]
            path_s = 0
            path_t = 0
            number = 0
            for memory in list_W:
                for n1 in range(network.N):
                    opt_model.add_constraint(x[user_pair,memory,n1,n1] == 0 ) 
#                     opt_model.add_constraint(x[user_pair,memory,n1,s] == 0 ) 
#                     opt_model.add_constraint(x[user_pair,memory,t,n1] == 0 ) 

                number +=  x[user_pair,memory,s,t]+ sum(x[user_pair,memory,s,i]+x[user_pair,memory,i,t] for i in list_R)
                for n1 in R_list:
                    n_list = np.array(list(set(R_list)-{n1}))
                    number += sum(x[user_pair,memory,n1,i] for i in n_list)

                    if nx.has_path(network.G, source=s, target=n1):
                        l_uv = nx.shortest_path_length(network.G, source=s, target=n1, weight="weight")
                        path_length += l_uv*  x[user_pair,memory,s,n1]
                        #model.constraints.add(model.minCapacity[user_pair] >=  l_uv*model.x[user_pair,memory,s,n1] )
                        opt_model.add_constraint(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,s,n1] )

                    if nx.has_path(network.G, source=n1, target=t):
                        l_uv = nx.shortest_path_length(network.G, source=n1, target=t, weight="weight")
                        path_length += l_uv* x[user_pair,memory,n1,t]
                        #model.constraints.add(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,n1,t] )
                        opt_model.add_constraint(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,n1,t] )
                    path = x[user_pair,memory,n1,t] - x[user_pair,memory,s,n1]
                    for n2 in R_list:
                        if n2 != n1:
                            path += x[user_pair,memory,n1,n2] - x[user_pair,memory,n2,n1]
                            if nx.has_path(network.G, source=n1, target=n2):
                                l_uv = nx.shortest_path_length(network.G, source=n1, target=n2, weight="weight")
                                path_length += l_uv* x[user_pair,memory,n1,n2]
                                #model.constraints.add(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,n1,n2] )
                                opt_model.add_constraint(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,n1,n2])
                    #model.constraints.add( path == 0 ) 
                    opt_model.add_constraint(path == 0)
                l_uv = nx.shortest_path_length(network.G, source=s, target=t, weight="weight")
                #model.constraints.add(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,s,t] )
                opt_model.add_constraint(minCapacity[user_pair] >=  l_uv*x[user_pair,memory,s,t])
                path_s += x[user_pair,memory,s,t] + sum(x[user_pair,memory,s,i] for i in list_R) 
                path_t += x[user_pair,memory,s,t] + sum(x[user_pair,memory,i,t] for i in list_R)
                #model.constraints.add( whichmemory[user_pair,memory]== x[user_pair,memory,s,t] + sum(x[user_pair,memory,s,i] for i in list_R) )
                opt_model.add_constraint(whichmemory[user_pair,memory]== x[user_pair,memory,s,t] + sum(x[user_pair,memory,s,i] for i in list_R))

            opt_model.add_constraint(hops[user_pair] == number-1)
    #         model.constraints.add( path_s== 1)
            opt_model.add_constraint(path_s== 1)
    #         model.constraints.add( path_t== 1)
            opt_model.add_constraint(path_t== 1)
        
        if network.checking_memory_decoherence_flag:
            for user_pair in network.user_pairs:
                opt_model.add_constraint(
                    (1.44*minCapacity[user_pair])/ 299792
                         <= network.repeater_memory_decoherence_time, 
                    ctname="repeater_life_time_constraint_{0}_{1}".format(user_pair,p))


                
        
        
        α = 0.02*np.log(10)
        # model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+ np.log(q)*model.hops[i] +pe.log(model.fidelity[i]) for i in model.C), sense=pe.maximize)
        # model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+ pe.log(model.memory[i])+ np.log(q)*model.hops[i] +pe.log(3*((4*F-1)/3)**model.hops[i]-1) for i in model.C), sense=pe.maximize)
#         model.objective = pe.Objective(expr= sum(-α*model.minCapacity[i]+ sum(model.whichmemory[i,m]*np.log(m) for m in model.W) +
#                                                  np.log(q)*model.hops[i] for i in model.C), sense=pe.maximize)
        objective = opt_model.sum(-α*minCapacity[i]+ sum(whichmemory[i,m]*np.log(m) for m in list_W) + np.log(network.q)*hops[i] for i in list_C 
                                  )
        # λ = 1e-2
        # model.objective = pe.Objective(expr= (sum(pe.log(model.minCapacity[i]*fidelity[i]) for i in model.C) - λ* path_length), sense=pe.maximize)
        # model.objective = pe.Objective(expr= model.minCapacity[0]*fidelity[0] , sense=pe.maximize)

        # opt = pe.SolverFactory("baron",tee=True)
#         opt = pe.SolverFactory("cplex",tee=True,executable="/modules/apps/cplex/2210/cplex/bin/x86-64_linux/cplex")
#         results = opt.solve(model)
#         results.write()
        objective_value = network.utility_default_value
        # for maximization
        opt_model.maximize(objective)
    #     opt_model.solve()
#         opt_model.print_information()
        try:
            start_time = time.time()
            print("start solving.....")
            opt_model.solve()
            
            
        except ValueError:
            print("ValueError",ValueError)
            objective_value = network.utility_default_value
        try:
            if opt_model.solution:
                objective_value =opt_model.solution.get_objective_value()
                current_time = time.time()
                print("solved in %s seconds and utility is %s "%(current_time-start_time,objective_value))
        except ValueError:
            print(ValueError)
            objective_value=network.utility_default_value
#         opt_model.objective.display()
#         opt_model.display()
#         opt_model.pprint()
        opt_model.print_information()
        print('docplex.mp.solution',opt_model.solution)
        x_opt = np.zeros((self.C,len(W_list),network.N,network.N))
        y_opt = np.zeros(network.N)
        for i in R_list:
            y_opt[i] = y[i].solution_value

        for i in range(network.N):
            for j in range(network.N):
                for user_pair in range(self.C):
                    for i_m, memory in enumerate(W_list):
                        x_opt[user_pair,i_m,i,j] =  x[user_pair,memory,i,j].solution_value
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
            t = self.t_list[user_pair]

        

        y_opt = np.zeros(network.N)
        for n1 in R_list:
            n_list = np.array(list(set(R_list)-{n1}))
            sum1_const = 0 
            for user_pair in range(self.C):
                s = self.s_list[user_pair]
                t = self.t_list[user_pair]
                for memory in range(len(W_list)):
                    sum1_const += x_opt[user_pair,memory,n1,t] 
                    for n2 in n_list:
                        sum1_const += x_opt[user_pair,memory,n1,n2]
            if sum1_const>0.1:
                y_opt[n1] = 1
                used_repeaters.add(n1)

        network.solution_repeater_counter = len(list(used_repeaters))
        self.plot_output_lp(x_opt,y_opt,W_list,network) 
        return objective_value
        
        
        
    def plot_output_lp(self,x_opt,y_opt,W_list,network):
        
        used_edges = []
        network.solution_path_counter = 0
        network.solution_used_nodes_by_paths = []
        network.solution_each_user_pair_used_paths = {}
   
        colors = ["violet","orange","gray"]
        for user_pair in range(self.C):
            s = self.s_list[user_pair]
#             plt.plot(network.pos[s][0],network.pos[s][1],"s", color = colors[user_pair])
#             plt.text(network.pos[s][0],network.pos[s][1],"%d" % s)
            t = self.t_list[user_pair]
#             plt.plot(network.pos[t][0],network.pos[t][1],"s", color = colors[user_pair])
#             plt.text(network.pos[t][0],network.pos[t][1],"%d" % t)
            for i_m, memory in enumerate(W_list):

                edges = np.argwhere(x_opt[user_pair,i_m,:,:]>0.5)
                if len(edges)>0:
                    print("optimal memory for ", (s,t)," :", memory)
                    path_output = []
                    this_user_path_W = memory
                    for e in edges:
                        if network.G.has_edge(e[0],e[1]):
#                             plt.plot([ network.pos[e[0]][0],network.pos[e[1]][0] ], [ network.pos[e[0]][1],network.pos[e[1]][1] ], color = colors[user_pair], linewidth=1)
                            path_output.append(list(e))
                        else:
                            path = nx.shortest_path(network.G, source=e[0], target=e[1], weight="weight")
                            path_output.append(path)
                            for i in range(len(path)-1):
                                e1 = path[i]
                                e2 = path[i+1]
#                                 plt.plot([ network.pos[e1][0],network.pos[e2][0] ], [ network.pos[e1][1],network.pos[e2][1] ], color = colors[user_pair], linewidth=1)

            print((self.s_list[user_pair],self.t_list[user_pair]),":",path_output)
            
            
            user_pair_src_dst = (self.s_list[user_pair],self.t_list[user_pair])
            
            
            
            path_edges = []
            for edge in path_output:
                path_edges.append((edge[0],edge[len(edge)-1]))
            network.solution_each_user_pair_used_paths[user_pair_src_dst,network.path_counter] = [path_edges]
            
            
            W = this_user_path_W
            E_p =network.compute_e2e_rate(network.path_counter,path_edges,W)
            e2e_F = network.compute_e2e_fidleity(path_edges)
            network.each_path_width[network.path_counter] = W
            network.each_path_e2e_F_value[network.path_counter] = e2e_F
            network.each_path_e_value[network.path_counter] =E_p
            network.path_counter+=1
            network.solution_path_counter+=1
            
            for edge in path_output:
                if edge[0] not in network.solution_used_nodes_by_paths and edge[0] not in user_pair_src_dst:
                    network.solution_used_nodes_by_paths.append(edge[0])
                if edge[1] not in network.solution_used_nodes_by_paths and edge[1] not in user_pair_src_dst:
                    network.solution_used_nodes_by_paths.append(edge[1])
        
        plt.show()
    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




