# Overview
This repository contains the source code and supplementary materials for the paper https://arxiv.org/pdf/2308.16264v2, which explores resource allocation for quantum network planning.

# Paper Details
Title: Resource Allocation for Rate and Fidelity Maximization in Quantum Networks
Authors: Shahrooz Pouryousef, Hassan Shapourian, Alireza Shabani, Ramana Kompella, Don Towsley
Published: August 2023

Abstract: Existing classical optical network infrastructure cannot be immediately used for quantum network applications due to photon loss. The first step towards enabling quantum networks is the integration of quantum repeaters into optical networks. However, the expenses and intrinsic noise inherent in quantum hardware underscore the need for an efficient deployment strategy that optimizes the allocation of quantum repeaters and memories. In this paper, we present a comprehensive framework for network planning, aiming to efficiently distributing quantum repeaters across existing infrastructure, with the objective of maximizing quantum network utility within an entanglement distribution network. We apply our framework to several cases including a preliminary illustration of a dumbbell network topology and real-world cases of the SURFnet and ESnet. We explore the effect of quantum memory multiplexing within quantum repeaters, as well as the influence of memory coherence time on quantum network utility. We further examine the effects of different fairness assumptions on network planning, uncovering their impacts on real-time network performance.



<img src="https://github.com/pooryousefshahrooz/q_net_planning/blob/main/data/esnet.png" align="center" height="400" width="700"/>



# Requirements:
To run the code in this repository, you need to have the following libraries installed on your machine:

* IBM CPLEX

* NetworkX:
  
You can install NetworkX using pip:

```pip install networkx```

Make sure you have the IBM CPLEX academic edition installed, as the community edition may not support the problem size.

# Usage

* Various scripts corresponding to different experiments and analyses described in the paper are located in the root directory. The primary script is ```main.py```.

* The ```config.py``` file includes the configuration and hyperparameters for the experiments. You can modify the hyperparamter values in this file to conduct various experiments with different network topologies, such as a repeater chain or SURFnet and ESnet.

For any questions or issues, please open an issue on the GitHub repository or contact the authors via email.

  
# Running an Experiment (figure 2)
Here we show how to run the repeater placement on a link for a dumbbell topology. We have shown the most important hyperparameters in the config.py file here to run this experiment:

* general hyperparameters
  ```repeating_times = 1 # repeatting the experiment 
  include_fidelity_in_utility = True 
  q_values = [0.5]# Different values of q (swap success probability) that we want to evaluate in the experiment
  utility_type = "NGTV" # Nagativity utility function```

  
* solver setup
  ```schemes = ["exhaustive"]#["link_based","exhaustive","linear_link_based"] # different ways of solving the optimization problem
  link_based_solver = "CPLEX"#"exhaustive"#"Bonmin","Baron
  relaxing_QCAST_formulation = True```
    
    
* network topology
  ```network_topology = "Dumbbell"# Dumbbell,Random,"ESnet2.gml",Repeater_chain,SurfnetCore.gml,"Random
  edge_F_values = [1.0]# set of values to experiment for link level fidelity. we assume all links have this fidelity
  set_of_number_of_user_pairs = [3] # number of user pairs in the network
  number_of_paths_values = [6000] # number of paths between each pair of user pairs
  K_values = [1] # set of number of paths allowed to be used for each user pair. One means each user pair can use at most one path.
  lengths_of_middle_link = [100] # in km for dumbell shape topology```
    
    
* network planning assumptions experiment hyperparameters
  
  ```dynamic_system_flag = False```
  


    
* repeaters and end nodes hyperparameters
  ```D_values= [10] # set of values for repeaters memory budget
  checking_repeater_memory_life_time_flag = False # set to True of we want to restrict paths 
  checking_end_node_memory_life_time_flag = False 
  R_values = [10,4,6,8] # set of values for number of repetaers budget
  end_user_memory_set = [10] # the memory budget of end nodes```
  
# Interpreting the Results
We report some or all of the following metrics in our experiments. Check function ```save()``` in ```network.py``` script to see the key for each column of the results ```csv``` file.

* maximum utility of the network
* processing time in seconds
* edges of each path used by each user pair in the optimal solution
* end-to-end fidelity of each path
* end-to-end rate on each path used by each user pair in the optimal solution

  
# Citing

Please cite our paper if you're using any part of this code for your project.

```@inproceedings{pouryousef2023quantum,
title={Quantum Network Planning for Utility Maximization},
  author={Pouryousef, Shahrooz and Shapourian, Hassan and Shabani, Alireza and Towsley, Don},
  booktitle={Proceedings of the 1st Workshop on Quantum Networks and Distributed Quantum Computing},
  pages={13--18},
  year={2023}
}```

```@inproceedings{pouryousef2023quantum,
  title={Resource Allocation for Rate and Fidelity Maximization in Quantum Networks},
  author={Pouryousef, Shahrooz and Shapourian, Hassan and Shabani, Alireza, Kompella, Ramana and Towsley, Don},
  journal={arXiv preprint arXiv:2308.16264v2},
  year={2024}
}```



