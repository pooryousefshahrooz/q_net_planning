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

