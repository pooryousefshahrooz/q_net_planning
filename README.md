# Overview
This repository contains the source code and supplementary materials for the paper https://arxiv.org/pdf/2308.16264v2, which explores methods for fault-tolerant distributed quantum computing through advanced network planning techniques.

# Paper Details
Title: Resource Allocation for Rate and Fidelity Maximization in Quantum Networks
Authors: Shahrooz Pouryousef, Hassan Shapourian, Alireza Shabani, Ramana Kompella, Don Towsley
Published: August 2023

Abstract: Existing classical optical network infrastructure cannot be immediately used for quantum network applications due to photon loss. The first step towards enabling quantum networks is the integration of quantum repeaters into optical networks. However, the expenses and intrinsic noise inherent in quantum hardware underscore the need for an efficient deployment strategy that optimizes the allocation of quantum repeaters and memories. In this paper, we present a comprehensive framework for network planning, aiming to efficiently distributing quantum repeaters across existing infrastructure, with the objective of maximizing quantum network utility within an entanglement distribution network. We apply our framework to several cases including a preliminary illustration of a dumbbell network topology and real-world cases of the SURFnet and ESnet. We explore the effect of quantum memory multiplexing within quantum repeaters, as well as the influence of memory coherence time on quantum network utility. We further examine the effects of different fairness assumptions on network planning, uncovering their impacts on real-time network performance.

# Requirements:
To run the code in this repository, you need to have the following libraries installed on your machine:

IBM CPLEX:

IBM CPLEX is an optimization solver for linear programming, mixed integer programming, and quadratic programming. You can download and install CPLEX from the IBM website.
Ensure that you have a valid license for CPLEX. Free academic licenses are available for students and researchers.
NetworkX:

NetworkX is a Python package for the creation, manipulation, and study of complex networks.
You can install NetworkX using pip:
pip install networkx

Ensure that you have IBM CPLEX installed as per the instructions above.

# Usage
Prepare your environment:

Set up your IBM CPLEX environment variables as required.
Ensure NetworkX is installed and can be imported in your Python environment.
Run the scripts:

You can find various scripts in the src directory that correspond to different experiments and analyses described in the paper.
For example, to run the main quantum network planning script, use:

# Configuration:

config.py contains the parameters such as the topology, assumptions. You can modify these files to run different scenarios and network configurations.

For any questions or issues, please open an issue on the GitHub repository or contact the authors via email.

