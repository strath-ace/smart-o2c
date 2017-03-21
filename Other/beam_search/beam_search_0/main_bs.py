from __future__ import print_function
from environment_node import *
from beam_search import beam_search
#!/usr/env python


# set the parameters of the problem and initialise environment
budget = 16
p = 8
branching_factor = 0        # the max. number of children that a node can bear, filtered after computing attributes, 0 for no filter
beaming_factor = 10000      # the max. number of members in a level that will bear children in next generation, 0 for no filter
verbosity = 1               # now is just 0 for nothing or 1 for some basic stuff
initial_env = environment(p,budget,branching_factor,beaming_factor,verbosity)

# set the parameters of the node and initialise node
sequence = [(0,0)]
fitness = 0
initial_node = node(sequence,fitness,initial_env)

# do the beam search
env = beam_search([initial_node],initial_env)

# print the best
print('=======================')
print('fitness : ', env.best_fitness)
print('n. sol. : ', len(env.best))
# for sol in env.best:
#     print(sol.visited)