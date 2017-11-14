from environment_node import *

def beam_search(level,environment):
    if isinstance(level,list):
        b = level                                   # b is this generation, level is the initial set or initial node
    else:
        error('need to pass alist of nodes to initialise beam_search')
    depth = 0

    environment.update_best(level)                  # update best sequence
    
    while len(b)>0:                                 # while there is a next generation
        depth+=1
        l=[]
        if environment.verbosity>0: print('depth : ',depth)
        for n in b:                                 # take te nodes
            l+=branch(n, environment)               # make them bear children
        
        if len(l)>0 : environment.update_best(l)    # update the environment
        environment.update_environment(depth)       

        if environment.verbosity>0: print('size of level before beaming : ', len(l))

        b = environment.beam(l)                     # select some beams, here is where the termination condition dwindles

        if environment.verbosity>0: print('size of level after beaming : ', len(b))
        if environment.verbosity>0: print('best fitness  : ', environment.best_fitness)

    return environment

def branch(parent,environment):
    children=[]   
    child_args_list = parent.give_birth(environment)            # returns a list of arguments to initialise new nodes
    for child_args in child_args_list:
        new_node = node(*child_args,environment =environment)   # initialise the node and do the costly computations
        children.append(new_node)                               # appends new node object to the next level
    parent.filter_children(children,environment)                # here is where the branching factor comes into play
    return children