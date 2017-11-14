
class environment:                          # represents the problem
    def __init__(self,p,initial_budget,branching_factor=0,beaming_factor=0,verbosity=1):
        self.best = []                      # it stores the best
        self.best_fitness = []
        self.p = p                          # and attributes that never change
        self.initial_budget = initial_budget
        self.branching_factor = branching_factor
        self.beaming_factor = beaming_factor
        self.verbosity = verbosity
        self.depth = 0                      # and attributes that change according to depth
        self.budget = initial_budget

    def update_best(self,level):            # this method updates the best
        # in this case the complete list of sequences with maximum fitness
        maximum = max( level , key=lambda x : x.fitness)
        best_fitness = maximum.fitness
        best = [node for node in level if node.fitness == best_fitness]
        self.best = best                    # stores best
        self.best_fitness = best_fitness    # and best fitness

    def update_environment(self,depth):     # this method updates the environment according to depth
        self.depth = depth
        self.budget = self.initial_budget - depth

    def beam(self,level):
        if self.beaming_factor>0:
            # print min(self.beaming_factor,len(level))
            return sorted(level,key= lambda x:-x.fitness)[:min(self.beaming_factor,len(level))]
        else:
            return level


class node:
    def __init__(self,visited,fitness,environment):
        self.visited = visited               # store the basic attributes passed as arguments
        self.coordinates = visited[-1]
        self.old_fitness = fitness
        self.compute_attributes(environment)  # do the costly computations

    def compute_attributes(self,environment): # do the costly computations
        self.fitness = self.old_fitness + abs(self.coordinates[0])**2 + abs(self.coordinates[1])**2

    def give_birth(self,environment):         # returns a list of next level node argument lists necessary for initialisation
        children = []
        # in this case we return coordinates possible by moving in one of the directions of the keyboard
        visited = self.visited
        fitness = self.fitness
        coordinates = self.visited[-1]
        # for delta in [[-1,0],[1,0],[0,-1],[0,1]]:
        for delta in [[-1,0],[1,0],[0,-1],[0,1],[-1,-1],[1,-1],[-1,1],[1,1]]:
            candidate_coordinates = (coordinates[0]+delta[0],coordinates[1]+delta[1]);
            child_args = [visited+[candidate_coordinates],fitness]; # construct the child arguments will be unpacked as *child_args
            if self.is_child_legal(child_args,environment):         # check if the child is legal before costly computations
                children.append(child_args)                         # apend the arguments to the list and a new node will be created 
        
        # print('yo', children) 
        return children

    def is_child_legal(self,child_args,environment):                # check if the child is legal before costly computations
        child_coordinates = child_args[0][-1]
        # print('cc', child_coordinates )
        return (( 0 <= child_coordinates[0] < environment.p ) and ( 0<= child_coordinates[1] < environment.p) and (environment.budget >= 1) and (child_coordinates not in self.visited ))

    def filter_children(self,children,environment):                 # here is where the branching factor comes into play
        if 0 < environment.branching_factor < len(children):
            return sorted(children,key= lambda x:-x.fitness)[environment.branching_factor]
        else:
            return children