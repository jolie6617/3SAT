from random import sample, choice, shuffle
from itertools import product
from numpy import arange
from time import time


class SAT():
    def __init__(self):
        # To initialize an empty 3SAT instance: phi
        self.variables = []  # Names of the variables
        self.clauses = []  # Clauses
        self.v = 0  # Number of variables = len(self.variables)
        self.c = 0  # Number of clauses  = len(self.clauses)
        self.arrangement = []  # Current value assignment of the variables

    def __repr__(self):
        # To print out the SAT instance in readable form
        r = "SAT instance with " + str(self.v) + " variables and " + str(self.c) + " clauses: \n"
        for c in self.clauses:
            r += str(c) + "\n"
        return r

    def set_variable(self, variable, value):
        # Set the given 'variable' to the given  'value'
        self.arrangement[abs(variable) - 1] = value

    def flip_variable(self, variable):
        # Flip the value of 'variable'
        self.arrangement[abs(variable) - 1] = not self.arrangement[abs(variable) - 1]

    def get_variable(self, variable):
        # Return the value assigned to 'variable'
        return self.arrangement[abs(variable) - 1]

    #    Evaluation methods

    def clause_value(self, clause):
        # Return the value of 'clause'
        x1 = self.get_variable(clause[0])
        x2 = self.get_variable(clause[1])
        x3 = self.get_variable(clause[2])
        # if the value is negative, negate boolean value of the variable
        if clause[0] < 0: x1 = not x1
        if clause[1] < 0: x2 = not x2
        if clause[2] < 0: x3 = not x3
        return (x1 or x2 or x3)

    def value(self):
        # To evaluate phi. It will stop as soon as a clause is False
        for c in self.clauses:
            # Any false clause will make phi unsatisfied
            # Else, it is satisfied
            if self.clause_value(c) == False:
                return False
        return True

    def value_full(self):
        # To evaluate phi. It does not stop until all clauses are evaluated
        val = True
        for c in self.clauses:
            val &= self.clause_value(c)
        return val

    def unsatisfied_ratio(self):
        # To calculate the ratio of satisfied clauses
        count = 0
        for c in self.clauses:
            if self.clause_value(c) == False:
                count += 1
        return count / self.c

    #    Random construction methods

    def randomized_clause(self):
        # To generate a clause randomly. Acts as helper function
        c = sample(self.variables, 3)
        c[0] *= choice([-1, 1])
        c[1] *= choice([-1, 1])
        c[2] *= choice([-1, 1])
        return c

    def randomized_instance(self, v, c):
        # To build a random 3SAT instance with v variables and c clauses
        self.v = v
        self.variables = list(range(1, v + 1))
        self.arrangement = [False] * v
        self.c = c
        # To generate c clauses: (... or ... or ...)
        self.clauses = []
        for c in range(c):
            self.clauses.append(self.randomized_clause())

    def randomized_arrangement(self):
        # To build a random arrangement (= assignment of the variables)
        # Return value of objective function.
        self.arrangement = [choice([True, False]) for v in range(self.v)]

    def randomized_yes_instance(self, v, c):
        # To build a random 'yes' 3SAT instance with v variables and c clauses
        # choose a random arrangement to satisfy
        self.v = v
        self.variables = list(range(1, v + 1))
        self.randomized_arrangement()

        # choose satisfying clauses
        self.c = c
        self.clauses = []
        for nc in range(self.c):
            # sample(self.variables,3)
            c = self.randomized_clause()
            # To alter c until it becomes satisfying
            while self.clause_value(c) == False:
                c[choice([0, 1, 2])] *= -1  # by flipping terminals (negating)
            self.clauses.append(c)
        if self.value() != True:
            print("something has gone wrong!")

    #    Solution methods

    def exhaustive_search(self):
        # Solve in the 3SAT problem instance using exhaustive search
        # iterate over all the possible Boolean variable assignments
        for self.arrangement in product([True, False], repeat=self.v):
            if self.value() == True:
                return True
        return False

    def full_exhaustive_search(self):
        # Solve in the 3SAT problem instance using exhaustive search
        # iterate over all the possible Boolean variable assignments
        self.decision = False
        for self.arrangement in product([True, False], repeat=self.v):
            self.decision |= self.value_full()
        return self.decision

    def randomized_greedy(self, block_size=2):

        # Solve in the 3SAT problem instance using a randomized greedy search.
        # Find the variable that appears most often and assign it accordingly to maximize
        # Randomize the order of literals_occurance

        literals_occurance = [[i, 0] for i in range(-self.v, self.v + 1)]
        for c in self.clauses:
            for v in c:
                literals_occurance[v + self.v][1] += 1
        literals_occurance.sort(key=lambda a: a[1])

        # Randomize list block by block
        for i in range(0, len(literals_occurance), block_size):
            tmp = literals_occurance[i:i + block_size]
            shuffle(tmp)
            literals_occurance[i:i + block_size] = tmp
        for v in literals_occurance:
            self.set_variable(v[0], v[0] > 0)
        return self.unsatisfied_ratio()

    def GRASP(self, repetition_max=100):
        # To solve the 3SAT problem instance using the GRASP
        # start with the worst solution that satisfies no clauses
        best = 1
        for i in range(repetition_max):
            candidate = self.randomized_greedy(choice([2, 3, 4]))

            # Flip variable assignment and the improvement by local search
            # Select the best candidate from the flipped attempts
            best_v = 0
            for v in self.variables:
                self.flip_variable(v)  # flipping variable
                a = self.unsatisfied_ratio()
                if a < candidate:
                    candidate = a
                    best_v = v
                self.flip_variable(v)  # undo, to try flipping next variable
            if candidate < best:
                best = candidate
        return best

# Testing


instance = SAT()  # global variable... [TODO?]


def output_print(file, line):
    print(line)
    file.write(line + '\n')


def exhaustive_test():
    with open("data_exhaustive.csv", "w") as f:
        # exhaustive search
        output_print(f, "n\tExhaustive")
        repeats_max = 50
        n = 10
        t0 = t1 = 0
        while t1 - t0 < 3600:  # in seconds; if it takes too long then stop testing
            t0 = time()
            for repeats in range(repeats_max):  # e.g. average over 100 instances
                instance.randomized_yes_instance(n, 3 * n)  # rho=3
                decision = instance.exhaustive_search()
            t1 = time()
            # record average time
            output_print(f, str(n) + "\t" + str((t1 - t0) / repeats_max))
            n += 1


def grasp_test():
    with open("data_approx.csv","w") as f:
        # test greedy search
        output_print(f, "rho\tGRASP")
        repeats_max = 100
        v = 10
        # rho = c/v
        for rho in arange(1,7,0.2):
            c = int(rho*v)
            # quality measurement
            q = [0]*3
            for repeat in range(repeats_max):
                instance.randomized_yes_instance(v, c)
                q[0] += instance.GRASP()
            for i in range(3):
                q[i] = '{:1.2f}'.format(1-q[i]/repeats_max)
            output_print(f, str(rho)+'\t'+"\t".join(q))


# exhaustive_test()
grasp_test()