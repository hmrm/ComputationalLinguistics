from __future__ import division
import sys
import argparse
import random
parser = argparse.ArgumentParser(description='HMM, input expected from STDIN')
parser.add_argument('-n', action='store', dest='n_states', type=int, default=2, help='number of states')
parser.add_argument('-d', action='store', dest='max_delta', type=float, default=.01, help='maximum distance between values for stopping iteration')
parsed = parser.parse_args()
n_states = parsed.n_states
max_delta = parsed.max_delta

def memoize(f):
    cache = {}
    def mf(*args):
        if args not in cache:
            cache[args] = f(*args)
        return cache[args]
    return mf

def main():
    data = ['#' + x.strip() + '#' for x in sys.stdin.read().split()]
    a = distribute_transitions(data) #transition matrix
    b = distribute_emissions(data)
    pi = distribute_initial_states(data)
    sys.stderr.write("Initial distributions complete\n")
    hmm = HMM(data, a, b, pi)
    next_hmm = hmm.next_model()
    while not hmm.close_enough(next_hmm, max_delta):
        hmm = next_hmm
        next_hmm = next_hmm.next_model()
    next_hmm.print_groups()
    sys.stderr.write('DONE\n')

def safe_divide(a,b):
    try:
        return a / b
    except ZeroDivisionError:
        return 0

class HMM():
    def __init__(self, data, a, b, pi):
        self.data = data #input string
        self.a = a #state transition matrix
        self.b = b #emission matrix
        self.pi = pi #initial state matrix

    def close_enough(self, model, delta):
        max_dist = 0
        for i in xrange(n_states):
            max_dist = max(abs(self.pi[i] - model.pi[i]), max_dist)
            for j in xrange(n_states):
                max_dist = max(abs(self.a[i][j] - model.a[i][j]), max_dist)
                for k in self.b[i][j]:
                    max_dist = max(abs(self.b[i][j][k] - model.b[i][j][k]), max_dist)
        sys.stderr.write("distance between models: " + str(max_dist) + "\n")
        return max_dist < delta

    def group_letters(self):
        groups = [[] for x in xrange(n_states)]
        counts = self.letterwise_soft_counts()
        for x in counts:
            groups[counts[x].index(max(counts[x]))].append(x)
        return groups
            
    def print_groups(self):
        sys.stderr.write('\n')
        sys.stderr.write('Groups for this model:\n')
        groups = self.group_letters()
        for i in xrange(len(groups)):
            sys.stderr.write(str(i) + ': ' + ''.join(groups[i]) + '\n')
        sys.stderr.write('\n')

    def letterwise_soft_counts(self):
        sys.stderr.write("Initializing soft count\n")
        countdict = {x : [0 for i in xrange(n_states)] for x in set(''.join(self.data))}
        for string in self.data:
            for t in xrange(len(string)):
                for i in xrange(n_states):
                    countdict[string[t]][i] += sum(self.p(t, i, j, string) for j in xrange(n_states))
        return {x : [safe_divide(countdict[x][i] , sum(countdict[x][j] for j in xrange(n_states))) for i in xrange(n_states)] for x in countdict}

    @memoize
    def alpha(self, j, t, string):
        if t == 0:
            return self.pi[j]
        else:
            return sum(self.alpha(i,t-1, string) * self.a[i][j] * self.b[i][j][string[t-1]] for i in xrange(n_states))

    @memoize
    def beta(self, i, t, string):
        if t == len(string):
            return 1
        else:
            return sum(self.a[i][j] * self.b[i][j][string[t]] * self.beta(j, t+1, string) for j in xrange(n_states))

    @memoize
    def p(self, t, i, j, string):
        return (self.alpha(i,t, string) * self.a[i][j] * self.b[i][j][string[t]] * self.beta(j,t+1, string)) / (sum(sum(self.alpha(m, t, string) * self.a[m][n] * self.b[m][n][string[t]] * self.beta(n, t+1, string) for n in xrange(n_states)) for m in xrange(n_states)))

    @memoize
    def gamma(self, i, t, string):
        return (self.alpha(i,t,string) * self.beta(i,t,string)) / sum(self.alpha(j,t,string) * self.beta(j,t,string) for j in xrange(n_states))

    def a_hat(self, i, j):
        e_transitions_ij = 0
        e_transitions_i = 0
        for string in self.data:
            e_transitions_ij += sum(self.p(t, i, j, string) for t in xrange(len(string)))
            e_transitions_i += sum(self.gamma(i, t, string) for t in xrange(len(string)))
        return e_transitions_ij / e_transitions_i

    def b_hat(self, i, j, k):
        e_transitions_ijk = 0
        e_transitions_ij = 0
        for string in self.data:
            e_transitions_ijk += sum(self.p(t, i, j, string) for t in xrange(len(string)) if string[t] == k)
            e_transitions_ij += sum(self.p(t, i, j, string) for t in xrange(len(string)))
        return e_transitions_ijk / e_transitions_ij

    def pi_hat(self, i):
        s_sum = 0
        count = 0
        for string in self.data:
            s_sum += self.gamma(i, 0, string)
            count  += 1
        return s_sum / count

    def stringwise_pi_hat(self, i, string):
        return self.gamma(i, 1, string)

    def stringwise_a_hat(self, i,j, string):
        return sum(self.p(t,i,j, string) for t in xrange(len(string))) / sum(self.gamma(i, t, string) for t in xrange(len(string)))

    def stringwise_b_hat(self, i,j,k, string):
        return sum(self.p(t,i,j,string) for t in xrange(len(string)) if string[t] == k) / sum(self.p(t,i,j, string) for t in xrange(len(string)))

    def next_model(self):
        print "Determining next model"
        data = self.data
        a = [[self.a_hat(i, j) for j in xrange(n_states)] for i in xrange(n_states)]
        b = [[{k : self.b_hat(i, j, k) for k in self.b[i][j]} for j in xrange(n_states)] for i in xrange(n_states)]
        pi = [self.pi_hat(i) for i in xrange(n_states)]
        return HMM(data, a, b, pi)
'''
Generates transition matrix from data, these are currently extremely basic, and will be improved in Part 2
'''
def distribute_transitions(data):
    ret = [[random.random() for i in xrange(n_states)] for j in xrange(n_states)]
    for i in xrange(n_states):
        psum = sum(ret[i][j] for j in xrange(n_states))
        for j in xrange(n_states):
            ret[i][j] /= psum
    return ret

def distribute_initial_states(data):
    ret = [random.random() for i in xrange(n_states)] #equal chance to each state
    psum = sum(ret)
    return [x / psum for x in ret]

def distribute_emissions(data):
    alphabet = list(set('#'.join(data)))
    ret = [[{x : random.random() for x in alphabet} for i in xrange(n_states)] for j in xrange(n_states)]
    for lst in ret:
        for d in lst:
            psum = sum(d.values())
            for x in d:
                d[x] /= psum
    return ret

if __name__ == "__main__":
    main()
