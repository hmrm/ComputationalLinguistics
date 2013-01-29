from __future__ import division
import sys
import argparse as ap
parser = ap.ArgumentParser(description='HMM, input expected from STDIN')
parser.add_argument('-n', action='store', dest='n_states', type=int, default=2, help='number of states')
n_states = parser.parse_args().__dict__['n_states']

def memoize(f):
    cache = {}
    def mf(*args):
        if args not in cache:
            cache[args] = f(*args)
        return cache[args]
    return mf

def main():
    data = [x.strip() for x in sys.stdin.read().split()]
    a = distribute_transitions(data) #transition matrix
    b = distribute_emissions(data)
    pi = distribute_initial_states(data)
    hmm = HMM(data, a, b, pi)
    for l, cs in hmm.soft_counts().iteritems():
        print l + ',' + str(cs[0]) + ',' + str(cs[1])

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

    def soft_counts(self):
        countdict = {x : [0 for i in xrange(n_states)] for x in set('#'.join(self.data))}
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
        return (self.alpha(i,t, string) * self.a[i][j] * self.b[i][j][string[t]] * self.beta(j,t+1, string)) \
            / (sum( \
                sum(self.alpha(m, t, string) * self.a[m][n] * self.b[m][n][string[t]] * self.beta(n, t+1, string) for n in xrange(n_states)) \
                    for m in xrange(n_states)))

    def gamma(self, i, t):
        pass

    def pi_hat(self, i):
        pass

    def a_hat(self, i,j):
        pass

    def b_hat(self, i,j,k):
        pass


'''
Generates transition matrix from data, these are currently extremely basic, and will be improved in Part 2
'''
def distribute_transitions(data):
    return [[1 / n_states for i in xrange(n_states)] for j in xrange(n_states)] #assigns equal transition chance to all states

def distribute_initial_states(data):
    return [1 / n_states for i in xrange(n_states)] #equal chance to each state

def distribute_emissions(data):
    alphabet = set('#'.join(data))
    return [[{x: 1 / len(alphabet) for x in alphabet} for i in xrange(n_states)] for j in xrange(n_states)]

if __name__ == "__main__":
    main()
