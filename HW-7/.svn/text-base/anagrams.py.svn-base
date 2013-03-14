import sys
from collections import Counter

def overlap(a, b):
    def getadjset(x):
        return set(x[i-1:i+1] for i in xrange(1,len(x)))
    return len(getadjset(a).intersection(getadjset(b)))

def multioverlap(xs):
    return sum(overlap(a,b) for a in xs for b in xs if b != a)

lines = [line.strip() for line in sys.stdin.readlines()]

anagrams = {}
for line in lines:
    c = Counter(line)
    c = tuple([(a, c[a]) for a in c])
    if c in anagrams:
        anagrams[c].append(line)
    else:
        anagrams[c] = [line]

relevantKeys = filter(lambda x: len(anagrams[x][0]) >= 6 and len(anagrams[x]) > 1, anagrams.keys())
relevantKeys = sorted(relevantKeys, key=lambda x: multioverlap(anagrams[x]))
relevantKeys = sorted(relevantKeys, key=lambda x: len(anagrams[x][0]))
relevantKeys = sorted(relevantKeys, key=lambda x: len(anagrams[x]))
for key in relevantKeys[::-1]:
    print anagrams[key]
