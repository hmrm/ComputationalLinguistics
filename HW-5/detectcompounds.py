#This is an attempt at using a very direct probabilistic approach: counting up the occurences of a substring as a word versus not
from __future__ import division
import sys

def substrings(string):
    ret = []
    for i in xrange(len(string)):
        for j in xrange(1, len(string) - i + 1):
            ret.append(string[i:i + j])
    return ret

with open("final_gold_standard.txt") as f:
    goldlines = f.readlines()

goldlines = [item.strip() for sublist in (x.replace(' ', '-').strip().split() for x in goldlines) for item in sublist]

golddict = {}

for line in goldlines:
    justletters = filter(str.isalnum, line)
    golddict[justletters] = ('-' in line) #if it's a compound, true, else, false

substrdict = {}
print("Here\n")
corpus = sys.stdin.readlines()
#print("Here\n")
lines = [item for sublist in (x.strip().split('-') for x in corpus) for item in sublist]
for line in lines:
    
    data = substrings(line)
    for substring in data:
        if not substring in substrdict:
            substrdict[substring] = {True : 0, False : 0}
        if substring != line:
            substrdict[substring][False] += 1
        else:
            substrdict[substring][True] += 1

resultdict = {}


for line in (x.strip() for x in corpus):
    if '-' in line:
        resultdict[filter(str.isalnum, line)] = True
    else:
        notcompound = 1.0
        for i in xrange(1, len(line)):
            a = line[:i]
            b = line[i:]

            splithere = (substrdict[a][True]/(substrdict[a][True] + substrdict[a][False])) * (substrdict[b][True]/(substrdict[b][True] + substrdict[b][False]))
            notsplithere = (substrdict[a][False]/(substrdict[a][True] + substrdict[a][False])) * (substrdict[b][False]/(substrdict[b][True] + substrdict[b][False]))

            #normalization, notsplithere is now a probability
            notsplithere /= splithere + notsplithere
            
            notcompound *= notsplithere
        if notcompound > .5:
            resultdict[filter(str.isalnum, line)] = False
        else:
            resultdict[filter(str.isalnum, line)] = True

print "Results:"

for x in resultdict:
    print x,
    if resultdict[x]:
        print " : Compound"
    else:
        print " : Atomic"

true_positives = 0
total_positives = 0
true_data = 0
for x in resultdict:
    if x in golddict:
        if resultdict[x]:
            total_positives += 1
        if golddict[x]:
            true_data += 1
        if resultdict[x] and golddict[x]:
            true_positives += 1

print "precision: " + str(true_positives / total_positives)
print "recall: " + str(true_positives / true_data)
