from __future__ import division
import sys
from collections import deque

debug = True
# expects the two strings to be the first and second arguments to the script
str1 = sys.argv[1]
str2 = sys.argv[2]

#this reduces to a pathing problem on a directed graph, so I'm implementing A*
def getcost(a, b):
    vowels = ['a', 'e', 'i', 'o', 'u']
    if a == b:
        return 0
    if a in vowels and b in vowels:
        return .5
    if a not in vowels and b not in vowels:
        return .6
    else:
        return 1.2

def neighbors(a, b):
    if len(a) > 0 and len(b) > 0:
        return [((a[1:], b), 1, (a[0], ' ')), 
                ((a, b[1:]), 1, (' ', b[0])), 
                ((a[1:], b[1:]), getcost(a[0], b[0]), (a[0], b[0]))]
    elif len(a) > 0:
        return [((a[1:], b), 1, (a[0], ' '))]
    elif len(b) > 0:
        return [((a, b[1:]), 1, (' ', b[0]))]
    else:
        return []

consideration = [(str1, str2)]
costs = {(str1, str2) : 0}
parents = {(str1, str2) : None}
stringchange = {(str1, str2) : ('', '')}
while 1:
    consideration = sorted(consideration, key=lambda (a, b): len(set(a).symmetric_difference(set(b)))*.25 + costs[(a,b)])
    current = consideration[0]
    if current == ('',''):
        break
    consideration = consideration[1:]
    for neighbor, incrementalcost, strings in neighbors(*current):
        cost = costs[current] + incrementalcost
        if neighbor not in costs or cost < costs[neighbor]:
            costs[neighbor] = cost
            parents[neighbor] = current
            stringchange[neighbor] = strings    
            if neighbor not in consideration:
                consideration.append(neighbor)

path = [current]
while parents[current] != None:
    path.append(stringchange[current])
    current = parents[current]
path.reverse()

print "String Edit Distance: " + str(costs[('','')])
print "Matching:"
print "".join([a[0] for a in path])
print "".join([a[1] for a in path])
