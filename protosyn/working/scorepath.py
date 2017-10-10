
def DFS(node, scorelist, score=0, n=1, path=[]):
    path = path + [node.index]
    score += (n*11 + 8*node.mass)
    for child in node.connect:
        if child.index not in path:
            DFS(child, scorelist, score, n+1, path)
    if len(node.connect) == 1:
        scorelist.append(score/100.0)
        # print n, path, score


# def DFS2(node, connect, mass, scorelist, score=0, n=1, path=[]):
#     path = path + [node]
#     score += (n*11 + 8*mass[node])
#     for child in connect[node]:
#         if child not in path:
#             DFS2(child, connect, mass, scorelist, score, n+1, path)
#     if len(connect[node]) == 1:
#         scorelist.append(score/100.0)
#         print n, path, score

# mass = [0] + 4*[6] + 10*[1]
# con = [
#     [],             # 0
#     [2, 5,6,7],     # 1
#     [1,8,3,9],      # 2
#     [2,10,4,11],    # 3
#     [3,12,13,14],   # 4
#     [1],            # 5
#     [1],            # 6
#     [1],            # 7
#     [2],            # 8
#     [2],            # 9
#     # [3,12],            # 10
#     [3],            # 10
#     [3],            # 11
#     # [4,10],            # 12
#     [4],            # 12
#     [4],            # 13
#     [4],            # 14
# ]

# # scorepath(1,-1, [], con)
# # DFS(1, con)
# scores = [None]
# atoms = range(1,15)
# for k in atoms:
#     p = []
#     DFS2(k, con, mass, p)
#     scores.append(sorted(p))
#     #print scores[-1]

# # groups = []
# while atoms:
#     pivot = atoms.pop(0)
#     group = [pivot]
#     nongroup = []
#     # groups.append(group)
#     for atom in atoms:
#         if scores[pivot] == scores[atom]:
#             group.append(atom)
#         else:
#             nongroup.append(atom)
#     atoms = nongroup
#     print group

# DFS2(3, con, mass, p, path=[2])
