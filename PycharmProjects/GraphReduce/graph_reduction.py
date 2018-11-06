# This program takes a text file in the form of an adjacency matrix as input
# Cycles are found and reduced in the algorithm

import re


global edges, cyclenumber, graph

n = 100000
graph = [[] for i in range(n)]
cycles = [[] for i in range(n)]

edges = 0
cyclenumber = 0


def main():
    global edges, cyclenumber, graph

    #f = "testMatrix.txt"
    #f = "test2.txt"
    #f = "doubleBondTest.txt"
    f = "hTest.txt"

    matrix = create_matrix(f)

    color = [0] * n
    par = [0] * n

    mark = [0] * n

    cyclenumber = 0

    dfs_cycle(1, 0, color, mark, par)

    graph = matrix

    cycles = print_cycles(edges, mark)

    #f = open("cycles.txt", 'w')


    reduce_graph(matrix, cycles)



# cycle_exists(matrix)


# open file and create matrix based on textfile
def create_matrix(file_name):
    global edges, graph

    f = open(file_name, 'r')

    # number of rows and columns not including attachments to the protein
    firstLine = f.readline()
    count = 0
    for element in firstLine:
        if element == '0' or element == '1':
            count = count + 1

    # subtract 1 to account for bond to protein
    count = count - 1

    m = [[0 for x in range(count)] for y in range(count)]

    # reset file pointer
    f.seek(0)

    # Fill matrix with data from the input file
    i = 0

    for line in f:
        if i == count:
            break
        j = 0
        # print(line)
        while j < count:
            m[i][j] = line[j * 2]
            if m[i][j] == "1":
                add_edge(graph, i, j)
            j += 1
        i += 1

    f.close()
    edges = count

    return m


# --------

def dfs_cycle(u, p, color, mark, par):
    global cyclenumber
    # If you already completely visited a vertex

    if color[u] == 2:
        return

    # Seen vertex, but was not completely visited, cycle detected
    # backtrack based on parents to find the complete cycle

    if color[u] == 1:

        cyclenumber += 1
        cur = p
        mark[cur] = cyclenumber

        # backtrack the vertex which are in the current cycle thats found
        while cur != u:
            cur = par[cur]
            mark[cur] = cyclenumber
        return
    par[u] = p

    # partially visited
    color[u] = 1

    # simple dfs on graph
    for v in graph[u]:

        # if it has not been visited previously
        if v == par[u]:
            continue

        dfs_cycle(v, u, color, mark, par)

    # completely visited
    color[u] = 2



# Add edge to graph
def add_edge(g, u, v):
    g[u].append(v)

# Change edge
def replace_edge(g, u, v):
    g[u][v] = '1'
    g[v][v] = '1'

def fix_edge(g, u):
    g[u][u] = '0'


def print_cycles(edges, mark):
    global cyclenumber, graph

    cycle_file = open("cycle_file.txt", "w")


    # Push the edges into adjacency list

    i = 0
    for i in range(edges):
        if mark[i] != 0:
            cycles[mark[i]].append(i)
            #print 'Cycles:    ' + str(i)
            #print cycles[mark[i]].append(i)\

    # Combines cycles of 3 and 5 if they share 2 edges
    # These should be combined because it means a double peptide bond was found
    count = 0
    for j in range(len(cycles)):
        if cycles[j] != []:
            count += 1
            if len(cycles[j]) == 3:
                for k in range(len(cycles)):
                    if len(cycles[k]) == 5:
                        connections = 0
                        # If a size 3 and 5 cycle share two connections then combine them
                        for nC in cycles[j]:
                            for nC2 in cycles[k]:
                                if graph[nC][nC2] == '1':
                                    connections = connections + 1
                        if connections == 2:
                            cycles[k].append(cycles[j][0])
                            cycles[k].append(cycles[j][1])
                            cycles[k].append(cycles[j][2])

                            # Remove the cycle of 3 because it was appended to the 5
                            cycles[j] = []


    for m in range(count+1):
        if cycles[m] != []:
            cycle_file.write(str(cycles[m]) + "\n")

    return cycles
"""
    j = 0
    for c in cycles:
        if c != []:
            if len(c) == 3:
                k = 0
                for c2 in cycles:
                    if len(c2) == 5:
                        for e in c2:
                            cycles[]
            print "CYCLES"
            print(c)
        j += 1
"""


# Creates a new matrix by removing all but one node from each cycle
def reduce_graph(matrix, cycles):

    new_matrix = matrix
    first = 0
    list_to_remove = set()

    for c in cycles:
        if c != []:
            i = 0
            for node in c:
                if i != 0:
                    list_to_remove.add(node)
                    replace_edge(new_matrix, first, node)

                    j = 0
                    for j in range(len(matrix)):
                        if (matrix[node][j]) == '1':
                            replace_edge(new_matrix, first, j)
                else:
                    first = node
                    i = 1

    remove_columns(new_matrix, list_to_remove)



# removes the columns from the cycles
def remove_columns(matrix, list_to_remove):

    num_rows = len(matrix)
    num_columns = len(matrix)

    new_file = open("reduced_matrix.txt", "w")


    for i in range(num_rows):
        line = ''
        if not (i in list_to_remove):
            for j in range(num_columns):
                if (i == j):
                    #fix_edge(matrix)
                    if j == num_columns-1:
                        line += '0'
                    else:
                        line += '0,'
                else:
                    curr_element = matrix[i][j]
                    if not (j in list_to_remove):
                        if j == num_columns-1:
                            line += (curr_element)
                        else:
                            line += (curr_element + ',')
            #line = re.sub('[[]]')
            new_file.write(line + "\n")


if __name__ == "__main__":
    main()
