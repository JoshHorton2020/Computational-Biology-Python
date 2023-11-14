from graphviz import Digraph

internal_count = 0
output_list = []
def makeDPrime(distance_matrix): 
    output = []
    nDeDuo = len(distance_matrix)-2
    for outerdex, row in enumerate(distance_matrix): 
        temp = []
        for index, col in enumerate(row): 
            if outerdex == index: 
                num = 0
            else: 
                num = nDeDuo * col - sum(row) - sum(x[index] for x in distance_matrix)
            temp.append(num)
        output.append(temp)
    return output

def makeD(dPrime, distance_matrix): 
    #find min 
    h1 = [x for x in range(len(dPrime))]
    minVal = min([min(x) for x in dPrime])
    for index, l in enumerate(dPrime): 
        if minVal in l: 
            mindex = [index, l.index(minVal)]
            break
    header = [x for x in h1 if x not in mindex]
    nh = [x for x in h1 if x in mindex]
    #print(header)
    newMatrix = []
    for x in range(len(header)+1):
        newMatrix.append([0 for x in range(len(header)+1)])
    #normal distances
    for index, x in enumerate(header): 
        for innerdex, y in enumerate(header):
            newMatrix[index][innerdex] = distance_matrix[x][y]
    #distances to m 
    for index, x in enumerate(newMatrix):
        if index == len(newMatrix)-1:
            for innerdex, y in enumerate(x): 
                newMatrix[index][innerdex] = newMatrix[innerdex][index]
        else: 
             x[-1] = (int) ((distance_matrix[header[index]][mindex[0]] + distance_matrix[header[index]][mindex[1]] - distance_matrix[mindex[0]][mindex[1]])/2)
    dist = storeNodes(distance_matrix, mindex, header, nh)
    #dist = dist + newMatrix[-1]
    #for x in range(len(dist)-1):
    #   og[x].append(dist[x])
    #og.append(dist)
    #print(nh)
    return [newMatrix, nh]

def storeNodes(distance_matrix, mindex, header, nh):
    #calc distances from a and b to internal node 
    global internal_count
    internal_count += 1
    d1 = (int) ((distance_matrix[mindex[0]][mindex[1]] + distance_matrix[mindex[0]][header[0]] - distance_matrix[mindex[1]][header[0]])/2)
    d2 = (int) ((distance_matrix[mindex[0]][mindex[1]] + distance_matrix[mindex[1]][header[0]] - distance_matrix[mindex[0]][header[0]])/2)
    print(f'new internal node {internal_count} -> node {reference_list[nh[0]]} [weight {d1}] ' )
    print(f'new internal node {internal_count} -> node {reference_list[nh[1]]} [weight {d2}] ' )  
    global output_list
    output_list.append(f'new internal node {internal_count} -> node {reference_list[nh[0]]} [weight {d1}] ' )
    output_list.append(f'new internal node {internal_count} -> node {reference_list[nh[1]]} [weight {d2}] ')
    return [d1, d2]



#main
with open("msa.pir") as file:
    contents = file.read()
seqs = contents.strip().split(">")
seqs.pop(0)
for index, seq in enumerate(seqs): 
    seqs[index] = seq[6:-1]
distance_matrix = [[0,0,0,0,0,0], 
                   [0,0,0,0,0,0], 
                   [0,0,0,0,0,0], 
                   [0,0,0,0,0,0], 
                   [0,0,0,0,0,0], 
                   [0,0,0,0,0,0]]
main_header = [x for x in range(len(distance_matrix))]
reference_list = [x for x in range(len(distance_matrix)+len(distance_matrix)-2)]

#print(reference_list)
for index, seq1 in enumerate(seqs): 
    for innerdex, seq2 in enumerate(seqs):
        #got this from stack overflow because i was curious if there was a way to get the sum of differences 
        #via list comphrehension and apparently there is an even shorter way abusing the fact that True = 1 and False = 0
        #sum(a != b for a, b in zip(seq1, seq2)) + abs(len(seq1) - len(seq2)), pretty clever IMO
        count = sum(1 for a, b in zip(seq1, seq2) if a != b) + abs(len(seq1) - len(seq2))
        distance_matrix[index][innerdex] = count
og = distance_matrix

#until make D method returns a 2x2 which will be the base case and join the two nodes
while(len(distance_matrix) != 2): 
    dPrime = makeDPrime(distance_matrix)
    output = makeD(dPrime, distance_matrix)
    distance_matrix = output[0]
    not_header = output[1]
    del reference_list[not_header[0]]
    del reference_list[not_header[1]-1]

print(f'node {reference_list[not_header[1]]}  - node {reference_list[not_header[0]]} [weight {distance_matrix[0][1]}] ' )
output_list.append(f'node {reference_list[not_header[1]]}  - node {reference_list[not_header[0]]} [weight {distance_matrix[0][1]}] ')

file1 = open('tree.dot', 'w')
file1.writelines(output_list)
file1.close()

dot = Digraph()
dot.attr(rankdir='BT', size='8,5')

# Add nodes
dot.node('Seq1', 'Seq1')
dot.node('Seq2', 'Seq2')
dot.node('Seq3', 'Seq3')
dot.node('Seq4', 'Seq4')
dot.node('Seq5', 'Seq5')
dot.node('Seq6', 'Seq6')
dot.node('m1', 'm1')
dot.node('m2', 'm2')
dot.node('m3', 'm3')
dot.node('m4', 'm4')

# Add edges
dot.edge('Seq5', 'm1', label='75')
dot.edge('Seq6', 'm1', label='82')

dot.edge('m1', 'm2', label='77')
dot.edge('Seq4', 'm2', label='57')

dot.edge('Seq1', 'm3', label='75')
dot.edge('Seq2', 'm3', label='54')

dot.edge('m2', 'm4', label = '61')
dot.edge('Seq3', 'm4', label = '111')

dot.edge('m3', 'm4', label = '61')
# Save the graph to a file
dot.render('phylogenetic_tree')
