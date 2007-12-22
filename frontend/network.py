def read_undirected_network(fileName):

    # Read the raw data
    # -------------------------------------------------------------------------
    file = open(fileName)
    data = file.readlines()
    file.close()

    # Read nodes and links
    # -------------------------------------------------------------------------
    linkList = []
    nodeList = []
    for aLine in data:
        try:
            [node1, node2] = aLine.strip().split()
        except:
            return 0, 0, 0 # Not a valid network file!

        if not node1 in nodeList:
            nodeList.append(node1)
        if not node2 in nodeList:
            nodeList.append(node2)
        if (not [node1, node2] in linkList) and \
           (not[node2, node1] in linkList) :
            linkList.append([node1, node2])
        
    # Done: return 1, #nodes, #links
    # -------------------------------------------------------------------------
    return 1, len(nodeList), len(linkList)
