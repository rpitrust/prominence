#! /usr/bin/python

import sys

'''
For each community, create a node-induced sub-graph using nodes in the community.
'''

# Read cluster assignment
# Each line is <node id> <cluster id>
def readClusAssigns( fn ):
    f = open( fn, 'r' );
    ret = {};
    for line in f:
        a = line.strip().split();
        if ( len(a) == 0 ):
            continue;
        assert( len(a) >= 2 );
        name = a[0];
        # a node could belong to multiple clusters
        for clusid in a[1:]:
            if name not in ret:
                ret[name] = set()
            ret[name].add( clusid )
    f.close();
    return ret;

# Create node-induced subgraph
def createSubGraphs( fn, clusAssigns ):
    f = open( fn, 'r' );

    ret = {}
    for line in f:
        a = line.strip().split();
        assert( len(a) == 3 )
        name1 = a[0]
        name2 = a[1]
        weight = a[2]
        # a node may bleong to multiple clusters
        for clus1 in clusAssigns[name1]:
            for clus2 in clusAssigns[name2]:
                if clus1 == clus2:
                    if clus1 not in ret:
                        ret[clus1] = set();
                    ret[clus1].add(line);

    f.close();
    return ret;

def writeBack( clusGraphs ):
    for clusid in clusGraphs:
        w = open( "cluster%s_edges.txt"%clusid, 'w' );

        for line in clusGraphs[clusid]:
            w.write( line );
        w.close();
    return;

def usage( prog ):
    print( "Usage: %s <cluster assigns> <network>"%prog );
    return;

def main( argv ):

    if ( len( argv ) < 3 ):
        usage( argv[0] );
        exit( -1 );

    clusAssigns = readClusAssigns( argv[1] );
    clusGraphs = createSubGraphs( argv[2], clusAssigns );
    writeBack( clusGraphs );
    return;

if __name__ == "__main__":
    main( sys.argv );
