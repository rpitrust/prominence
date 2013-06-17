/*
 *
 * Given a connected network and a node partition (node communities), compute
 * community distances between two communities A and B by:
 * 1. construct (A union B)-node induced sub-graph.
 * 2. randomly sample k% nodes from community A and compute shortest distances from
 * sampled nodes to nodes in community B using the induced graph.
 * 3. randomly sample k% nodes from community B and compute shortest distances from
 * sampled nodes to nodes in community A using the induced graph.
 * 4. the average of these distances is defined as the distance between A and B.
 * 5. in case that we can not sample a node from a community (small community), increase
 * sample probability by k%.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * See http://www.gnu.org/copyleft/gpl.txt for more details.
 *
 * File:   main.cpp
 * Contributors: Dr. Sibel Adali (sibel@cs.rpi.edu), Dr. Malik Magdon-Ismail (magdon@cs.rpi.edu)
 * Code Author: Xiaohui Lu (lux3@cs.rpi.edu)
 *
 */

#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

typedef unsigned int uint;

// edge properties
typedef property< edge_weight_t, double > EdgeProperties;
// vertex properties
typedef property< vertex_name_t, std::string > VertexProperties;
// undirected simple graph
typedef adjacency_list<listS, vecS, undirectedS,
                     VertexProperties, EdgeProperties> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;

// Vertex name : vertex_descriptor map
// To enforce unique name:vertex_desc map
typedef typename std::map<std::string, Vertex> NameVertexMap;

/*
 * Read clus information and store them into unordered_map where key is clusid and
 * value is the set of nodes in the cluster.
 *
 * Parameters:
 * std::istream& in: input stream
 * std::unordered_map< uint, std::unordered_set<Vertex> > &clus: hash map for storing cluster informaiton
 * const NameVertexMap &vertex: vertext name
 */
void readClusInfo( std::istream& in, std::unordered_map< uint, std::unordered_set<Vertex> > &clus, const NameVertexMap &vertex ) {

    std::string actor;
    uint clusid;

    int numnodes = 0;
    while ( in >> actor >> clusid ) {
        NameVertexMap::const_iterator it = vertex.find( actor );
        if ( it != vertex.end() ) {
            clus[ clusid ].insert( it->second );
            numnodes++;
        } else {
            std::cerr << "Actor is not in the network " << actor << "." << std::endl;
            std::exit( -1 );
        }
    }

    //std::cout << "Total number clusters: " << clus.size() << std::endl;
    //std::cout << "Total number actors in cluster assigns: " << numnodes << std::endl;
    return;
}

/*
 * Read cluster assignment, i.e. for each node get its cluster id. Cluster id is
 * supposed to be an non-negative integer.
 *
 * Parameters:
 * std::istream& in: input stream
 * std::unordered_map< Vertex, uint > &clusAssigns: hash map for storing cluster assignments
 * const NameVertexMap &vertex: vertex name
 */
void readClusAssigns( std::istream& in, std::unordered_map< Vertex, uint > &clusAssigns, const NameVertexMap &vertex ) {
    std::string actor;
    uint clusid;

    while ( in >> actor >> clusid ) {
        NameVertexMap::const_iterator it = vertex.find( actor );
        if ( it != vertex.end() ) {
            clusAssigns[ it->second ] = clusid;
        } else {
            std::cerr << "Actor is not in the network " << actor << "." << std::endl;
            std::exit( -1 );
        }
    }

    return;
}

/*
 * Read graph into boost graph structure. Each line of the input must be
 * <node1> <node2> <weight>.
 * The weight is supposed to be positive. When it is a unweighted graph, the weight
 * in the file is ignored and replaced by 1.0.
 *
 * Parameters:
 * std::istream& in: input stream
 * Graph &g: boost graph structure
 * NameVertexMap &nodes: vertex name map
 * const bool weighted: weighted or unweighted network
 */
template <typename Graph>
void readGraph( std::istream& in, Graph &g, NameVertexMap &nodes, const bool weighted ) {

    // To store vertex name
    typedef typename property_map<Graph, vertex_name_t>::type VertexNameMap;
    VertexNameMap nodename = get(vertex_name, g);

    Vertex u, v;
    Edge e;
    std::string node1, node2;
    double w;
    NameVertexMap::iterator pos;
    bool inserted;

    while ( in >> node1 >> node2 >> w ) {
        // add first node
        tie(pos,inserted) = nodes.insert(std::make_pair(node1, Vertex()));
        if(inserted) { // new vertex
            u = add_vertex(g);
            nodename[u] = node1;
            pos->second = u;
        } else // vertex already exist
            u = pos->second;
        // add second node
        tie(pos,inserted) = nodes.insert(std::make_pair(node2, Vertex()));
        if(inserted) { // new vertex
            v = add_vertex(g);
            nodename[v] = node2;
            pos->second = v;
        } else // vertex already exist
            v = pos->second;
        // add edge
        tie(e,inserted) = add_edge(u,v,g);
        if( inserted ) { // edge weight
            if ( weighted )
                put(edge_weight, g, e, w);
            else
                // unweighted graph here
                put(edge_weight, g, e, 1.0);
        }
    }

    //std::cout << "Total number of nodes in graph: " << num_vertices( g ) << std::endl;
    //std::cout << "Total number of edges in graph: " << num_edges( g ) << std::endl;
    return;
}

/*
 * Find all connected communities by go through edges in the global graph. If end-points
 * of an edge are in different communities, then the two communities are connected, i.e.
 * two communities are connected if they share one edge.
 *
 * Parameters:
 * const Graph &g: boost graph structure
 * const std::unordered_map< Vertex, uint > &clusAssigns: cluster assignments
 * std::unordered_map< uint, std::unordered_set< uint > > &ret: hash map for connected communities
 */
template <typename Graph>
void connComm( const Graph &g, const std::unordered_map< Vertex, uint > &clusAssigns, std::unordered_map< uint, std::unordered_set< uint > > &ret ) {
    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

    for (edge_iterator e = edges(g).first; e != edges(g).second; ++e) {
        Vertex u = source(*e, g);
        Vertex v = target(*e, g);

        uint i = (*clusAssigns.find(u)).second;
        uint j = (*clusAssigns.find(v)).second;
        if ( i < j ) {
            ret[i].insert( j );
        } else if ( i > j ) {
            ret[j].insert( i );
        } else {
            // the two nodes are in the same community
            continue;
        }
    }
    return;
}

template <typename Graph>
void printSet( Graph &g, std::set< Vertex > &myset ) {

    // To access vertex name
    typedef typename property_map<Graph, vertex_name_t>::type VertexNameMap;
    const VertexNameMap nodename = get(vertex_name, g);
    // To iterate over vertices
    for ( std::set< Vertex >::iterator it = myset.begin(); it != myset.end(); ++it ) {
        std::cout << nodename[*it] << std::endl;
    }
    return;
}


template <typename Graph>
void printGraph( Graph &g ) {

    std::ofstream w("subgraph.txt");

    std::cout << "------Number of nodes in the subgraph " << num_vertices(g) << std::endl;
    std::cout << "------Number of edges in the subgraph " << num_edges(g) << std::endl;
    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
    for (edge_iterator e = edges(g).first; e != edges(g).second; ++e) {
        Vertex u = source(*e, g);
        Vertex v = target(*e, g);

        w << u << "\t" << v << "\t" << get( edge_weight, g, *e ) << std::endl;

    }

    w.close();
    return;
}
/*
 * Compute community distance of A and B.
 *
 * Parameters:
 * const Graph &g: boost graph structure
 * const std::unordered_set< Vertex > &setA: vertices in community A
 * const std::unordered_set< Vertex > &setB: vertices in community B
 * const double p: sample probability
 * double &maxdist: maximum distance
 * double &mindist: minimum distance
 */
template <typename Graph>
double commDist( const Graph &g, const std::unordered_set< Vertex > &setA, const std::unordered_set< Vertex > &setB, const double p, double &maxdist, double &mindist ) {
    // map global vertex descriptor to local ones
    // it is essential for constructing a subgraph using descriptors from the global graph g
    std::unordered_map< Vertex, Vertex > g2l;
    std::unordered_map< Vertex, Vertex >::iterator pos;
    bool inserted;
    Vertex v = 0;
    for ( std::unordered_set< Vertex >::const_iterator it = setA.begin(); it != setA.end(); it++ ) {
        tie(pos,inserted) = g2l.insert(std::make_pair(*it, Vertex()));
        if ( inserted ) {
            // vertex is not in the map
            pos->second = v;
            v++;
        }
    }
    for ( std::unordered_set< Vertex >::const_iterator it = setB.begin(); it != setB.end(); it++ ) {
        tie(pos,inserted) = g2l.insert(std::make_pair(*it, Vertex()));
        if ( inserted ) {
            // vertex is not in the map
            pos->second = v;
            v++;
        }
    }

    // subgraph induced by nodes
    Graph subg( g2l.size() );
    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
    for (edge_iterator e = edges(g).first; e != edges(g).second; ++e) {
        Vertex u = source(*e, g);
        Vertex v = target(*e, g);

        if ( g2l.find( u ) == g2l.end() || g2l.find( v ) == g2l.end() ) {
            // the two end-points both must be in the subgraph
            // otherwise skip
            continue;
        }

        // create subgraph
        Edge sube;
        tie( sube,inserted ) = add_edge( g2l[u], g2l[v], subg );
        if ( inserted )
            put( edge_weight, subg, sube, get(edge_weight, g, *e) );
    }

    int numComponents = 0;
    {
        std::vector<int> component( num_vertices(subg), 0 );
        numComponents = connected_components( subg, &component[0] );
    }

    // check for number of components
    if ( numComponents > 1 ) {
        std::cerr << "Number of components in the new subgraph: " << numComponents << std::endl;
        return 0.0;
    }

    // randomly sample vertices from set A, B
    std::unordered_set< Vertex > randSetA, randSetB;
    // initialize random seed
    std::srand ( time(NULL) );
    // generate random vertices and store in randSetA
    double prob = p;
    while ( randSetA.size() < 1 ) {
        for ( std::unordered_set< Vertex >::const_iterator it = setA.begin(); it != setA.end(); it++ ) {
            double r = (float)rand()/(float)RAND_MAX;
            if ( r <= prob ) {
                randSetA.insert( g2l[*it] );
            }
        }
        // the set is too small, we will not be able to sample enough acotrs
        // set p to a larger nubmer
        if ( (uint)(setA.size() * prob) < 1 )
            prob += p;
    }


    prob = p;
    while ( randSetB.size() < 1 ) {
        // generate random vertices and store in randSetB
        for ( std::unordered_set< Vertex >::const_iterator it = setB.begin(); it != setB.end(); it++ ) {
            double r = (float)rand()/(float)RAND_MAX;
            if ( r <= prob ) {
                randSetB.insert( g2l[*it] );
            }
        }
        // the set is too small, we will not be able to sample enough acotrs
        // set p to a larger nubmer
        if ( (uint)(setB.size() * prob) < 1 )
            prob += p;
    }

    // the random sets must contain at least one actor
    if ( randSetA.size() == 0 || randSetB.size() == 0 ) {
        std::cerr << "I can not sample enough actors from input sets." << std::endl;
        std::exit( -1 );
    }


    double ret = 0.0;
    // single source shortest distance
    uint V = num_vertices( subg );
    // compute distance from random nodes to others in the subgraph
    for ( std::unordered_set< Vertex >::iterator it = randSetA.begin(); it != randSetA.end(); it++ ) {
        // To store parents
        std::vector< Vertex > parents(V);
        // To store distances
        std::vector< double > distances(V, 0.0);
        // compute average distance for these random nodes.
        dijkstra_shortest_paths(subg, *it, boost::predecessor_map(&parents[0]).distance_map(&distances[0]));

        for ( std::unordered_set< Vertex >::const_iterator j = setB.begin(); j != setB.end(); j++ ) {
            ret += distances[ g2l[*j] ];
            if (distances[ g2l[*j] ] > maxdist) {
                maxdist = distances[ g2l[*j] ];
            }

            if (distances[ g2l[*j] ] < mindist) {
                mindist = distances[ g2l[*j] ];
            }
        }
    }

    // compute distance from random nodes to others in the subgraph
    for ( std::unordered_set< Vertex >::iterator it = randSetB.begin(); it != randSetB.end(); it++ ) {
        // To store parents
        std::vector< Vertex > parents(V);
        // To store distances
        std::vector< double > distances(V, 0.0);
        // compute average distance for these random nodes.
        dijkstra_shortest_paths(subg, *it, boost::predecessor_map(&parents[0]).distance_map(&distances[0]));

        for ( std::unordered_set< Vertex >::const_iterator j = setA.begin(); j != setA.end(); j++ ) {
            ret += distances[ g2l[*j] ];
            if (distances[ g2l[*j] ] > maxdist) {
                maxdist = distances[ g2l[*j] ];
            }

            if (distances[ g2l[*j] ] < mindist) {
                mindist = distances[ g2l[*j] ];
            }
        }
    }

    // normalization
    ret = (double)ret / (randSetA.size() * setB.size() + randSetB.size() * setA.size() );

    return ret;
}

/*
 * Save community distance results. Each line is a pair of community ids followed
 * by distance between them.
 *
 * Parameters:
 * std::ostream& out: output stream
 * std::unordered_map<uint, std::unordered_map<uint, double> >& results: community distance
 */
void writeResults( std::ostream& out, std::unordered_map<uint, std::unordered_map<uint, double> >& results ) {
    std::unordered_map< uint, std::unordered_map<uint, double> >::iterator outer;
    std::unordered_map< uint, double >::iterator inner;
    for ( outer = results.begin(); outer != results.end(); outer++ ) {
        for ( inner = (outer->second).begin(); inner != (outer->second).end(); inner++ ) {
            if ( 0.0 == inner->second )
                // special value to indicate disconnected communities
                continue;
            out << outer->first << "\t" << inner->first << "\t" << inner->second << std::endl;
        }
    }

    return;
}

/*
 * Print usage information
 */
void usage( std::string prog ) {
    std::cerr << "Usage: " << prog << " <network> <cluster assigns> [p] [weighted?]." << std::endl;
    std::cerr << "\tp = 0.05(default) or other values in (0, 1]." << std::endl;
    std::cerr << "\tweighted = 1: weighted graph (default). 0 = unweighted." << std::endl;
    return;
}

/*
 * Main entry of the program
 */
int main(int argc, char** argv) {

    if ( argc < 3 ) {
        usage( argv[0] );
        std::exit( -1 );
    }

    Graph g;
    const std::string network( argv[1] );
    const std::string clusters( argv[2] );
    NameVertexMap nodes;

    double p = 0.05;
    if ( argc > 3 ) {
        p = std::atof( argv[3] );

        if ( p <= 0.0 || p > 1 ) {
            std::cerr << "Probability must be in range (0, 1]: " << p << std::endl;
            std::cerr << "I'll use the default value 0.05." << std::endl;
            p = 0.05;
        }
    }

    bool weighted = true;
    if ( argc > 4 ) {
        if ( 0 == std::atoi( argv[4] ) )
            weighted = false;
    }

    {   // read network into boost graph structure
        std::ifstream in( network.c_str() );

        if ( in.is_open() ) {
            readGraph( in, g, nodes, weighted );
            in.close();
        } else {
            std::cerr << "Can not read network from file " << network << "." << std::endl;
            std::exit( -1 );
        }
    }

    std::unordered_map< uint, std::unordered_set<Vertex> > clus;
    {   // read cluster assigns
        std::ifstream in( clusters.c_str() );

        if ( in.is_open() ) {
            readClusInfo( in, clus, nodes );
            in.close();
        } else {
            std::cerr << "Can not read clusters from file " << clusters << "." << std::endl;
            std::exit( -1 );
        }
    }

    std::unordered_map< uint, std::unordered_set< uint > > connectedComm;
    {
        std::ifstream in( clusters.c_str() );
        std::unordered_map< Vertex, uint > clusAssigns;
        if ( in.is_open() ) {
            readClusAssigns( in, clusAssigns, nodes );
            in.close();
        } else {
            std::cerr << "Can not read clusters from file " << clusters << "." << std::endl;
            std::exit( -1 );
        }

        in.close();
        connComm( g, clusAssigns, connectedComm );
    }

    // compute distance
    std::unordered_map< uint, std::unordered_map<uint, double> > meanDistance;
    //std::unordered_map< uint, std::unordered_map<uint, double> > maxDistance;
    //std::unordered_map< uint, std::unordered_map<uint, double> > minDistance;
    for ( std::unordered_map< uint, std::unordered_set< uint > >::iterator i = connectedComm.begin(); i != connectedComm.end(); i++ ) {
        for ( std::unordered_set< uint >::iterator j = (i->second).begin(); j != (i->second).end(); j++ ) {
            double mindist = std::numeric_limits<double>::max( );
            double maxdist = std::numeric_limits<double>::min( );
            meanDistance[i->first][*j] = commDist( g, clus[i->first], clus[*j], p, maxdist, mindist );
            //minDistance[i->first][*j] = mindist;
            //maxDistance[i->first][*j] = maxdist;
        }
    }

    {
        // output results...
        std::string ofn( "clus_distances_mean.txt" );
        std::ofstream out( ofn.c_str() );

        if ( out.is_open() ) {
            writeResults( out, meanDistance );
            out.close();
        } else {
            std::cerr << "Can not write to file " << ofn << "." << std::endl;
            writeResults( std::cout, meanDistance );
        }

        out.close();
    }

    return 0;
}

