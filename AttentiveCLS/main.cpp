/*
 * Attentive closeness centrality (ACC) - infers individuals' efficiency of
 * spreading information.
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
 * Authors: Dr. Sibel Adali (sibel@cs.rpi.edu), Dr. Malik Magdon-Ismail (magdon@cs.rpi.edu), and Xiaohui Lu (lux3@cs.rpi.edu)
 *
 * Created on March 6, 2013 10:39 PM
 */

// c library
#include <cstdlib>
#include <cassert>
#include <cmath>

// c++ library
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <stack>

// bost library
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>


using namespace boost;
using namespace std;
// vertex color
typedef enum {WHITE = 0, RED, GRAY, BLACK} color;

/*
 * Read graph from file and add vertices and edges to the boost graph structure.
 *
 * Parameters:
 * const std::string fn: input file name
 * Graph &g: boost graph structure
 * const bool weighted: weighted or unweighted network
 *   - if weighted, edge weights are read from the input file
 *   - if unweighted, set all edge weights to 1.0
 */
template<typename Graph>
const bool readGraph( const std::string fn, Graph &g, const bool weighted = true ) {

    // vertex descriptor
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    // edge descriptor
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    // To store vertex name
    typedef typename property_map<Graph, vertex_name_t>::type VertexNameMap;
    VertexNameMap nodename = get(vertex_name, g);
    // Vertex name : vertex_descriptor map
    // To enforce unique name:vertex_desc map
    typedef typename std::map<std::string, vertex_descriptor> NameVertexDescMap;
    NameVertexDescMap nodes;

    std::ifstream f( fn.c_str() );

    if ( !f.is_open() ) {
        std::cerr << "Can not open file " << fn << "." << std::endl;
        return false;
    }

    vertex_descriptor u, v;
    edge_descriptor e;
    std::string node1, node2;
    double w;
    typename NameVertexDescMap::iterator pos;
    bool inserted;

    while ( f >> node1 >> node2 >> w ) {
        assert( w > 0.0 );
        // add first node
        tie(pos,inserted) = nodes.insert(std::make_pair(node1, vertex_descriptor()));
        if(inserted) { // new vertex
            u = add_vertex(g);
            nodename[u] = node1;
            pos->second = u;
        } else // vertex already exist
            u = pos->second;
        // add second node
        tie(pos,inserted) = nodes.insert(std::make_pair(node2, vertex_descriptor()));
        if(inserted) { // new vertex
            v = add_vertex(g);
            nodename[v] = node2;
            pos->second = v;
        } else // vertex already exist
            v = pos->second;
        // add edge
        tie(e,inserted) = add_edge(u,v,g);
        if(inserted) { // edge weight
            if ( weighted )
                // weighted graph, here the weight is the strengh of the edge
                // not the distance of nodes
                put(edge_weight, g, e, w);
            else
                // unweighted graph
                put(edge_weight, g, e, 1.0);
        }
    }

    f.close();

    return true;
}

/*
 * Compute normalization term, i.e. total child and same level edge weights in a
 * BFS traversal. This process computes normalization terms for parent-child(normP)
 * and sibling(normS) levels.
 * normP = total strength of parent-child level edges for a parent node.
 * normS = total strength of sibling level edges for a node u.
 *
 * WHITE: not discovered.
 * RED: discovered, but in immediate next level.
 * GRAY: discovered, in the current processing level.
 * BLACK: processed.
 *
 * Parameters:
 * Graph &g: boost graph structure
 * const typename graph_traits<Graph>::vertex_descriptor s: boost graph vertex descriptor for iterating vertices
 * std::vector<double >& normP: vector of doubles of normalization for parent-child flow and credit computation
 * std::vector<double >& normS: vector of doubles of normalization for same level flow and credit computation
 */
template <typename Graph>
void computeNorms( const Graph& g, const typename graph_traits<Graph>::vertex_descriptor s,
                                std::vector<double >& normP, std::vector<double >& normS, const bool weighted = true ) {

    // number of vertices in the graph
    const size_t V = boost::num_vertices(g);
    normP.resize(V, 0.0);
    normS.resize(V, 0.0);

    // prepare things for breadth-first traverse
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    // process vertices level by level
    std::queue<vertex_descriptor > currentLevel, nextLevel;
    // color all vertices to be white
    std::vector<color> vertexColor( V, WHITE );

    // init source node
    vertexColor[s] = GRAY;
    currentLevel.push( s );

    vertex_descriptor u;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    out_edge_iterator ei, ei_end;
    while ( !currentLevel.empty() ) {
        u = currentLevel.front();
        currentLevel.pop();
        for ( tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei ) {
            // scan vertices pointed to by u
            // v is a vertex pointed to by u
            vertex_descriptor v = target(*ei, g);
            if ( vertexColor[v] == WHITE ) {
                // v has not been discovered, white
                // put v into the immediate next level, red
                vertexColor[v] = RED;
                nextLevel.push(v);
                if ( weighted )
                    normP[u] += get(edge_weight, g, *ei);
                else
                    normP[u] += 1.0;
            } else {
                // v has been discovered, but has multiple parents (multiple contributors)
                if ( vertexColor[v] == RED ) {
                    if ( weighted )
                        normP[u] += get(edge_weight, g, *ei);
                    else
                        normP[u] += 1.0;
                }
            }

        }
        // u has been processed, black
        vertexColor[u] = BLACK;

        // currentLevel has been processed
        // check if there are more vertices needed to be processed, i.e. next
        // level has elements.
        if ( currentLevel.empty() ) {
            while( !nextLevel.empty() ) {
                u = nextLevel.front();
                nextLevel.pop();
                for ( tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei ) {
                    vertex_descriptor v = target(*ei, g);
                    if ( vertexColor[v] == RED || vertexColor[v] == GRAY ) {
                        if ( weighted )
                            normS[u] += get(edge_weight, g, *ei);
                        else
                            normS[u] += 1.0;
                    }
                }
                // move vertex u to current level
                currentLevel.push( u );
                // change color of u
                vertexColor[u] = GRAY;
            }
        }
    }

    return;
}

/*
 * Compute flow distances for a specific vertex.
 *
 * Parameters:
 * Graph &g: boost graph structure
 * const typename graph_traits<Graph>::vertex_descriptor s: boost graph vertex descriptor for iterating vertices
 * const double alpha: damping factor to control how much flow to ditribute in the next action
 * const double epsilon: a factor controls how much information flow to virtual path
 * const bool forwardOnly: internal test only, always fase
 * const bool weighted: weighted or unweighted network
 */
template<typename Graph>
const vector<double > singleSourceFlowDistance( const Graph &g, const typename graph_traits<Graph>::vertex_descriptor s,
                                        const double alpha, const double epsilon, const bool forwardOnly, const bool weighted )
{
    // number of vertices
    const size_t V = num_vertices( g );

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

    // compute normalization terms
    std::vector<double > normP, normS;
    computeNorms( g, s, normP, normS, weighted );


    // the amount of flow
    std::vector< double > flow( V, 0.0 );
    // the amount of flow from sibling
    std::vector< double > flowS( V, 0.0 );


    // final score
    std::vector< double > score( V, 0.0 );
    // score from siblings
    std::vector< double > scoreS( V, 0.0 );

    // level info
    std::vector< int > level( V, -1 );

    // prepare things for BFS
    // process vertices level by level
    std::queue<vertex_descriptor > currentLevel, nextLevel;

    // epsilon information goes to an imaginary route
    // the rest goes to classic
    flow[s] = 1.0 - epsilon;
    level[s] = 0;
    currentLevel.push( s );

    vertex_descriptor u;
    out_edge_iterator ei, ei_end;
    while ( !currentLevel.empty() ) {
        u = currentLevel.front();
        currentLevel.pop();

        double totalStrength = normP[u];
        if ( !forwardOnly ) {
            totalStrength += normS[u];
            // add flow from sibling before propagating to
            // next level
            flow[u] += flowS[u];

            score[u] += scoreS[u];
        }

        for ( tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei ) {
            // scan vertices pointed to by u
            // v is a vertex pointed to by u
            vertex_descriptor v = target(*ei, g);
            if ( level[v] < 0 ) {
                // put v into the immediate next level
                nextLevel.push(v);
                level[v] = level[u] + 1;
                assert( level[v] > 0 );
            }

            if ( level[v] == level[u] + 1 ) {
                // flow propagation from parent to child level
                double mystrength = 1.0;
                if ( weighted )
                    mystrength = get(edge_weight, g, *ei);
                const double norm = alpha * mystrength / totalStrength;
                flow[v] += flow[u] * norm;

		// we take 1/mystrength as the distance
		// therefore, there is no distance and mystrength in the second
		// part of the summation
                score[v] += score[u] * norm + flow[u] * alpha / totalStrength;
            }
        }

        // currentLevel has been processed
        // check if there are more vertices needed to be processed, i.e. next
        // level has elements.
        if ( currentLevel.empty() ) {
            while( !nextLevel.empty() ) {
                u = nextLevel.front();
                currentLevel.push( u );
                nextLevel.pop();
                if ( !forwardOnly ) {
                    double totalStrength = normP[u] + normS[u];
                    for ( tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei ) {
                        vertex_descriptor v = target(*ei, g);
                        if ( level[v] == level[u] ) {
                            assert( level[v] > 0 );
                            assert( level[u] > 0 );
                            // flow in the same level
                            double mystrength = 1.0;
                            if ( weighted )
                                mystrength = get(edge_weight, g, *ei);
                            const double norm = mystrength * alpha / totalStrength;
                            flowS[v] += flow[u] * norm;

			    // we take 1/mystrength as the distance
			    // therefore, there is no distance and mystrength in the second
			    // part of the summation
                            scoreS[v] += score[u] * norm + flow[u] * alpha / totalStrength;
                        }
                    }
                }
            }
        }
    }

    // the length of virtual path from source node to others
    // assuming the length of the path is 6 (six degree separation)
    const double epsilonPath = 6.0;
    const double f = epsilon * pow(alpha, epsilonPath);

    for ( size_t i = 0; i < score.size(); i++ ) {
        score[i] = (score[i] + f * epsilonPath) / (flow[i] + f);
    }
    score[s] = 0.0;
    
    return score;
}

/*
 * Print centrality score.
 *
 * Parameters:
 * Graph& g: boost graph structure
 * const std::vector<double> &centrality: centrality score fore each vertex
 */
template<typename Graph>
void printResults( Graph& g, const std::vector<double> &centrality ) {
    std::string fn( "results.txt" );
    std::ofstream w( fn.c_str() );

    if ( !w.is_open() ) {
        std::cerr << "Can not open file " << fn << "." << std::endl;
        return;
    }

    // To access vertex name
    typedef typename property_map<Graph, vertex_name_t>::type VertexNameMap;
    const VertexNameMap nodename = get(vertex_name, g);
    // To iterate over vertices
    typedef typename graph_traits< Graph >::vertex_iterator vertex_iterator;
    vertex_iterator vi, ve;
    for ( tie(vi, ve) = vertices(g); vi != ve; vi++ ) {
        w << nodename[*vi] << "\t" << centrality[*vi] << std::endl;
    }
    
    w.close();

    return;

}

/*
 * Main entry of attentative Closeness centrality. It computes flow distance for
 * each vertex.
 *
 * parameters:
 * Graph &g: boost graph structure
 * const double alpha: damping factor to control how much flow to ditribute in the next action
 * const double epsilon: a factor controls how much inoformation to virtual path
 * const bool forwardOnly: internal test only, do not use
 * const bool weighted: weighted or unweighted network
 */
template <typename Graph>
const vector<double > flowDistance( Graph &g, const double alpha, const double epsilon, const bool forwardOnly, const bool weighted ) {
    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    vertex_iterator vi, vi_end;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

    const size_t V = num_vertices(g);

    vector<double> ret(V, 0);

    for ( tie(vi, vi_end) = vertices(g); vi != vi_end; vi++ ) {
        vector<double> bu = singleSourceFlowDistance( g, *vi, alpha, epsilon, forwardOnly, weighted );

        for (size_t i = 0; i < V; i++)
            ret[*vi] += bu[i];
    }

    // normalization
    const double normalization = V-1.0;
    for ( vector<double>::iterator it = ret.begin(); it != ret.end(); ++it ) {
        *it = normalization / (*it);
    }
    return ret;
}


/*
 * Print usage information.
 */
void usage( const std::string &prog ) {
    std::cout << "Usage: " << prog << " <graph> [U|W] [alpha] [epsilon] [F]" << std::endl;
    std::cout << "U/W indicates a unweighted or weighted graph. Weighted is the default value." << std::endl;
    std::cout << "U = unweighted graph, other value = weighted graph." << std::endl;
    std::cout << "Each line of the graph file must be <source_node> <target_node> <edge_strength>" << std::endl;
    std::cout << "No matter weighted of unweighted graph." << std::endl;
    std::cout << "0 < alpha <= 1.0, alpha = 1.0 by default." << std::endl;
    std::cout << "0 <= epsilon < 1.0, epsilon = 0.05 by default." << std::endl;
    std::cout << "F = forward only, i.e. no flow among nodes in the same level of a BFS tree." << std::endl;
    std::cout << "Forward only is false by default." << std::endl;
    return;
}

/*
 * main entry of the program
 */
int main(int argc, char** argv) {
    if ( argc < 2 ) {
        usage( argv[0] );
        std::exit( -1 );
    }

    const std::string fn(argv[1]);
    double alpha = 1.0;
    double epsilon = 0.05;

    bool weighted = true;
    if ( argc > 2 ) {
        if ( 'U' == argv[2][0] )
            weighted = false;
    }

    if ( argc > 3 ) {
        alpha = atof( argv[3] );
    }
    assert( alpha > 0.0 && alpha <= 1.0 );

    if ( argc > 4 ) {
        epsilon = atof( argv[4] );
    }
    assert( epsilon >= 0.0 && epsilon < 1.0 );

    bool forwardOnly = false;
    if ( argc > 5 ) {
        if ( 'F' == argv[5][0] )
            forwardOnly = true;
    }

    // edge properties
    typedef property< edge_weight_t, double > EdgeProperties;
    // vertex properties
    typedef property< vertex_name_t, std::string > VertexProperties;

    // undirected simple graph
    typedef adjacency_list<listS, vecS, undirectedS,
                         VertexProperties, EdgeProperties> Graph;

    Graph g;

    bool readSucc = readGraph( fn, g, weighted );
    if ( !readSucc ) {
        std::cerr << "Something wrong while reading " << fn << "." << std::endl;
        return ( -1 );
    }

    std::vector<double> importance = flowDistance( g, alpha, epsilon, forwardOnly, weighted );

    printResults( g, importance );

    return 0;
}
