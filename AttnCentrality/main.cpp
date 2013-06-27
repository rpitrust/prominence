/*
 * Attentive betweenness centrality (ABC) - infers importance of actors in a social network.
 *
 * Copyright (c) 2013, Xiaohui Lu, Sibel Adali and Malik Magdon-Ismail
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *	• Redistributions of source code must retain the above copyright notice, 
 *  this list of conditions and the following disclaimer.
 *	• Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Any publication resulting from the use of this work must cite the
 * following publication :
 *	S. Adali, X. Lu and M. Magdon-Ismail, "Attentive Betweenness Centrality
 *	(ABC): Considering Options and Bandwidth When Measuring Criticality".
 *	SocialCom/PASSAT 2012: 358-367.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// assertion switch
// comment out the following line if need assertion to work
// #define NDEBUG

// c library
#include <cstdlib>
#include <cassert>

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

/* read graph from file and add vertices and edges to the graph
 * parameters:
 * const std::string fn: input file name
 * Graph &g: boost graph structure
 *
 * Note: weight will be replace by 1.0 no matter what value is as the current
 * version of ABC does not support weighted network. This may be expanded in the
 * future.
 */
template<typename Graph>
const bool readGraph( const std::string fn, Graph &g ) {

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
            // unweighted graph
            put(edge_weight, g, e, 1.0);
        }
    }

    f.close();

    return true;
}

/*
 * Perform a breadth-first-search-like traversal of the graph. This process
 * computes normalization terms for parent-child(normP) and sibling(normS) levels.
 * normP = total strength of parent-child level edges for a parent node.
 * normS = total strength of sibling level edges for a node u.
 *
 * WHITE: not discovered.
 * RED: discovered, but in immediate next level.
 * GRAY: discovered, in the current processing level.
 * BLACK: processed.
 *
 * parameters:
 * Graph &g: boost graph structure
 * const typename graph_traits<Graph>::vertex_descriptor s: boost graph vertex descriptor for iterating vertices
 * std::vector<double >& normP: vector of doubles of normalization for parent-child flow and credit computation
 * std::vector<double >& normS: vector of doubles of normalization for same level flow and credit computation
 */
template <typename Graph>
void computeNorms( const Graph& g, const typename graph_traits<Graph>::vertex_descriptor s,
                                std::vector<double >& normP, std::vector<double >& normS ) {

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
                normP[u] += 1.0;
            } else {
                // v has been discovered, but has multiple parents (multiple contributors)
                if ( vertexColor[v] == RED ) {
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
 * Perform a breadth-first-search-like traversal of the graph. This process
 * assigns flow to each node and gives credit back to its source (parent and sibling)
 *
 * parameters:
 * Graph &g: boost graph structure
 * const typename graph_traits<Graph>::vertex_descriptor s: boost graph vertex descriptor for iterating vertices
 * const double alpha: damping factor to control how much flow to ditribute in the next action
 */
template<typename Graph>
const vector<double > singleSourceAttnCentrality( const Graph &g, const typename graph_traits<Graph>::vertex_descriptor s, const double alpha )
{
    // number of vertices
    const size_t V = num_vertices( g );

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

    // compute normalization terms
    std::vector<double > normP, normS;
    computeNorms( g, s, normP, normS );


    // the amount of flow
    std::vector< double > flow( V, 0.0 );
    // the amount of flow from sibling
    std::vector< double > flowS( V, 0.0 );
    // level info
    std::vector< int > level( V, -1 );

    // prepare things for BFS
    // process vertices level by level
    std::queue<vertex_descriptor > currentLevel, nextLevel;
    // stack of vertices for bottom-up operation
    std::stack< std::queue<vertex_descriptor > > S;

    flow[s] = 1.0;
    level[s] = 0;
    currentLevel.push( s );
    S.push( currentLevel );

    vertex_descriptor u;
    out_edge_iterator ei, ei_end;
    while ( !currentLevel.empty() ) {
        u = currentLevel.front();
        currentLevel.pop();

        double totalStrength = normP[u];
        
        totalStrength += normS[u];
        // add flow from sibling before propagating to
        // next level
        flow[u] += flowS[u];


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
                flow[v] += flow[u] * alpha * mystrength / totalStrength;
            }
        }

        // currentLevel has been processed
        // check if there are more vertices needed to be processed, i.e. next
        // level has elements.
        if ( currentLevel.empty() ) {
            if ( !nextLevel.empty() ) {
                S.push( nextLevel );
                currentLevel = nextLevel;
            }

            while( !nextLevel.empty() ) {
                u = nextLevel.front();
                nextLevel.pop();
                double totalStrength = normP[u] + normS[u];
                for ( tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei ) {
                    vertex_descriptor v = target(*ei, g);
                    if ( level[v] == level[u] ) {
                        assert( level[v] > 0 );
                        assert( level[u] > 0 );
                        // flow in the same level
                        double mystrength = 1.0;
                        flowS[v] += flow[u] * mystrength * alpha / totalStrength;
                    }
                }
            }
        }
    }

    std::vector<double> scoreP( V, 0.0 );
    std::vector<double> scoreS( V, 0.0 );
    while( !S.empty() ) {
        currentLevel = S.top();
        S.pop();

        std::vector<vertex_descriptor > aLevel( currentLevel.size(), vertex_descriptor() );
        size_t indx = 0;
        while ( !currentLevel.empty() ) {
            aLevel[indx++] = currentLevel.front();
            currentLevel.pop();
        }

        typename std::vector<vertex_descriptor >::iterator it;
        vertex_descriptor w, u;

        typedef typename graph_traits<Graph>::in_edge_iterator in_edge_iterator;
        in_edge_iterator ei, ei_end;

        // give credit back to actors in the same level, because getting information
        // from them
        for ( it = aLevel.begin(); it != aLevel.end(); ++it ) {
            w = *it;
            assert(flow[w] > 0.0);
            for ( tie(ei, ei_end) = in_edges(w, g); ei != ei_end; ++ei ) {
                u = source( *ei, g );
                if ( level[u] == level[w] ) {
                    double mystrength = 1.0;
                    scoreS[u] += (1.0 + scoreP[w]) * (mystrength * (flow[u] - flowS[u]) * alpha) / ((normP[u] + normS[u]) * flow[w]);
                }

            }
        }

        // update credit before giving back to parent level
        for ( it = aLevel.begin(); it != aLevel.end(); ++it ) {
            w = *it;
            scoreP[w] += scoreS[w];
        }

        // give credit back to actors in the parent level
        for ( it = aLevel.begin(); it != aLevel.end(); ++it ) {
            w = *it;
            assert(flow[w] > 0.0);
            for ( tie(ei, ei_end) = in_edges(w, g); ei != ei_end; ++ei ) {
                    u = source( *ei, g );
                    if ( level[u] + 1 == level[w] ) {
                        double mystrength = 1.0;
                        double totalStrength = normP[u];
                        totalStrength += normS[u];
                        scoreP[u] += (1.0 + scoreP[w]) * (flow[u] * mystrength * alpha) / (totalStrength * flow[w]) ;
                    }
            }
        }
    }

    scoreP[s] = 0.0;

    return scoreP;
}

/*
 * Main entry of attentative betweenness centrality.
 *
 * Perform single cource ABC for each node and normalized that value.
 *
 * parameters:
 * Graph &g: boost graph structure
 * const double alpha: damping factor to control how much flow to ditribute in the next action
 */
template <typename Graph>
const vector<double > attnCentrality( Graph &g, const double alpha ) {
    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    vertex_iterator vi, vi_end;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

    const size_t V = num_vertices(g);

    vector<double> ret(V, 0);

    for ( tie(vi, vi_end) = vertices(g); vi != vi_end; vi++ ) {
        vector<double> bu = singleSourceAttnCentrality( g, *vi, alpha );
        for (size_t i = 0; i < V; i++)
            ret[i] += bu[i];
    }
    
    // normalization
    const double normalization = (V - 1) * (V - 2) * 0.5;
    for ( vector<double>::iterator it = ret.begin(); it != ret.end(); ++it ) {
        *it = *it / normalization;
    }
    return ret;
}

/*
 * Print results.
 *
 * parameters:
 * Graph &g: boost graph structure
 * const std::vector<double> &centrality: centrality resluts
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
 * Print usage message.
 *
 * parameters:
 * const std::string &prog: the executable name
 */
void usage( const std::string &prog ) {
    std::cout << "Usage: " << prog << " <graph> [alpha]" << std::endl;
    std::cout << "Each line of the graph file must be <source_node> <target_node> <edge_strength>" << std::endl;
    std::cout << "0 < alpha <= 1.0, alpha = 1.0 by default." << std::endl;
    return;
}

/*
 * Main entry of the program.
 */
int main(int argc, char** argv) {
    if ( argc < 2 ) {
        usage( argv[0] );
        std::exit( -1 );
    }

    const std::string fn(argv[1]);
    double alpha = 1.0;

    if ( argc > 2 ) {
        alpha = atof( argv[2] );
    }

    if ( alpha < 0.0 || alpha > 1.0 ) {
        std::cerr << "alpha must be in range (0, 1.0]";
        std::exit( -1 );
    }

    // edge properties
    typedef property< edge_weight_t, double > EdgeProperties;
    // vertex properties
    typedef property< vertex_name_t, std::string > VertexProperties;

    // undirected simple graph
    typedef adjacency_list<listS, vecS, undirectedS,
                         VertexProperties, EdgeProperties> Graph;

    Graph g;

    bool readSucc = readGraph( fn, g );
    if ( !readSucc ) {
        std::cerr << "Something wrong while reading " << fn << "." << std::endl;
        std::exit( -1 );
    }
    
    std::vector<double> centrality = attnCentrality( g, alpha );

    printResults( g, centrality );

    return 0;
}
