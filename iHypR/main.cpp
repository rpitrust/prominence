/*
 * iterative hyperedge ranking (iHypR) - infers prominence of actors in heterogeneous networks by
 * incoporating hyperedge structure.
 *
 * Copyright (c) 2013, Xiaohui Lu, Sibel Adali and Malik Magdon-Ismail
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *	• Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimer.
 *	• Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Any publication resulting from the use of this work must cite the
 * following publication :
 *	S. Adali, Xiaohui Lu and Malik Magdon-Ismail, "iHypR: Prominence Ranking
 *	in Networks of Collaborations with Hyperedges", accepted to appear in
 *	Transactions on Knowledge Discovery from Data, 2013.
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
 *
 */

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <fstream>
#include <iostream>
#include <vector>

#include <unistd.h>
#include <stdint.h>
using namespace std;

#include "basic.h"
#include "quicksort.h"

// maximum number of iterations
const int MAXITERATIONS	= 1000;

// output message when new operator fails
void no_memory () {
    cout << "Failed to allocate memory!\n";
    exit (1);
}


/*
 * read file into corresponding data structure
 *
 * parameters:
 * const std::string &fn: input file name
 * Entity **actor, Entity **object, Entity **hyperedge: pointers to actor, object, hyperedge arrays
 * const bool ah: true => computing actor value based on actor-hyperedge, otherwise computing actor
 * value based on actor-object-hyperedge
 *
 * Note: input file must be in format of:
 * 1. the first line must be the max ids for actors, objects, and hyperedges respectively
 * 2. each line after the first one is actor_id, object_id, hyperedge_id
 * 3. all ids are integer only and must be consecutive start from 1
 * 
 */
const vector<uint32_t > read_file( const string &fn, Entity **actor, Entity **object, Entity **hyperedge, const bool ah ) {
    uint32_t num_actor = 0, num_object = 0, num_hyperedge = 0;
    vector<uint32_t > ret;

    ifstream file( fn.c_str() );

    if ( file.is_open() ) {
        uint32_t id1 = 0, id2 = 0, id3 = 0;

        file >> num_actor >> num_object >> num_hyperedge;

        if ( num_actor > 0 )
            *actor = new Entity[ num_actor + 1 ];
        else {
            cerr << "Error found in input file " << fn << endl;
            abort();
        }

        if ( false == ah ) {
            if ( num_object > 0 )
                *object = new Entity[ num_object + 1 ];
            else {
                cerr << "Error found in input file " << fn << endl;
                abort();
            }
        } else {
            num_object = 0;
        }

        if ( num_hyperedge > 0 )
            *hyperedge = new Entity[ num_hyperedge + 1 ];
        else {
            cerr << "Error found in input file " << fn << endl;
            abort();
        }

        ret.push_back( num_actor );
        ret.push_back( num_object );
        ret.push_back( num_hyperedge );

        while ( file >> id1 >> id2 >> id3 ) {
            if ( false == ah ) {
                (*actor)[ id1 ].insert( id2 );
                (*object)[ id2 ].insert( id1 );
                // insert into the 2nd list, list2
                (*object)[ id2 ].insert( id3, 2 );
                (*hyperedge)[ id3 ].insert( id2 );
            } else {
                (*actor)[ id1 ].insert( id3, 2 );
                (*hyperedge)[ id3 ].insert( id1 );
            }
        }

        file.close();

    } else {
        cerr << "Can not open file " << fn << endl;
        abort();
    }

    return ret;
}

/* 
 * reorganize actor, object, and hyperedge structures for efficiency purpose
 *
 * parameters:
 * Entity **actor, Entity **object, Entity **hyperedge: structure pointer of actor, object, and hyeredge respectively
 * uint32_t size1, uint32_t size2, uint32_t size3: the size of actor, object, and hyperedge respectively
 *
 */
void init_entities ( Entity **actor, uint32_t size1, Entity **object, uint32_t size2, Entity **hyperedge, uint32_t size3 ) {
    uint32_t i = 0;

    for ( i = 1; i <= size1; i++ ) {
        (*actor)[i].re_org();
    }

    for ( i = 1; i <= size2; i++ ) {
        (*object)[i].re_org();
    }

    for ( i = 1; i<= size3; i++ ) {
        (*hyperedge)[i].re_org();
    }
}

/*
 * initialize values for computation
 *
 * parameters:
 * Entity **entity: structure pointer to one entity
 * const uint32_t size: the size of the entity
 * const uint32_t target: target edges for initialization
 *
 */
void init_one_entity_value( Entity **entity, const uint32_t size, const uint32_t target = 1 ) {

    uint32_t i = 0, num1 = 0, num2 = 0;

    for ( i = 0; i < size; i++ ) {
        if ( 0 != (*entity)[ i ].list1.size() )
            num1++;
        if ( 0 != (*entity)[ i ].list2.size() )
            num2++;
    }

    if ( 1 == target )
        for ( i = 0; i < size; i++ ) {
            if ( 0 != (*entity)[ i ].list1.size() )
                (*entity)[ i ].newval = 1.0 / num1;
        }
    else
        for ( i = 0; i < size; i++ ) {
            if ( 0 != (*entity)[ i ].list2.size() )
                (*entity)[ i ].newval = 1.0 / num2;
        }

}

/*
 * compute actor value from object
 *
 * parameters:
 * const double topk: the beta in the paper
 * Entity * const actor, const Entity * const object: structure pointer to actor and object entities
 * const uint32_t size1: the size of the actor entities
 *
 */
void compute_actor_values( const double topk, Entity * const actor, const uint32_t size1, const Entity * const object ) {
    uint32_t i = 0, j = 0;
    double Tk = 0;

    double total = 0;

    // store old values
    // and apply AAPN and sorting if applicable
    for ( i = 1; i <= size1; i++ ) {
        if ( 0 != actor[i].norm1 ) {
            actor[ i ].oldval = actor[ i ].newval;
            actor[ i ].newval = 0;
            for ( j = 0; j < actor[i].norm1; j++ ) {
                actor[i].array1[j].val = object[actor[i].array1[j].idx].newval \
                                              /  object[actor[i].array1[j].idx].norm1;
            }

            // tpk requires the array to be sorted, in decreasing order
            randomized_quicksort( actor[i].array1, 0, actor[i].norm1 - 1 );
        }
    }

    // compute new values
    for ( i = 1; i <= size1; i++ ) {
        // for simplicity reason, consider the first list only
        if ( 0 != actor[i].norm1 ) {

            // set Tk zero for each author value computation
            Tk = 0;

            // compute topk
            uint32_t size = (int)(topk * actor[i].norm1);
            size = (size > 0) ? size:1;

            for ( j = 0; j < size; j++ )
                Tk += actor[i].array1[j].val;

            actor[i].newval = Tk;

            total += actor[i].newval;
        }
    }

    // normalization
    for ( i = 1; i <= size1; i++ )  {
        if ( 0 != actor[i].norm1)
            // normalization
            actor[i].newval = actor[i].newval / total;
    }

}

/*
 * compute object value from actor
 *
 * parameters:
 * Entity * const object, const Entity * const actor: structure pointer to object and actor entities
 * const uint32_t size1: the size of the object entities
 *
 */
void compute_object_value( Entity * const object, const uint32_t size1, const Entity * const actor ) {
    uint32_t i = 0, j = 0;
    double total = 0;

    // store old values
    for ( i = 1; i <= size1; i++ ) {

        if ( 0 != object[i].norm1 ) {
            object[i].oldval = object[i].newval;
            object[i].newval = 0;

            for ( j = 0; j < object[i].norm1; j++ ) {
                object[i].array1[j].val = actor[object[i].array1[j].idx].newval;
            }
        }
    }

    // compute new values
    for ( i = 1; i <= size1; i++ ) {

        if ( 0 != object[i].norm1 ) {
            for ( j = 0; j < object[i].norm1; j++ )
                object[i].newval += object[i].array1[j].val;

            total += object[i].newval;
        }
    }

    // normalization
    for ( i = 1; i <= size1; i++ ) {
        if ( 0 != object[i].norm1 ) {
            object[i].newval = object[i].newval / total;
        }
    }
}

/*
 * compute hyperedge value from object
 *
 * parameters:
 * const double topk: the beta in the paper
 * Entity * const hyperedge, const Entity * const object: structure pointer to hyperedge and object entities
 * const uint32_t size1: the size of the hyperedge entities
 *
 */
void compute_hyperedge_value( const double topk, Entity * const hyperedge, const uint32_t size1, const Entity * const object ) {
    uint32_t i = 0, j = 0;
    double Tk = 0;

    double total = 0;

    // store old values
    // sort if applicable
    for ( i = 1; i <= size1; i++ ) {
        if ( 0 != hyperedge[i].norm1 ) {
            hyperedge[i].oldval = hyperedge[i].newval;
            hyperedge[i].newval = 0;

            for ( j = 0; j < hyperedge[i].norm1; j++ ) {
                // no PTN or something like AAPN applied here
                hyperedge[i].array1[j].val = object[hyperedge[i].array1[j].idx].newval;
            }

           randomized_quicksort( hyperedge[i].array1, 0, hyperedge[i].norm1 - 1 );
           
        }
    }

    // compute new values
    for ( i = 1; i <= size1; i++ ) {
        if ( 0 != hyperedge[i].norm1 ) {
            
            // set Tk zero for each run
            Tk = 0;

            // compute topk & btmk
            uint32_t size = (int)(topk * hyperedge[i].norm1);
            size = (size > 0) ? size:1;

            // summation of topk most valued papers' value
            for ( j = 0; j < size; j++ )
                // PVN is always ON here
                // normalized by size instead of norm1
                Tk += hyperedge[i].array1[j].val / size;


            // hyperedge value
            hyperedge[i].newval = Tk;
           

            total += hyperedge[i].newval;
        }
    }

    // normalization
    for (i = 1; i <= size1; i++) {
        if ( 0 != hyperedge[i].norm1 ) {
            // total normalization
            hyperedge[i].newval = hyperedge[i].newval / total;
        }
    }

}

/*
 * compute object value from hyperedge
 *
 * parameters:
 * Entity * const object, const Entity * const hyperedge: structure pointer to object and hyperedge entities
 * const uint32_t size1: the size of the object entities
 *
 */
void reassign_object_value( Entity * const object, const uint32_t size1, const Entity * const hyperedge) {
    uint32_t i = 0, j = 0;
    double avg = 0;

    for ( i = 1; i <= size1; i++ ) {
        avg = 0;


        if ( 0 != object[i].norm1 ) {
            for ( j = 0; j < object[i].norm2; j++ )
                avg += hyperedge[object[i].array2[j].idx].newval;

            object[i].newval = avg / object[i].norm2;

        }

    }

}

/*
 * convergence or max number of iterations check
 *
 * parameters:
 * const double threshold: threshold value for stopping the program
 * const int iterations: max number of iterations for stopping the program
 * const Entity * const actor: structure pointer to entities
 * const uint32_t size1: the size of the actor entities
 *
 */
const bool satisfied( const double threshold, const int iterations, const Entity * const actor, const uint32_t size1 ) {
    uint32_t i = 0;
    double avg = 0;

    if ( iterations >= MAXITERATIONS )
    	return true;

    int size = 0;
    for ( i = 1; i <= size1; i++ ) {
        if ( 0 != actor[i].norm1 ) {
            avg += fabs( actor[i].oldval - actor[i].newval );
            size ++;
        }
    }

    avg = avg / size;

    return ( avg < threshold );

}

/*
 * print results
 *
 * parameters:
 * const int iterations: real number of iterations
 * const Entity * const actor: structure pointer to entities
 * const uint32_t size1: the size of the actor entities
 *
 */
void print_actor_result( const uint32_t iteration, const Entity * const actor, const uint32_t size1, const bool ah ) {
    // save computation results to files
    uint32_t i = 0, count = 0;
    char fn[255];

    // allocate memory
    Aux *ret = new Aux[ size1 + 1 ];

    // copy value into AUX array
    for ( i = 1; i <= size1; i++ ) {
        if ( false == ah ) {
            if ( 0 != actor[i].norm1 ) {
                ret[count].idx = i;
                ret[count].val = actor[i].newval;
                count++;
            }
        } else {
            if ( 0 != actor[i].norm2 ) {
                ret[count].idx = i;
                ret[count].val = actor[i].newval;
                count++;
            }
        }
    }

    // sorting author value, in-place
    randomized_quicksort( ret, 0, count - 1 );

    // open file for writing
    sprintf( fn, "actor_%d.txt", iteration);

    ofstream file( fn );
    // count = count < lines ? count:lines;
    if ( file.is_open() ) {
        // print_header_info( file, optstr, iteration, alpha, topk );

        for ( i = 0; i < count; i++ )
            file << ret[i].idx << "\t" << ret[i].val << endl;

    } else {
        cerr << "Can not open file " << fn << endl;
    }

    delete [] ret;

    return;

}

/*
 * compute actor value from hyperedge
 *
 * parameters:
 * Entity * const actor, const Entity * const hyperedge: structure pointer to actor and hyperedge entities
 * const uint32_t size1: the size of the actor entities
 *
 */
void reassign_actor_value( Entity * const actor, const uint32_t size1, const Entity * const hyperedge ) {

    uint32_t i = 0, j = 0;
    double avg = 0;

    for ( i = 1; i <= size1; i++ ) {
        avg = 0;


        if ( 0 != actor[i].norm2 ) {
            for ( j = 0; j < actor[i].norm2; j++ )
                avg += hyperedge[actor[i].array2[j].idx].newval;

            actor[i].newval = avg;

        }

    }

    return;
}

/*
 * run test
 *
 * parameters:
 * const double threshold: threshold value for stopping the program
 * const double topk: the beta in the paper
 * Entity * const actor,  Entity * const object, const uint32_t size2, Entity * const hyperedge: pointers to actor, object, hyperedge arrays respectively
   const uint32_t size1, const uint32_t size2, const uint32_t size3: size of actor, object and hyperedge structure respectively
 *
 */
const uint32_t run_test(  const double threshold, const double topk, Entity * const actor, const uint32_t size1, \
        Entity * const object, const uint32_t size2, Entity * const hyperedge, const uint32_t size3, const bool ah ) {

    uint32_t i = 0;

    while ( true ) {
        if ( false == ah ) {
            compute_actor_values( topk, actor, size1, object );
            compute_object_value( object, size2, actor );
            compute_hyperedge_value( topk, hyperedge, size3, object );

            if ( satisfied( threshold, ++i, actor, size1 ) ) {
                break;
            }

            reassign_object_value( object, size2, hyperedge );
        } else {
            compute_hyperedge_value( topk, hyperedge, size3, actor );
            reassign_actor_value( actor, size1, hyperedge );
            if (satisfied(threshold, ++i, actor, size1)) {
                break;
            }
        }
    }

    cout << "Total iterations = " << i << endl;

    return i;
}

/*
 * usage information
 * 
 */
void usage( const char * const arg ){
    cout << "Usage:" << endl;
    cout << arg << " <author_paper_venue_triple_file> [AH]" << endl;
    cout << "where the first line of the input file must indicate the max ids for actors, objects, and hyperedges respectively" << endl;
    cout << "and, followed by positive interger actor_id, object_id, hyperedge_id starting from 1." << endl;
    cout << "furthermore, each entity id must be consecutive." << endl;
    cout << "AH = 0 (default) uses A-O-H network. AH = 1 uses A-H network instead." << endl;
    cout << endl;
    exit( -1 );
}

/*
 * main entry of the program
 *
 */
int main(int argc, char **argv) {
    // default value of alpha, topk, threshold, etc.
    const double topk = 0.50, threshold =  0.000000000001;

    // A-H or A-O-H
    bool ah = false;

    // entities
    Entity *actor = NULL, *object = NULL, *hyperedge = NULL;
    vector<uint32_t> size;
    uint32_t size1, size2, size3;

    uint32_t iteration;

    // init random generator
    srand (19L);

    if ( argc < 2 )
        usage( argv[0] );

    const string fn( argv[1] );

    if ( argc > 2 )
        if ( 1 == atoi(argv[2]) )
            ah = true;

    // set new handler for message output while new failed
    set_new_handler( no_memory );

    cout << "Reading file " << fn << endl;
    size = read_file( fn, &actor, &object, &hyperedge, ah );
    size1 = size[0];
    size2 = size[1];
    size3 = size[2];

    cout << "Initiating entities..." << endl;
    init_entities( &actor, size1, &object, size2, &hyperedge, size3 );

    if ( false == ah ) {
        cout << "Initiating object values..." << endl;
        init_one_entity_value( &object, size2 );
    } else {
        cout << "Initiating actor values..." << endl;
        init_one_entity_value( &actor, size1, 2 );
    }

    cout << "Computing..." << endl;
    iteration = run_test( threshold, topk, actor, size1, object, size2, hyperedge, size3, ah );

    cout << "Saving results..." << endl;
    print_actor_result( iteration, actor, size1, ah );

    if ( actor ) delete []actor;
    if ( object ) delete []object;
    if ( hyperedge ) delete []hyperedge;

    return 0;
}

