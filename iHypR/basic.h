// -*- C++ -*-
/*
 * File:   basic.h
 * Author: xiaohui
 *
 * Created on March 20, 2011, 2:48 PM
 */

#include <stdint.h>

#include <cstdlib>

#include <iostream>
#include <set>

#ifndef BASIC_H
#define	BASIC_H

typedef struct {
    uint32_t idx;
    double val;
} Aux;


typedef struct _Entity_ {
    // old and new values for comparison
    double oldval;
    double newval;
    // normalization term
    uint32_t norm1;
    uint32_t norm2;
    // list of objects who is pointed to by the entity
    // or, list of objects who point to the entity
    std::set< uint32_t > list1;
    std::set< uint32_t > list2;
    // store entity id and values, two list are corresponding to
    // the lists above
    // since we need to sort the elements in the list, so we want to
    // store them in a separate space
    Aux *array1;
    Aux *array2;

    // construction, called when the object is created
    _Entity_() {
        oldval = newval = 0;
        norm1 = norm2 = 0;
        array1 = array2 = NULL;
    }

    // release memory
    ~_Entity_() {
        if ( array1 )
            delete [] array1;
        if ( array2 )
            delete [] array2;
        array1 = NULL;
        array2 = NULL;
    }

    // target must be 1 or 2 corresponding to list1 or list2
    // duplicate values won't be inserted since the list is a set
    void insert( uint32_t id, int target = 1 ) {
        if ( 1 == target ) {
            list1.insert( id );
        } else if ( 2 == target ) {
            list2.insert( id );
        } else {
            std::cerr << "Invalid link list target " << target << std::endl;
            std::abort();
        }
    }

    // copy elements in list to arrays
    // update normalization term correspondingly
    void re_org( ) {

        uint32_t idx = 0;
        std::set< uint32_t >::iterator it;

        if ( list1.size() > 0 ) {
            array1 = new Aux[ list1.size() ];
            norm1 = list1.size();
        }
        if ( list2.size() > 0 ) {
            array2 = new Aux[ list2.size() ];
            norm2 = list2.size();
        }

        idx = 0;
        for ( it = list1.begin(); it != list1.end(); it++ ) {
            array1[ idx ].idx = *it;
            array1[ idx ].val = 0;
            idx++;
        }

        idx = 0;
        for ( it = list2.begin(); it != list2.end(); it++ ) {
            array2[ idx ].idx = *it;
            array2[ idx ].val = 0;
            idx++;
        }

    }

} Entity;

#endif	/* BASIC_H */

