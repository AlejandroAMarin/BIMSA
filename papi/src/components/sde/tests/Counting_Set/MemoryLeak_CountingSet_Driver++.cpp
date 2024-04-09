#include <stdio.h>
#include <string.h>
#include "sde_lib.h"
#include "papi.h"
#include "papi_test.h"
#include "cset_lib.hpp"

void setup_PAPI(int *event_set);

int main(int argc, char **argv){
    int cnt, ret, event_set = PAPI_NULL;
    long long counter_values[1];
    CSetLib LibCSetCPP;

    (void)argc;
    (void)argv;

    setup_PAPI(&event_set);

    // --- Start PAPI
    if((ret=PAPI_start(event_set)) != PAPI_OK){
        test_fail( __FILE__, __LINE__, "PAPI_start", ret );
    }

    LibCSetCPP.do_memory_allocations();

    // --- Stop PAPI
    if((ret=PAPI_stop(event_set, counter_values)) != PAPI_OK){
        test_fail( __FILE__, __LINE__, "PAPI_stop", ret );
    }

    if( (argc > 1) && !strcmp(argv[1], "-verbose") ){
        LibCSetCPP.dump_set( (cset_list_object_t *)counter_values[0] );
    }

    cnt = LibCSetCPP.count_set_elements( (cset_list_object_t *)counter_values[0] );
    if( 56 == cnt )
        test_pass(__FILE__);
    else
        test_fail( __FILE__, __LINE__, "CountingSet contains wrong number of elements", ret );

    return 0;
}


void setup_PAPI(int *event_set){
    int ret;

    if((ret=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT){
        test_fail( __FILE__, __LINE__, "PAPI_library_init", ret );
    }

    if((ret=PAPI_create_eventset(event_set)) != PAPI_OK){
        test_fail( __FILE__, __LINE__, "PAPI_create_eventset", ret );
    }

    if((ret=PAPI_add_named_event(*event_set, "sde:::CPP_CSET_LIB::malloc_tracking")) != PAPI_OK){
        test_fail( __FILE__, __LINE__, "PAPI_add_named_event", ret );
    }

    return;
}

