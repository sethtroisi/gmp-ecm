/*
 ============================================================================
 Name        : LucasChainGen.c
 Author      : P. B. McLaughlin
 Version     : 1.0
 Copyright   : Copyright 2023 by Philip McLaughlin. All rights reserved.
 Description : Minimal-length Lucas chain generator for prime numbers
 ============================================================================
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <time.h>
#include "LucasChainGen.h"
#include "LCG_macros.h"


mem_struct thread_mem[MAX_THREADS];

u_int8_t thread_count;

/* target prime list, used by all threads */
target_prime tgt_prime_list[MAX_CODE_OR_PRIME_COUNT];


pthread_cond_t my_cond = PTHREAD_COND_INITIALIZER;
pthread_mutex_t my_cond_m = PTHREAD_MUTEX_INITIALIZER;

int all_work_done = 0;
int num_waiting_threads = 0;
int *work_request[MAX_THREADS];

/* work assignment array */
work_struct work_assignment[TOTAL_WORK_COUNT];

/* returns current clock count in microseconds */
u_int64_t cputime(void)
{
    return (u_int64_t)clock();
}

void *recursive_work(void *io)
{
	void (*gen_ptr)(void);
	void (*work_ptr)(u_int8_t);
	thread_io_struct *thrd_io;
    int my_work = -1;
	u_int8_t thrd_indx, work_indx;

	thrd_io = (thread_io_struct *)io;
	thrd_indx = thrd_io->thrd_indx;

	switch(thrd_indx)
	{
	case 0:
		gen_ptr = &generate_and_process_candidate_list;
		work_ptr = &copy_work_assignment_to_thread;
		break;
#if MAX_THREADS > 1
	case 1:
		gen_ptr = &generate_and_process_candidate_list_01;
		work_ptr = &copy_work_assignment_to_thread_01;
		break;
#endif
#if MAX_THREADS > 2
	case 2:
		gen_ptr = &generate_and_process_candidate_list_02;
		work_ptr = &copy_work_assignment_to_thread_02;
		break;
#endif
#if MAX_THREADS > 3
	case 3:
		gen_ptr = &generate_and_process_candidate_list_03;
		work_ptr = &copy_work_assignment_to_thread_03;
		break;
#endif
#if MAX_THREADS > 4
	case 4:
		gen_ptr = &generate_and_process_candidate_list_04;
		work_ptr = &copy_work_assignment_to_thread_04;
		break;
#endif
#if MAX_THREADS > 5
	case 5:
		gen_ptr = &generate_and_process_candidate_list_05;
		work_ptr = &copy_work_assignment_to_thread_05;
		break;
#endif
#if MAX_THREADS > 6
	case 6:
		gen_ptr = &generate_and_process_candidate_list_06;
		work_ptr = &copy_work_assignment_to_thread_06;
		break;
#endif
#if MAX_THREADS > 7
	case 7:
		gen_ptr = &generate_and_process_candidate_list_07;
		work_ptr = &copy_work_assignment_to_thread_07;
		break;
#endif
#if MAX_THREADS > 8
	case 8:
		gen_ptr = &generate_and_process_candidate_list_08;
		work_ptr = &copy_work_assignment_to_thread_08;
		break;
#endif
#if MAX_THREADS > 9
	case 9:
		gen_ptr = &generate_and_process_candidate_list_09;
		work_ptr = &copy_work_assignment_to_thread_09;
		break;
#endif
#if MAX_THREADS > 10
	case 10:
		gen_ptr = &generate_and_process_candidate_list_10;
		work_ptr = &copy_work_assignment_to_thread_10;
		break;
#endif
#if MAX_THREADS > 11
	case 11:
		gen_ptr = &generate_and_process_candidate_list_11;
		work_ptr = &copy_work_assignment_to_thread_11;
		break;
#endif
#if MAX_THREADS > 12
	case 12:
		gen_ptr = &generate_and_process_candidate_list_12;
		work_ptr = &copy_work_assignment_to_thread_12;
		break;
#endif
#if MAX_THREADS > 13
	case 13:
		gen_ptr = &generate_and_process_candidate_list_13;
		work_ptr = &copy_work_assignment_to_thread_13;
		break;
#endif
#if MAX_THREADS > 14
	case 14:
		gen_ptr = &generate_and_process_candidate_list_14;
		work_ptr = &copy_work_assignment_to_thread_14;
		break;
#endif
#if MAX_THREADS > 15
	case 15:
		gen_ptr = &generate_and_process_candidate_list_15;
		work_ptr = &copy_work_assignment_to_thread_15;
		break;
#endif
	default:
		printf("ERROR: invalid thread index (%u) in recursive_work\n", thrd_indx);
		return NULL;
	}

	while( !all_work_done)
	{
		pthread_mutex_lock(&my_cond_m);

       	++num_waiting_threads;

       	/* Signal producer to give work */
        work_request[thrd_indx] = &my_work;
        pthread_cond_broadcast(&my_cond);

        /* Wait for work to arrive
           It is wrapped in a while loop because the condition
           might be triggered by another worker thread intended
           to wake up the producer */
        while( (my_work == -1) && !all_work_done)
        	pthread_cond_wait(&my_cond, &my_cond_m);

        --num_waiting_threads;

        /* Work has arrived OR we are done */
        pthread_mutex_unlock(&my_cond_m);

        if(my_work != -1)
        {
            work_indx = (u_int8_t)my_work;
           	(*work_ptr)(work_indx);
           	(*gen_ptr)();
           	my_work = -1;
        }
	}

	pthread_exit(NULL);
}


/* set up variables & variable arrays for multi-routine access */
void init_working_chains(void)
{
	chain_element *working_chain;
	u_int8_t i;

	for( i = 0; i < thread_count; i++ )
	{
		working_chain = thread_mem[i].working_chain;

		working_chain[0].value = 1;

		working_chain[1].value = 2;
		working_chain[1].comp_offset_1 = 0;
		working_chain[1].comp_offset_2 = 0;
		working_chain[1].dif_offset = 0; /* 2 is the only doubled chain element with component offsets both = 0 */
		working_chain[1].chain_dbl_count = 1;

		working_chain[2].value = 3;
		working_chain[2].comp_offset_1 = 0;
		working_chain[2].comp_offset_2 = 1;
		working_chain[2].dif_offset = 1;
		working_chain[2].chain_dbl_count = 1;
	}
}

/* chain start possibilities:
 1 2 3 4 5  c's = 6, 7, 8\, 9*
 1 2 3 4 6  c's = 7*, 8*, 9*, 10X
 1 2 3 4 7  c's = 8*, 10*, 11*
 1 2 3 5 6  c's = 7\, 8\, 9*, 10*, 11*
 1 2 3 5 7  c's = 8\, 9*, 10*, 12*
 1 2 3 5 8  c's = 10*, 11*, 13*
 \ ==> has appeared in a minimal length chain
 * ==> has appeared in a minimal length chain with maximum # of doubled elements
 X ==> only composite continuations
*/
void set_work_assignments(void)
{
	chain_element *working_chain;
	int i;

	i = 0;

	/* CHAIN_START_4_6 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 6;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 0;
	working_chain[4].chain_dbl_count = 3;

	work_assignment[i++].current_partial_length = 4;

	/* CHAIN_START_4_7_8 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 3;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 8;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 0;
	working_chain[5].chain_dbl_count = 3;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_4_7_11_14 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 3;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 2;

	working_chain[6].value = 14;
	working_chain[6].comp_offset_1 = 1;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 0;
	working_chain[6].chain_dbl_count = 3;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_4_7_11_18 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 3;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 2;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_4_7_11_15 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 3;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 2;

	working_chain[6].value = 15;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_4_7_10 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 3;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 10;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_4_5 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 4;
	working_chain[3].comp_offset_1 = 1;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 0;
	working_chain[3].chain_dbl_count = 2;

	working_chain[4].value = 5;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 3;
	working_chain[4].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 4;

	/* CHAIN_START_5_6_10 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 6;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 0;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 10;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 0;
	working_chain[5].chain_dbl_count = 3;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_6_11 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 6;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 0;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 4;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_6_9 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 6;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 0;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 9;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_6_8 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 6;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 0;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 8;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_6_7 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 6;
	working_chain[4].comp_offset_1 = 1;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 0;
	working_chain[4].chain_dbl_count = 2;

	working_chain[5].value = 7;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 3;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_8_10 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 10;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 0;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_8_13_16 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 16;
	working_chain[6].comp_offset_1 = 1;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 0;
	working_chain[6].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_5_8_13_21_26 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 26;
	working_chain[7].comp_offset_1 = 1;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 0;
	working_chain[7].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_13_21_34_42 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 42;
	working_chain[8].comp_offset_1 = 1;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 0;
	working_chain[8].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_21_34_55_68 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 68;
	working_chain[9].comp_offset_1 = 1;
	working_chain[9].comp_offset_2 = 1;
	working_chain[9].dif_offset = 0;
	working_chain[9].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 9;

	/* CHAIN_START_5_8_13_21_34_55_89_110 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 89;
	working_chain[9].comp_offset_1 = 0;
	working_chain[9].comp_offset_2 = 1;
	working_chain[9].dif_offset = 2;
	working_chain[9].chain_dbl_count = 1;

	working_chain[10].value = 110;
	working_chain[10].comp_offset_1 = 1;
	working_chain[10].comp_offset_2 = 1;
	working_chain[10].dif_offset = 0;
	working_chain[10].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 10;

	/* CHAIN_START_5_8_13_21_34_55_89_144_178 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 89;
	working_chain[9].comp_offset_1 = 0;
	working_chain[9].comp_offset_2 = 1;
	working_chain[9].dif_offset = 2;
	working_chain[9].chain_dbl_count = 1;

	working_chain[10].value = 144;
	working_chain[10].comp_offset_1 = 0;
	working_chain[10].comp_offset_2 = 1;
	working_chain[10].dif_offset = 2;
	working_chain[10].chain_dbl_count = 1;

	working_chain[11].value = 178;
	working_chain[11].comp_offset_1 = 1;
	working_chain[11].comp_offset_2 = 1;
	working_chain[11].dif_offset = 0;
	working_chain[11].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 11;

	/* CHAIN_START_5_8_13_21_34_55_89_144_233 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 89;
	working_chain[9].comp_offset_1 = 0;
	working_chain[9].comp_offset_2 = 1;
	working_chain[9].dif_offset = 2;
	working_chain[9].chain_dbl_count = 1;

	working_chain[10].value = 144;
	working_chain[10].comp_offset_1 = 0;
	working_chain[10].comp_offset_2 = 1;
	working_chain[10].dif_offset = 2;
	working_chain[10].chain_dbl_count = 1;

	working_chain[11].value = 233;
	working_chain[11].comp_offset_1 = 0;
	working_chain[11].comp_offset_2 = 1;
	working_chain[11].dif_offset = 2;
	working_chain[11].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 11;

	/* CHAIN_START_5_8_13_21_34_55_89_144_199 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 89;
	working_chain[9].comp_offset_1 = 0;
	working_chain[9].comp_offset_2 = 1;
	working_chain[9].dif_offset = 2;
	working_chain[9].chain_dbl_count = 1;

	working_chain[10].value = 144;
	working_chain[10].comp_offset_1 = 0;
	working_chain[10].comp_offset_2 = 1;
	working_chain[10].dif_offset = 2;
	working_chain[10].chain_dbl_count = 1;

	working_chain[11].value = 199;
	working_chain[11].comp_offset_1 = 0;
	working_chain[11].comp_offset_2 = 2;
	working_chain[11].dif_offset = 1;
	working_chain[11].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 11;

	/* CHAIN_START_5_8_13_21_34_55_89_123 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 89;
	working_chain[9].comp_offset_1 = 0;
	working_chain[9].comp_offset_2 = 1;
	working_chain[9].dif_offset = 2;
	working_chain[9].chain_dbl_count = 1;

	working_chain[10].value = 123;
	working_chain[10].comp_offset_1 = 0;
	working_chain[10].comp_offset_2 = 2;
	working_chain[10].dif_offset = 1;
	working_chain[10].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 10;

	/* CHAIN_START_5_8_13_21_34_55_76 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 55;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	working_chain[9].value = 76;
	working_chain[9].comp_offset_1 = 0;
	working_chain[9].comp_offset_2 = 2;
	working_chain[9].dif_offset = 1;
	working_chain[9].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 9;

	/* CHAIN_START_5_8_13_21_34_47 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 34;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 47;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 2;
	working_chain[8].dif_offset = 1;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_21_29_42 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 29;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 42;
	working_chain[8].comp_offset_1 = 1;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 0;
	working_chain[8].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_21_29_50 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 29;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 50;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 3;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_21_29_37 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 29;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 37;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 3;
	working_chain[8].dif_offset = 1;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_21_29_34 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 21;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 29;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 34;
	working_chain[8].comp_offset_1 = 1;
	working_chain[8].comp_offset_2 = 2;
	working_chain[8].dif_offset = 3;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_18_26 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 26;
	working_chain[7].comp_offset_1 = 1;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 0;
	working_chain[7].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_13_18_31_36 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 31;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 3;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 36;
	working_chain[8].comp_offset_1 = 1;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 0;
	working_chain[8].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_18_31_49 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 31;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 3;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 49;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_18_31_44 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 31;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 3;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 44;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 2;
	working_chain[8].dif_offset = 1;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_8_13_18_23 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 23;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 3;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_13_18_21 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 13;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 2;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 18;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 21;
	working_chain[7].comp_offset_1 = 1;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 3;
	working_chain[7].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_11_16 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 16;
	working_chain[6].comp_offset_1 = 1;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 0;
	working_chain[6].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_5_8_11_19_22 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 3;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 22;
	working_chain[7].comp_offset_1 = 1;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 0;
	working_chain[7].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_11_19_30 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 3;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 30;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_11_19_27 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 3;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 27;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_8_11_14 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 14;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 3;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_5_8_11_13 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 8;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 1;
	working_chain[4].dif_offset = 2;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 11;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 13;
	working_chain[6].comp_offset_1 = 1;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 3;
	working_chain[6].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_5_7_10 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 10;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 0;
	working_chain[5].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_7_12_14 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 14;
	working_chain[6].comp_offset_1 = 1;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 0;
	working_chain[6].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_5_7_12_19_24 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 24;
	working_chain[7].comp_offset_1 = 1;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 0;
	working_chain[7].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_7_12_19_31_38 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 31;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 38;
	working_chain[8].comp_offset_1 = 1;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 0;
	working_chain[8].chain_dbl_count = 2;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_7_12_19_31_50 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 31;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 50;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 1;
	working_chain[8].dif_offset = 2;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_7_12_19_31_43 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 31;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 1;
	working_chain[7].dif_offset = 2;
	working_chain[7].chain_dbl_count = 1;

	working_chain[8].value = 43;
	working_chain[8].comp_offset_1 = 0;
	working_chain[8].comp_offset_2 = 2;
	working_chain[8].dif_offset = 1;
	working_chain[8].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 8;

	/* CHAIN_START_5_7_12_19_26 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 19;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 1;
	working_chain[6].dif_offset = 2;
	working_chain[6].chain_dbl_count = 1;

	working_chain[7].value = 26;
	working_chain[7].comp_offset_1 = 0;
	working_chain[7].comp_offset_2 = 2;
	working_chain[7].dif_offset = 1;
	working_chain[7].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 7;

	/* CHAIN_START_5_7_12_17 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 12;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 1;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	working_chain[6].value = 17;
	working_chain[6].comp_offset_1 = 0;
	working_chain[6].comp_offset_2 = 2;
	working_chain[6].dif_offset = 1;
	working_chain[6].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 6;

	/* CHAIN_START_5_7_9 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 9;
	working_chain[5].comp_offset_1 = 0;
	working_chain[5].comp_offset_2 = 3;
	working_chain[5].dif_offset = 1;
	working_chain[5].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 5;

	/* CHAIN_START_5_7_8 */
	working_chain = work_assignment[i].working_chain;

	working_chain[3].value = 5;
	working_chain[3].comp_offset_1 = 0;
	working_chain[3].comp_offset_2 = 1;
	working_chain[3].dif_offset = 2;
	working_chain[3].chain_dbl_count = 1;

	working_chain[4].value = 7;
	working_chain[4].comp_offset_1 = 0;
	working_chain[4].comp_offset_2 = 2;
	working_chain[4].dif_offset = 1;
	working_chain[4].chain_dbl_count = 1;

	working_chain[5].value = 8;
	working_chain[5].comp_offset_1 = 1;
	working_chain[5].comp_offset_2 = 2;
	working_chain[5].dif_offset = 3;
	working_chain[5].chain_dbl_count = 1;

	work_assignment[i++].current_partial_length = 5;

	printf("Total number of work asssignments = %d\n", i );
}

/* Fibonacci numbers */
void init_Fib_sequence(void)
{
	u_int64_t *Fib;
	int32_t i, k;

	for( k = 0; k < thread_count; k++ )
	{
		Fib = thread_mem[k].Fib;

		Fib[0] = 0;
		Fib[1] = 1;
		for(i=2;i<FIB_LIMIT;i++)
		{
			Fib[i] = Fib[i-1] + Fib[i-2];
		}
	}
}

/* Lucas numbers */
void init_Luc_sequence(void)
{
	u_int64_t *Luc;
	int32_t i, k;

	for( k = 0; k < thread_count; k++ )
	{
		Luc = thread_mem[k].Luc;

		Luc[0] = 2;
		Luc[1] = 1;
		for(i=2;i<FIB_LIMIT;i++)
		{
			Luc[i] = Luc[i-1] + Luc[i-2];
		}
	}
}

u_int8_t *get_dif_table_ptr(void)
{
	static u_int8_t dif_table[5760];
	return dif_table;
}

u_int8_t *get_sieve_space_ptr(void)
{
	static u_int8_t sieve_space[SIEVE_SPACE_SIZE];
	return sieve_space;
}

sieve_params *get_sieve_primes_ptr(void)
{
	static sieve_params sieve_primes[SIEVE_PRIME_COUNT];
	return sieve_primes;
}


/* sieving & prime generation routines */
u_int32_t sieve_init(void)
{
	u_int8_t *dif_table, *newsieve;
	sieve_params *sieve_primes;
	u_int32_t i, j;
	u_int32_t s_index, indx_p, last_index, dif_index, dif_sum, max_dif, max_dif_sum = 0;
	u_int32_t small_primes[5] = {3, 5, 7, 11, 13};
	u_int32_t p = 0;
	u_int64_t tmp;

	dif_table = get_dif_table_ptr();
	newsieve = get_sieve_space_ptr();
	sieve_primes = get_sieve_primes_ptr();

	/* set up the difference table, sieve primes, and associated table indices */

	/* construct a table of differences for integers not divisible by
	any of the primes 2, 3, 5, 7, 11, or 13. Note that sieve arrays
	are half size since we will only be dealing with odd integers.
	Each byte in the sieve array will represent a prime/composite
	flag for the integer 2*s_index+1. s_index=0 corresponds to 1, s_index=1 to 3, etc. */

	/* for reference: pi(10^3) = 168, pi(10^4) = 1229, pi(10^5) = 9592
	pi(10^6) = 78498, pi(10^7) = 664579, pi(10^8) = 5761455,
	pi(10^9) = 50847534, pi(10^10) = 455052511, pi(2*10^10) = 882206716
	 */

	/* initialize array */
	for(j=0;j<NEWSIEVE_SIZE;j++) newsieve[j] = 1;

	/* cross out all (odd) multiples of the small primes */
	for(i=0;i<5;i++)
	{
		s_index = (small_primes[i] - 1)/2;
		do
		{
			newsieve[s_index] = 0;
			s_index += small_primes[i];
		}
		while(s_index < NEWSIEVE_SIZE);
	}

	/* populate the difference table */
	last_index = 0;
	s_index = 1;
	j = 0;
	do
	{
		while(newsieve[s_index] == 0) s_index++;
		dif_table[j] = (u_int8_t)(s_index - last_index);
		last_index = s_index;
		j++;
		s_index++;
	}
	while( j < 5760 );

	printf("dif[0] = %d, dif[5758] = %d, dif[5759] = %d\n", dif_table[0], dif_table[5758], dif_table[5759]);

	dif_sum = 0;
	max_dif = 0;
	for(j=0;j<5760;j++)
	{
		dif_sum += dif_table[j];
		if( dif_table[j] > max_dif )
		{
			max_dif = dif_table[j];
			max_dif_sum = dif_sum;
		}
	}
	printf("sum of differences = %d, max dif = %d at sum = %d\n", dif_sum, max_dif, max_dif_sum);

	/* starting with 17, find all primes < SBL and complete sieving to SBL^2 = SPL */
	s_index = 0;
	dif_index = 0;
	while( (s_index = s_index + dif_table[dif_index]) < SBL/2)
	{
		dif_index++;
		if(dif_index == 5760) dif_index = 0;
		if( newsieve[s_index] )
		{
			p = 2*s_index + 1;
			indx_p = (p*p - 1)/2;
			do
			{
				newsieve[indx_p] = 0;
				indx_p += p;
			}
			while(indx_p < NEWSIEVE_SIZE);
		}
	}
	printf("Final bootstrap prime = %u\n", p);

	/* populate the sieve prime struct array */
	i = 0;
	s_index = 8;
	dif_index = 1;
	do
	{
		if( newsieve[s_index] )
		{
			p = 2*s_index + 1;
			sieve_primes[i].prime = p;
			tmp = (u_int64_t)p;
			sieve_primes[i].sieve_space_start_index = (tmp*tmp - 1)/2;
			sieve_primes[i].dif_table_start_index = dif_index;
			i++;
		}
		s_index += dif_table[dif_index];
		dif_index++;
		if(dif_index == 5760) dif_index = 0;
	}
	while( s_index < NEWSIEVE_SIZE );

	printf(" check 10^3: p(168) = %u, p(169) = %u\n",sieve_primes[161].prime, sieve_primes[162].prime);
	printf(" check 10^4: p(1229) = %u, p(1230) = %u\n",sieve_primes[1222].prime, sieve_primes[1223].prime);
	printf(" check 10^5: p(9592) = %u, p(9593) = %u\n",sieve_primes[9585].prime, sieve_primes[9586].prime);
	printf("Total sieve primes = %d, pi(%d) = %d, last p = %u, last sieve_space_start_index = %lu\n",
			i, 2*NEWSIEVE_SIZE, i+6, sieve_primes[i-1].prime, sieve_primes[i-1].sieve_space_start_index);

	return i;
}

void standard_sieve( u_int32_t sieve_prime_count )
{
	static u_int8_t *dif_table, *sieve_space, init = 0;
	static sieve_params *sieve_primes;

	u_int32_t i, j, k, dif_index;
	u_int64_t indx;
	u_int32_t p, p_multiples[11];

	if(init == 0)
	{
		dif_table = get_dif_table_ptr();
		sieve_space = get_sieve_space_ptr();;
		sieve_primes = get_sieve_primes_ptr();
		init = 1;
	}

	/* initialize the sieve space */
	/* note that it is not necessary to clear any multiples of 3,5,7,11, or 13 */
	for(i=0;i<SIEVE_SPACE_SIZE;i++) sieve_space[i] = 1;

	/* the following code should handle any SIEVE_SPACE_SIZE */
	for(i=0;i<sieve_prime_count;i++)
	{
		indx = sieve_primes[i].sieve_space_start_index;
		if(indx < SIEVE_SPACE_SIZE)
		{
			/* sieve */
			p = sieve_primes[i].prime;
			j = p;
			for(k=0;k<11;k++)
			{
				p_multiples[k] = j;
				j += p;
			}
			dif_index = sieve_primes[i].dif_table_start_index;
			do
			{
				sieve_space[indx] = 0;
				indx += p_multiples[dif_table[dif_index] - 1];
				dif_index++;
				if( dif_index == 5760 ) dif_index = 0;
			}
			while( indx < SIEVE_SPACE_SIZE );
			sieve_primes[i].dif_table_start_index = dif_index;
		}
		/* reset start index for next interval */
		sieve_primes[i].sieve_space_start_index = indx - SIEVE_SPACE_SIZE;
	}
}

u_int32_t prime_count( u_int32_t *sieve_space_start_index, u_int32_t *dif_table_start_index )
{
	static u_int8_t *dif_table, *sieve_space, init = 0;

	u_int32_t indx, dif_index, p_count;

	if(init == 0)
	{
		dif_table = get_dif_table_ptr();
		sieve_space = get_sieve_space_ptr();
		init = 1;
	}

	indx = *sieve_space_start_index;
	dif_index = *dif_table_start_index;
	p_count = 0;
	do
	{
		if( sieve_space[indx] ) p_count++;
		indx += dif_table[dif_index];
		dif_index++;
		if( dif_index == 5760 ) dif_index = 0;
	}
	while( indx < SIEVE_SPACE_SIZE);

	/* return indices for next interval */
	*sieve_space_start_index = indx - SIEVE_SPACE_SIZE;
	*dif_table_start_index = dif_index;

	return p_count;
}

/* routines to re-generate Lucas chain from code. Used for test only */
#if 0
/* generate Lucas chains for primes >= 5. Returns chain length */
u_int8_t generate_Lchain( u_int64_t prime, u_int64_t chain_code, chain_element *Lchain, u_int8_t *dbl_count, u_int8_t *max_index, u_int32_t *start_frag_count )
{
	static u_int8_t init = 0;
	u_int8_t code_fragment, chain_length, i, k;
	u_int64_t dif;
	u_int64_t chain_code_save;

	if( init == 0 )
	{
		/* the first 3 chain elements are always value = 1, 2, and 3, respectively */
		Lchain[0].value = 1;

		Lchain[1].value = 2;
		Lchain[1].comp_offset_1 = 0;
		Lchain[1].comp_offset_2 = 0;
		Lchain[1].dif_offset = 0;

		Lchain[2].value = 3;
		Lchain[2].comp_offset_1 = 0;
		Lchain[2].comp_offset_2 = 1;
		Lchain[2].dif_offset = 1;

		init = 1;
	}
	*dbl_count = 1;
	chain_code_save = chain_code;

	/* the code file starts at p = 11, so handle 5 and 7 separately */
	if( prime < 11 )
	{
		if( prime == 5 )
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			chain_length = 3;
			return chain_length;
		}
		else if( prime == 7 )
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;
			(*dbl_count)++;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 3;

			chain_length = 4;
			return chain_length;
		}
		else
		{
			printf("ERROR: generate_Lchain entered with prime < 11 but != 5 or 7\n");
			return 0;
		}
	}

	/* first 3 bits of code give the next two or three chain components */
	code_fragment = (u_int8_t)(chain_code & 0x7);
	chain_code >>= 3;
	switch( code_fragment )
	{
	case CHAIN_START_5_8_13:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 13;
		Lchain[5].comp_offset_1 = 0;
		Lchain[5].comp_offset_2 = 1;
		Lchain[5].dif_offset = 2;

		chain_length = 5;
		start_frag_count[CHAIN_START_5_8_13]++;
		break;

	case CHAIN_START_5_8_11:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 11;
		Lchain[5].comp_offset_1 = 0;
		Lchain[5].comp_offset_2 = 2;
		Lchain[5].dif_offset = 1;

		chain_length = 5;
		start_frag_count[CHAIN_START_5_8_11]++;
		break;

	case CHAIN_START_5_8_10:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 8;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 2;

		Lchain[5].value = 10;
		Lchain[5].comp_offset_1 = 1;
		Lchain[5].comp_offset_2 = 1;
		Lchain[5].dif_offset = 0;

		(*dbl_count)++;
		chain_length = 5;
		start_frag_count[CHAIN_START_5_8_10]++;
		break;

	case CHAIN_START_5_7:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 7;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 2;
		Lchain[4].dif_offset = 1;

		chain_length = 4;
		start_frag_count[CHAIN_START_5_7]++;
		break;

	case CHAIN_START_5_6:
		Lchain[3].value = 5;
		Lchain[3].comp_offset_1 = 0;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 2;

		Lchain[4].value = 6;
		Lchain[4].comp_offset_1 = 1;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 0;
		(*dbl_count)++;

		chain_length = 4;
		start_frag_count[CHAIN_START_5_6]++;
		break;

	case CHAIN_START_4_7:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;
		(*dbl_count)++;

		Lchain[4].value = 7;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 3;

		chain_length = 4;
		start_frag_count[CHAIN_START_4_7]++;
		break;

	case CHAIN_START_4_6:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;

		Lchain[4].value = 6;
		Lchain[4].comp_offset_1 = 1;
		Lchain[4].comp_offset_2 = 1;
		Lchain[4].dif_offset = 0;

		*dbl_count += 2;
		chain_length = 4;
		start_frag_count[CHAIN_START_4_6]++;
		break;

	case CHAIN_START_4_5:
		Lchain[3].value = 4;
		Lchain[3].comp_offset_1 = 1;
		Lchain[3].comp_offset_2 = 1;
		Lchain[3].dif_offset = 0;
		(*dbl_count)++;

		Lchain[4].value = 5;
		Lchain[4].comp_offset_1 = 0;
		Lchain[4].comp_offset_2 = 3;
		Lchain[4].dif_offset = 1;

		chain_length = 4;
		start_frag_count[CHAIN_START_4_5]++;
		break;

	default:
		printf("ERROR: bad chain code start value in generate_Lchain = %u\n", code_fragment);
		return 0;
	}

	/* rebuild chain from code fragments */
	while( chain_code != 0 )
	{
		code_fragment = (u_int8_t)( chain_code & 0xF );
		chain_code >>= 4;
		switch( code_fragment )
		{
			case 0: /* step type 1 or 4 */
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				if( i == 3 )
				{
					 /* a rare double, but does occur */
					Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-2 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 2;
					Lchain[ chain_length+1 ].dif_offset = 0;
					(*dbl_count)++;
					chain_length++;
				}
				else
				{
					i = (u_int8_t)(12*(i + 1));
					max_continuation( Lchain, &chain_length, i );
				}
				break;
			}
			case 1:
			{
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				dif = Lchain[ chain_length ].value - Lchain[ chain_length-2 ].value;
				k = 1;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
				if( dif != Lchain[ chain_length-k ].value )
				{
					printf("ERROR: invalid code fragment 1, dif not in chain!\n");
					return 0;
				}
				Lchain[ chain_length+1 ].dif_offset = k;
				chain_length++;
				break;
			}
			case 2:
			{
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				(*dbl_count)++;
				chain_length++;
				break;
			}
			case 3:
			{
				i = 1;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 4:
			{
				i = 1;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				(*dbl_count)++;
				chain_length++;
				break;
			}
			case 5:
			{
				i = 2;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 6:
			{
				i = 2;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				(*dbl_count)++;
				chain_length++;
				break;
			}
			case 7:
			{
				i = 3;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 8:
			{
				i = 3;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				(*dbl_count)++;
				chain_length++;
				break;
			}
			case 9:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 4;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 10:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 4;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				(*dbl_count)++;
				chain_length++;
				break;
			}
			case 11:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 8;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-2 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 2;
				Lchain[ chain_length+1 ].dif_offset = 1;
				chain_length++;
				break;
			}
			case 12:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;
				i += 8;
				max_continuation( Lchain, &chain_length, i );
				Lchain[ chain_length+1 ].value = 2*Lchain[ chain_length-1 ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 1;
				Lchain[ chain_length+1 ].dif_offset = 0;
				(*dbl_count)++;
				chain_length++;
				break;
			}
			case 13:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				switch( i )
				{
				case 0:
				{
					Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-7 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 0;
					Lchain[ chain_length+1 ].comp_offset_2 = 7;
					dif = Lchain[ chain_length ].value - Lchain[ chain_length-7 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					if( dif != Lchain[ chain_length-k ].value )
					{
						printf("ERROR: invalid code fragment 0x0D, dif not in chain!\n");
						return 0;
					}
					Lchain[ chain_length+1 ].dif_offset = k;

					chain_length++;
					printf("Info: code fragment 0x0D in chain for p = %lu, code = %016lX\n", prime, chain_code_save);
					break;
				}
				case 1:
				{
					Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-8 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 0;
					Lchain[ chain_length+1 ].comp_offset_2 = 8;
					dif = Lchain[ chain_length ].value - Lchain[ chain_length-8 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					if( dif != Lchain[ chain_length-k ].value )
					{
						printf("ERROR: invalid code fragment 0x1D, dif not in chain!\n");
						return 0;
					}
					Lchain[ chain_length+1 ].dif_offset = k;

					chain_length++;
					printf("Info: code fragment 0x1D in chain for p = %lu, code = %016lX\n", prime, chain_code_save);
					break;
				}
				case 2: /* can occur but don't have an example (yet) */
				{
					Lchain[ chain_length+1 ].value = Lchain[ chain_length - 2 ].value + Lchain[ chain_length - 3 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 3;
					dif = Lchain[ chain_length - 2 ].value - Lchain[ chain_length - 3 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					if( dif != Lchain[ chain_length-k ].value )
					{
						printf("ERROR: invalid code fragment 0x2D, dif not in chain!\n");
						return 0;
					}
					Lchain[ chain_length+1 ].dif_offset = k;

					chain_length++;
					printf("Info: code fragment 0x2D in chain for p = %lu, code = %016lX\n", prime, chain_code_save);
					break;
				}
				case 3: /* can occur but don't have an example (yet) */
				{
					Lchain[ chain_length+1 ].value = Lchain[ chain_length - 2 ].value + Lchain[ chain_length - 4 ].value;
					Lchain[ chain_length+1 ].comp_offset_1 = 2;
					Lchain[ chain_length+1 ].comp_offset_2 = 4;
					dif = Lchain[ chain_length - 2 ].value - Lchain[ chain_length - 4 ].value;
					k = 1;
					while( dif < Lchain[ chain_length-k ].value )
						k++;
					if( dif != Lchain[ chain_length-k ].value )
					{
						printf("ERROR: invalid code fragment 0x2D, dif not in chain!\n");
						return 0;
					}
					Lchain[ chain_length+1 ].dif_offset = k;

					chain_length++;
					printf("Info: code fragment 0x3D in chain for p = %lu, code = %016lX\n", prime, chain_code_save);
					break;
				}
				default:
					{
						printf("ERROR: unimplemented code fragment 0xiD, i = %u\n", i);
						return 0;
					}
				}
				break;
			}
			case 14:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				Lchain[ chain_length+1 ].value = Lchain[ chain_length ].value + Lchain[ chain_length-3-i ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 0;
				Lchain[ chain_length+1 ].comp_offset_2 = 3+i;
				dif = Lchain[ chain_length ].value - Lchain[ chain_length-3-i ].value;
				k = 1;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
				if( dif != Lchain[ chain_length-k ].value )
				{
					printf("ERROR: invalid code fragment 14, dif not in chain!\n");
					return 0;
				}
				Lchain[ chain_length+1 ].dif_offset = k;

				chain_length++;
				break;
			}
			case 15:
			{
				i = (u_int8_t)( chain_code & 0x3 );
				chain_code >>= 2;

				Lchain[ chain_length+1 ].value = Lchain[ chain_length-1 ].value + Lchain[ chain_length-2-i ].value;
				Lchain[ chain_length+1 ].comp_offset_1 = 1;
				Lchain[ chain_length+1 ].comp_offset_2 = 2+i;
				dif = Lchain[ chain_length-1 ].value - Lchain[ chain_length-2-i ].value;
				k = 2;
				while( dif < Lchain[ chain_length-k ].value )
					k++;
				if( dif != Lchain[ chain_length-k ].value )
				{
					printf("ERROR: invalid code fragment 15, dif not in chain!\n");
					return 0;
				}
				Lchain[ chain_length+1 ].dif_offset = k;

				chain_length++;
				break;
			}
			default: /* placeholder - should never get here */
			{
				printf("ERROR: reached default #1\n");
				break;
			}
		}
	}

	/* finish the chain with maximum continuations until prime is reached */
	i = 1;
	while( Lchain[chain_length].value < prime )
	{
		max_continuation( Lchain, &chain_length, i );
	}

	if( Lchain[ chain_length ].value != prime )
	{
		printf("ERROR: chain code mismatch with prime or invalid chain code!\n");
		return 0;
	}

	/* see if any dif_offset indices exceed the current maximum */
	for( i = 6; i <= chain_length; i++ )
		if( Lchain[ i ].dif_offset > *max_index )
		{
			*max_index = Lchain[ i ].dif_offset;
			printf("# of older x,z states to save = %u for prime = %lu, code = %016lX\n", *max_index, prime, chain_code_save);
			for(k = 0; k <= chain_length; k++)
				printf(" %lu", Lchain[ k ].value);
			printf("\n chain length = %u, doubles count = %u\n", chain_length, *dbl_count);
		}

	if( prime == 710559673 )
	{
		printf("PRAC performance worst for prime = %lu, optimal code = %016lX\n", prime, chain_code_save);
		for(k = 0; k <= chain_length; k++)
			printf(" %lu", Lchain[ k ].value);
		printf("\n chain length = %u, doubles count = %u\n", chain_length, *dbl_count);
	}

	return chain_length;
}

/* extend the chain with maximum elements for i steps */
void max_continuation( chain_element *Lchain, u_int8_t *chain_length, u_int8_t i )
{
	u_int8_t k;
	u_int64_t dif;

	Lchain[ *chain_length+1 ].value = Lchain[ *chain_length ].value + Lchain[ *chain_length-1 ].value;
	Lchain[ *chain_length+1 ].comp_offset_1 = 0;
	Lchain[ *chain_length+1 ].comp_offset_2 = 1;
	dif = Lchain[ *chain_length ].value - Lchain[ *chain_length-1 ].value;
	k = 2;
	while( dif < Lchain[ *chain_length-k ].value )
		k++;
	if( dif != Lchain[ *chain_length-k ].value )
	{
		printf("ERROR: invalid code fragment, max continuation dif not in chain!\n");
		return;
	}
	Lchain[ *chain_length+1 ].dif_offset = k;
	(*chain_length)++;
	if( i > 1 )
	{
		for( k = 1; k < i; k++)
		{
			Lchain[ *chain_length+1 ].value = Lchain[ *chain_length ].value + Lchain[ *chain_length-1 ].value;
			Lchain[ *chain_length+1 ].comp_offset_1 = 0;
			Lchain[ *chain_length+1 ].comp_offset_2 = 1;
			Lchain[ *chain_length+1 ].dif_offset = 2;
			(*chain_length)++;
		}
	}
}
#endif

/* the workhorse recursive routine to find Lucas chains for primes */
/*	1. generate candidates for the current working chain.
	2. Filter out candidates whose maximal continuation is less than tgt_prime_list[0].prime.
	3. Filter out candidates whose continuation is always composite.
	4. Once the length of the current working chain has reached *w_chain_length - 1, compare
	   valid candidates to primes in the target prime list and record any matches.
	4. TBD If certain constants related to extending the partial chain containing a particular candidate
	   fall below tgt_prime_list[0].prime, then the final values for all chains with length == test_length,
	   final value >= tgt_prime_list[0].prime, and containing this partial chain can be computed immediately.
	   All such final values will be compared to the list of target primes to find matches, and the
	   optimal full chains with matching prime values saved. */

void generate_and_process_candidate_list(void)
GEN_AND_PROCESS_C_LIST_1(0)
	c_count[r_level] = gen_candidate_list();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain(void)
COPY_C_TO_W_CHAIN(0)

u_int8_t extract_chain_values(void)
EXTRACT_CHAIN_VALUES(0)

/* to be as short as possible, chain codes are "incomplete," i.e.
 * we assume that the final prime value will be available when the
 * Lucas chain is reconstructed from the code. and any "missing" steps
 * will consist of a(i+1) = a(i) + a(i-1) for each i from the length
 * reached with the code until the final prime value is obtained.
 * Note that we may enter this routine with a partial chain as long as
 * the continuation of that chain consists of maximal steps only.
 *
 *  step types for next chain element
 *  indexed here from largest a(0) to smallest
 *  Note that no other step type has been observed (yet)
 *  in a shortest possible chain. Let

 Type 1 = a(0) + a(1)
 Type 2 = a(0) + a(2)
 Type 3 = a(1) + a(1) (dbl)
 Type 4 = a(2) + a(2) (dbl)
 Type 5 = a(0) + a(c), c = 3, 4, 5, 6 or 7 (2 bit count + 3; count = next 2 bits after E code)
 Type 6 = a(1) + a(c), c = 2, 3, 4, or 5 (2 bit count + 2; count = next 2 bits after F code)
 Type 7 = a(2) + a(3)

 Chain steps are encoded from largest to smallest elements, but
 chains are reconstructed smallest to largest from the code.
 4-bit and 6-bit codes (xx is binary) applied right to left:

For chain elements ai , i >= 6, we use the following:

0xk0.  k = 0, 1, or 2 ==> string of 12*(k + 1) type 1 elements.
           k = 3 ==> type 4.
0x1.  Type 2.
0x2.  Type 3.
0x3.  One type 1 followed by type 2.
0x4.  One type 1 followed by type 3.
0x5.  Two type 1 followed by type 2.
0x6.  Two type 1 followed by type 3.
0x7.  Three type 1 followed by type 2.
0x8.  Three type 1 followed by type 3.
0xk9.  (k + 4)  type 1 followed by type 2.
0xkA.  (k + 4)  type 1 followed by type 3.
0xkB.  (k + 8)  type 1 followed by type 2.
0xkC.  (k + 8)  type 1 followed by type 3.
0x0D.  Type 5 with k = i - 7.
0x1D   Type 7
0xkD, k > 1.  Reserved.
0xjE.  Type 5 with k = i – j - 3.
0xjF.  Type 6 with k = i – j - 2.
*/
u_int64_t encode_Lchain(void)
ENCODE_LCHAIN_TEMPLATE(0)

/* This routine returns 0 if arg % 3 = 0, 1 otherwise */
/* NOTE: arg must be <= 3*2^32 + 2 */
u_int8_t not_divisible_by_3( u_int64_t arg )
NOT_DIVISIBLE_3

/* This routine returns 0 if arg % 5 = 0, 1 otherwise */
/* NOTE: arg must be <= 5*2^32 + 4 */
u_int8_t not_divisible_by_5( u_int64_t arg )
NOT_DIVISIBLE_5

/* First, throw out a candidate if the maximum possible continuation
 * of the chain including it cannot reach the smallest target prime.
 *
 * for steps_to_go > 0, check that neither of the divisibility theorems
 * apply to the chain (including the candidate). Returns validity flag = 0
 * if either theorem applies, otherwise returns validity flag = 1
 *
 * For steps_to_go == 0, check if the candidate == any target prime.
 * Of course, candidates with c % 2 = 0 or gcd(c, chain_values[0]) > 1
 * can be bypassed immediately. Always returns 0 since we are done with
 * this particular chain */
void check_candidate(void)
CHECK_CAN_TEMPLATE_1(0)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3( max_val )) && (not_divisible_by_5( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain();
CHECK_CAN_TEMPLATE_9

/* find all candidates to extend the working chain that pass
 * the filters in check_candidate(). If steps_to_go == 1, valid candidates
 * will be compared to the list of target primes and all successes will be recorded */
u_int16_t gen_candidate_list(void)
GEN_C_LIST_1(0)
	doubles_count = extract_chain_values();
GEN_C_LIST_2
			check_candidate(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3(c)) && (not_divisible_by_5(c)) )) )
GEN_C_LIST_3A
							check_candidate(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(0)


#if (MAX_THREADS > 1)
void generate_and_process_candidate_list_01(void)
GEN_AND_PROCESS_C_LIST_1(1)
	c_count[r_level] = gen_candidate_list_01();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_01();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_01(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_01(void)
COPY_C_TO_W_CHAIN(1)

u_int8_t extract_chain_values_01(void)
EXTRACT_CHAIN_VALUES(1)

u_int64_t encode_Lchain_01(void)
ENCODE_LCHAIN_TEMPLATE(1)

u_int8_t not_divisible_by_3_01( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_01( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_01(void)
CHECK_CAN_TEMPLATE_1(1)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_01( max_val )) && (not_divisible_by_5_01( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_01();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_01();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_01();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_01();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_01(void)
GEN_C_LIST_1(1)
	doubles_count = extract_chain_values_01();
GEN_C_LIST_2
			check_candidate_01(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_01(c)) && (not_divisible_by_5_01(c)) )) )
GEN_C_LIST_3A
							check_candidate_01(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_01( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(1)
#endif


#if (MAX_THREADS > 2)
void generate_and_process_candidate_list_02(void)
GEN_AND_PROCESS_C_LIST_1(2)
	c_count[r_level] = gen_candidate_list_02();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_02();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_02(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_02(void)
COPY_C_TO_W_CHAIN(2)

u_int8_t extract_chain_values_02(void)
EXTRACT_CHAIN_VALUES(2)

u_int64_t encode_Lchain_02(void)
ENCODE_LCHAIN_TEMPLATE(2)

u_int8_t not_divisible_by_3_02( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_02( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_02(void)
CHECK_CAN_TEMPLATE_1(2)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_02( max_val )) && (not_divisible_by_5_02( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_02();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_02();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_02();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_02();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_02(void)
GEN_C_LIST_1(2)
	doubles_count = extract_chain_values_02();
GEN_C_LIST_2
			check_candidate_02(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_02(c)) && (not_divisible_by_5_02(c)) )) )
GEN_C_LIST_3A
							check_candidate_02(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_02( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(2)
#endif


#if (MAX_THREADS > 3)
void generate_and_process_candidate_list_03(void)
GEN_AND_PROCESS_C_LIST_1(3)
	c_count[r_level] = gen_candidate_list_03();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_03();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_03(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_03(void)
COPY_C_TO_W_CHAIN(3)

u_int8_t extract_chain_values_03(void)
EXTRACT_CHAIN_VALUES(3)

u_int64_t encode_Lchain_03(void)
ENCODE_LCHAIN_TEMPLATE(3)

u_int8_t not_divisible_by_3_03( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_03( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_03(void)
CHECK_CAN_TEMPLATE_1(3)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_03( max_val )) && (not_divisible_by_5_03( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_03();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_03();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_03();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_03();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_03(void)
GEN_C_LIST_1(3)
	doubles_count = extract_chain_values_03();
GEN_C_LIST_2
			check_candidate_03(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_03(c)) && (not_divisible_by_5_03(c)) )) )
GEN_C_LIST_3A
							check_candidate_03(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_03( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(3)
#endif


#if (MAX_THREADS > 4)
void generate_and_process_candidate_list_04(void)
GEN_AND_PROCESS_C_LIST_1(4)
	c_count[r_level] = gen_candidate_list_04();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_04();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_04(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_04(void)
COPY_C_TO_W_CHAIN(4)

u_int8_t extract_chain_values_04(void)
EXTRACT_CHAIN_VALUES(4)

u_int64_t encode_Lchain_04(void)
ENCODE_LCHAIN_TEMPLATE(4)

u_int8_t not_divisible_by_3_04( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_04( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_04(void)
CHECK_CAN_TEMPLATE_1(4)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_04( max_val )) && (not_divisible_by_5_04( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_04();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_04();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_04();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_04();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_04(void)
GEN_C_LIST_1(4)
	doubles_count = extract_chain_values_04();
GEN_C_LIST_2
			check_candidate_04(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_04(c)) && (not_divisible_by_5_04(c)) )) )
GEN_C_LIST_3A
							check_candidate_04(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_04( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(4)
#endif


#if (MAX_THREADS > 5)
void generate_and_process_candidate_list_05(void)
GEN_AND_PROCESS_C_LIST_1(5)
	c_count[r_level] = gen_candidate_list_05();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_05();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_05(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_05(void)
COPY_C_TO_W_CHAIN(5)

u_int8_t extract_chain_values_05(void)
EXTRACT_CHAIN_VALUES(5)

u_int64_t encode_Lchain_05(void)
ENCODE_LCHAIN_TEMPLATE(5)

u_int8_t not_divisible_by_3_05( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_05( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_05(void)
CHECK_CAN_TEMPLATE_1(5)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_05( max_val )) && (not_divisible_by_5_05( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_05();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_05();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_05();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_05();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_05(void)
GEN_C_LIST_1(5)
	doubles_count = extract_chain_values_05();
GEN_C_LIST_2
			check_candidate_05(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_05(c)) && (not_divisible_by_5_05(c)) )) )
GEN_C_LIST_3A
							check_candidate_05(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_05( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(5)
#endif


#if (MAX_THREADS > 6)
void generate_and_process_candidate_list_06(void)
GEN_AND_PROCESS_C_LIST_1(6)
	c_count[r_level] = gen_candidate_list_06();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_06();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_06(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_06(void)
COPY_C_TO_W_CHAIN(6)

u_int8_t extract_chain_values_06(void)
EXTRACT_CHAIN_VALUES(6)

u_int64_t encode_Lchain_06(void)
ENCODE_LCHAIN_TEMPLATE(6)

u_int8_t not_divisible_by_3_06( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_06( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_06(void)
CHECK_CAN_TEMPLATE_1(6)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_06( max_val )) && (not_divisible_by_5_06( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_06();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_06();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_06();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_06();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_06(void)
GEN_C_LIST_1(6)
	doubles_count = extract_chain_values_06();
GEN_C_LIST_2
			check_candidate_06(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_06(c)) && (not_divisible_by_5_06(c)) )) )
GEN_C_LIST_3A
							check_candidate_06(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_06( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(6)
#endif


#if (MAX_THREADS > 7)
void generate_and_process_candidate_list_07(void)
GEN_AND_PROCESS_C_LIST_1(7)
	c_count[r_level] = gen_candidate_list_07();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_07();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_07(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_07(void)
COPY_C_TO_W_CHAIN(7)

u_int8_t extract_chain_values_07(void)
EXTRACT_CHAIN_VALUES(7)

u_int64_t encode_Lchain_07(void)
ENCODE_LCHAIN_TEMPLATE(7)

u_int8_t not_divisible_by_3_07( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_07( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_07(void)
CHECK_CAN_TEMPLATE_1(7)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_07( max_val )) && (not_divisible_by_5_07( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_07();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_07();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_07();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_07();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_07(void)
GEN_C_LIST_1(7)
	doubles_count = extract_chain_values_07();
GEN_C_LIST_2
			check_candidate_07(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_07(c)) && (not_divisible_by_5_07(c)) )) )
GEN_C_LIST_3A
							check_candidate_07(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_07( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(7)
#endif


#if (MAX_THREADS > 8)
void generate_and_process_candidate_list_08(void)
GEN_AND_PROCESS_C_LIST_1(8)
	c_count[r_level] = gen_candidate_list_08();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_08();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_08(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_08(void)
COPY_C_TO_W_CHAIN(8)

u_int8_t extract_chain_values_08(void)
EXTRACT_CHAIN_VALUES(8)

u_int64_t encode_Lchain_08(void)
ENCODE_LCHAIN_TEMPLATE(8)

u_int8_t not_divisible_by_3_08( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_08( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_08(void)
CHECK_CAN_TEMPLATE_1(8)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_08( max_val )) && (not_divisible_by_5_08( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_08();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_08();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_08();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_08();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_08(void)
GEN_C_LIST_1(8)
	doubles_count = extract_chain_values_08();
GEN_C_LIST_2
			check_candidate_08(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_08(c)) && (not_divisible_by_5_08(c)) )) )
GEN_C_LIST_3A
							check_candidate_08(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_08( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(8)
#endif


#if (MAX_THREADS > 9)
void generate_and_process_candidate_list_09(void)
GEN_AND_PROCESS_C_LIST_1(9)
	c_count[r_level] = gen_candidate_list_09();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_09();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_09(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_09(void)
COPY_C_TO_W_CHAIN(9)

u_int8_t extract_chain_values_09(void)
EXTRACT_CHAIN_VALUES(9)

u_int64_t encode_Lchain_09(void)
ENCODE_LCHAIN_TEMPLATE(9)

u_int8_t not_divisible_by_3_09( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_09( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_09(void)
CHECK_CAN_TEMPLATE_1(9)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_09( max_val )) && (not_divisible_by_5_09( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_09();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_09();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_09();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_09();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_09(void)
GEN_C_LIST_1(9)
	doubles_count = extract_chain_values_09();
GEN_C_LIST_2
			check_candidate_09(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_09(c)) && (not_divisible_by_5_09(c)) )) )
GEN_C_LIST_3A
							check_candidate_09(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_09( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(9)
#endif


#if (MAX_THREADS > 10)
void generate_and_process_candidate_list_10(void)
GEN_AND_PROCESS_C_LIST_1(10)
	c_count[r_level] = gen_candidate_list_10();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_10();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_10(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_10(void)
COPY_C_TO_W_CHAIN(10)

u_int8_t extract_chain_values_10(void)
EXTRACT_CHAIN_VALUES(10)

u_int64_t encode_Lchain_10(void)
ENCODE_LCHAIN_TEMPLATE(10)

u_int8_t not_divisible_by_3_10( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_10( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_10(void)
CHECK_CAN_TEMPLATE_1(10)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_10( max_val )) && (not_divisible_by_5_10( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_10();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_10();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_10();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_10();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_10(void)
GEN_C_LIST_1(10)
	doubles_count = extract_chain_values_10();
GEN_C_LIST_2
			check_candidate_10(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_10(c)) && (not_divisible_by_5_10(c)) )) )
GEN_C_LIST_3A
							check_candidate_10(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_10( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(10)
#endif


#if (MAX_THREADS > 11)
void generate_and_process_candidate_list_11(void)
GEN_AND_PROCESS_C_LIST_1(11)
	c_count[r_level] = gen_candidate_list_11();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_11();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_11(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_11(void)
COPY_C_TO_W_CHAIN(11)

u_int8_t extract_chain_values_11(void)
EXTRACT_CHAIN_VALUES(11)

u_int64_t encode_Lchain_11(void)
ENCODE_LCHAIN_TEMPLATE(11)

u_int8_t not_divisible_by_3_11( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_11( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_11(void)
CHECK_CAN_TEMPLATE_1(11)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_11( max_val )) && (not_divisible_by_5_11( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_11();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_11();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_11();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_11();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_11(void)
GEN_C_LIST_1(11)
	doubles_count = extract_chain_values_11();
GEN_C_LIST_2
			check_candidate_11(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_11(c)) && (not_divisible_by_5_11(c)) )) )
GEN_C_LIST_3A
							check_candidate_11(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_11( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(11)
#endif


#if (MAX_THREADS > 12)
void generate_and_process_candidate_list_12(void)
GEN_AND_PROCESS_C_LIST_1(12)
	c_count[r_level] = gen_candidate_list_12();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_12();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_12(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_12(void)
COPY_C_TO_W_CHAIN(12)

u_int8_t extract_chain_values_12(void)
EXTRACT_CHAIN_VALUES(12)

u_int64_t encode_Lchain_12(void)
ENCODE_LCHAIN_TEMPLATE(12)

u_int8_t not_divisible_by_3_12( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_12( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_12(void)
CHECK_CAN_TEMPLATE_1(12)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_12( max_val )) && (not_divisible_by_5_12( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_12();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_12();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_12();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_12();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_12(void)
GEN_C_LIST_1(12)
	doubles_count = extract_chain_values_12();
GEN_C_LIST_2
			check_candidate_12(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_12(c)) && (not_divisible_by_5_12(c)) )) )
GEN_C_LIST_3A
							check_candidate_12(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_12( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(12)
#endif


#if (MAX_THREADS > 13)
void generate_and_process_candidate_list_13(void)
GEN_AND_PROCESS_C_LIST_1(13)
	c_count[r_level] = gen_candidate_list_13();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_13();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_13(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_13(void)
COPY_C_TO_W_CHAIN(13)

u_int8_t extract_chain_values_13(void)
EXTRACT_CHAIN_VALUES(13)

u_int64_t encode_Lchain_13(void)
ENCODE_LCHAIN_TEMPLATE(13)

u_int8_t not_divisible_by_3_13( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_13( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_13(void)
CHECK_CAN_TEMPLATE_1(13)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_13( max_val )) && (not_divisible_by_5_13( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_13();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_13();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_13();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_13();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_13(void)
GEN_C_LIST_1(13)
	doubles_count = extract_chain_values_13();
GEN_C_LIST_2
			check_candidate_13(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_13(c)) && (not_divisible_by_5_13(c)) )) )
GEN_C_LIST_3A
							check_candidate_13(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_13( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(13)
#endif


#if (MAX_THREADS > 14)
void generate_and_process_candidate_list_14(void)
GEN_AND_PROCESS_C_LIST_1(14)
	c_count[r_level] = gen_candidate_list_14();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_14();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_14(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_14(void)
COPY_C_TO_W_CHAIN(14)

u_int8_t extract_chain_values_14(void)
EXTRACT_CHAIN_VALUES(14)

u_int64_t encode_Lchain_14(void)
ENCODE_LCHAIN_TEMPLATE(14)

u_int8_t not_divisible_by_3_14( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_14( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_14(void)
CHECK_CAN_TEMPLATE_1(14)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_14( max_val )) && (not_divisible_by_5_14( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_14();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_14();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_14();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_14();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_14(void)
GEN_C_LIST_1(14)
	doubles_count = extract_chain_values_14();
GEN_C_LIST_2
			check_candidate_14(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_14(c)) && (not_divisible_by_5_14(c)) )) )
GEN_C_LIST_3A
							check_candidate_14(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_14( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(14)
#endif


#if (MAX_THREADS > 15)
void generate_and_process_candidate_list_15(void)
GEN_AND_PROCESS_C_LIST_1(15)
	c_count[r_level] = gen_candidate_list_15();
GEN_AND_PROCESS_C_LIST_2
	copy_candidate_to_working_chain_15();
GEN_AND_PROCESS_C_LIST_3
		generate_and_process_candidate_list_15(); /* recursive call */
GEN_AND_PROCESS_C_LIST_4

void copy_candidate_to_working_chain_15(void)
COPY_C_TO_W_CHAIN(15)

u_int8_t extract_chain_values_15(void)
EXTRACT_CHAIN_VALUES(15)

u_int64_t encode_Lchain_15(void)
ENCODE_LCHAIN_TEMPLATE(15)

u_int8_t not_divisible_by_3_15( u_int64_t arg )
NOT_DIVISIBLE_3

u_int8_t not_divisible_by_5_15( u_int64_t arg )
NOT_DIVISIBLE_5

void check_candidate_15(void)
CHECK_CAN_TEMPLATE_1(15)
			if( (max_c_flag == 1) && (max_val & 1) && (not_divisible_by_3_15( max_val )) && (not_divisible_by_5_15( max_val )) )
			{
				/* check if max_val is on the target list */
				if( ( max_val % 7 ) && ( max_val % 11 ) && ( max_val % 13 ) )
CHECK_CAN_TEMPLATE_4
										chain_code_list[code_index] = encode_Lchain_15();
CHECK_CAN_TEMPLATE_5
												temp_var = encode_Lchain_15();
CHECK_CAN_TEMPLATE_6
								chain_code_list[code_index] = encode_Lchain_15();
CHECK_CAN_TEMPLATE_8
									temp_var = encode_Lchain_15();
CHECK_CAN_TEMPLATE_9

u_int16_t gen_candidate_list_15(void)
GEN_C_LIST_1(15)
	doubles_count = extract_chain_values_15();
GEN_C_LIST_2
			check_candidate_15(); /* results written to check_result[] array */
GEN_C_LIST_3
						if( c_flag && ( (steps_to_go > 1) || ( (c & 0x1) && (not_divisible_by_3_15(c)) && (not_divisible_by_5_15(c)) )) )
GEN_C_LIST_3A
							check_candidate_15(); /* results written to check_result[] array */
GEN_C_LIST_4

void copy_work_assignment_to_thread_15( u_int8_t wrk_indx )
COPY_WORK_TO_THREAD(15)
#endif


void init_thread_memory(void)
{
	u_int8_t i;
	u_int32_t k, lim;

	lim = tgt_prime_list[thread_mem[0].tgt_p_count - 1].save_index - tgt_prime_list[0].save_index;

	for( i = 1; i < thread_count; i++ )
	{
		thread_mem[i].tgt_p_count = thread_mem[0].tgt_p_count;
		thread_mem[i].index_count_per_val = thread_mem[0].index_count_per_val;
		thread_mem[i].chain_code_list_start_index = thread_mem[0].chain_code_list_start_index;
		thread_mem[i].w_chain_length = thread_mem[0].w_chain_length;

		for(k = 0; k < thread_mem[0].tgt_p_count; k++)
		{
			thread_mem[i].chain_count[k] = 0;
			thread_mem[i].chain_max_dbl_count[k] = 0;
			thread_mem[i].chain_count_max_dbls[k] = 0;
			thread_mem[i].tgt_prime_code_length[k] = 0;
		}

		for(k = 0; k <= lim; k++)
			thread_mem[i].chain_code_list[k] = 0;
	}
}

void consolidate_results(void)
{
	u_int32_t i, k, chain_sum;
	u_int32_t code_index, max_dbls_thrd_indx;
	u_int16_t max_dbls_chain_sum;
	u_int8_t max_dbls, min_code_length;

	/* find and store best chain for each target prime */
	for( i = 0; i < thread_mem[0].tgt_p_count; i++ )
	{
		/* add up total chains found */
		chain_sum = thread_mem[0].chain_count[i];
		for( k = 1; k < thread_count; k++ )
			chain_sum += thread_mem[k].chain_count[i];

		if( chain_sum > 0 )
		{
			thread_mem[0].chain_count[i] = chain_sum;

			/* find thread with maximum # of doubled elements */
			max_dbls = thread_mem[0].chain_max_dbl_count[i]; /* max # of doubled elements in a chain */
			max_dbls_thrd_indx = 0;
			min_code_length = thread_mem[0].tgt_prime_code_length[i];
			max_dbls_chain_sum = thread_mem[0].chain_count_max_dbls[i]; /* number of chains with max doubles */
			for( k = 1; k < thread_count; k++ )
			{
				if( thread_mem[k].chain_count[i] > 0 )
				{
					if( thread_mem[k].chain_max_dbl_count[i] >= max_dbls )
					{
						if( thread_mem[k].chain_max_dbl_count[i] > max_dbls )
						{
							max_dbls = thread_mem[k].chain_max_dbl_count[i]; /* max # of doubled elements in a chain */
							max_dbls_thrd_indx = k;
							min_code_length = thread_mem[k].tgt_prime_code_length[i];
							max_dbls_chain_sum = thread_mem[k].chain_count_max_dbls[i]; /* number of chains with max doubles */
						}
						else /* doubled element counts are equal */
						{
							max_dbls_chain_sum += thread_mem[k].chain_count_max_dbls[i];
							if( thread_mem[k].tgt_prime_code_length[i] < min_code_length )
							{
								max_dbls_thrd_indx = k;
								min_code_length = thread_mem[k].tgt_prime_code_length[i];
							}
						}
					}
				}
			}
			thread_mem[0].chain_count_max_dbls[i] = max_dbls_chain_sum;
			if( max_dbls_thrd_indx > 0 )
			{
				thread_mem[0].chain_max_dbl_count[i] = max_dbls;
				thread_mem[0].tgt_prime_code_length[i] = min_code_length;
				/* move chain code to thread_mem[0].chain_code_list */
				code_index = tgt_prime_list[i].save_index - thread_mem[0].chain_code_list_start_index;
				thread_mem[0].chain_code_list[code_index] = thread_mem[max_dbls_thrd_indx].chain_code_list[code_index];
			}
		}
	}
}

int32_t main( int argc, char *argv[])
{
	u_int8_t *dif_table;
	u_int8_t *sieve_space;
	u_int32_t dif_table_start_index, sieve_prime_count, p_count, total_p_count;
	u_int32_t sieve_space_start_index;
	u_int32_t indx, dif_index;
	u_int64_t true_indx, max_indx_value, max_odd_val, i64;

/*	target_prime *tgt_prime_list; */
	u_int64_t *chain_code_list, clock_start, clock_stop, temp_var;
	u_int32_t *chain_code_list_start_index;
	u_int32_t code_save_index, last_save_index = 0xFFFFFFFF;
	u_int32_t j, last_j;
	u_int32_t *tgt_p_count, *chain_count;
	u_int64_t *Fib, *Luc;
	int32_t i;
	u_int32_t k;
#if 0
	u_int16_t exception_count, exception_index, prime_exception_count;
	u_int64_t excp_1_val, exception_list_1_step[40];
	u_int8_t on_list_flag;
	u_int32_t unique_chain_count;
#endif
	u_int16_t *c_list_start_index;
	u_int8_t *current_partial_length, *w_chain_length;
	u_int16_t *chain_count_max_dbls;
	u_int8_t *chain_max_dbl_count;
	u_int64_t total_prime_chain_count, interval_chain_count;
	u_int64_t total_chain_count_max_dbls;
	u_int8_t test_length, min_test_length, max_test_length, test_length_restart;
	u_int8_t restart_flag, new_length_init, final_interval = 0;
    u_int8_t truncating_for_B1, gen_exit_flag, reached_last_prime = 0;
	u_int32_t code_save_count;
	double c_value_range, *index_count_per_val;

	thread_io_struct thrd_io[MAX_THREADS];
	int rc;
	pthread_t tid[MAX_THREADS];

	double B1_in, thread_count_in;
	u_int64_t B1;
	u_int32_t old_tgt_prime_list_count, old_pending_code_list_count;
	u_int32_t old_tgt_p_file_read_count, old_pending_code_file_read_count;
	u_int32_t new_tgt_prime_list_count, new_pending_code_list_count;
	u_int32_t chain_list_zero_count, old_smallest_unsaved_code_index;
	u_int32_t smallest_unsaved_code_index;
	u_int64_t smallest_target_prime_next_list;
	u_int32_t pending_list_read_count, pending_list_zero_count;
	u_int32_t chain_code_array_space_remaining, chain_code_array_count;
	u_int64_t largest_target_prime_next_list = 0;
	u_int8_t max_code_length;
	u_int8_t *tgt_prime_code_length;
	u_int32_t code_index;
	FILE *chain_code_file;
	FILE *old_tgt_p_list_read_file, *new_tgt_p_list_write_file;
	FILE *old_pending_code_list_read_file, *new_pending_code_list_write_file;
	FILE *current_status_file;
	static const char *file_1 = "Lchain_codes.dat";
	static const char *file_2 = "Pending_Lchain_code_list_1.dat";
	static const char *file_3 = "Tgt_prime_list_1.dat";
	static const char *file_4 = "Pending_Lchain_code_list_2.dat";
	static const char *file_5 = "Tgt_prime_list_2.dat";
	static const char *file_6 = "current_status.dat";
	static const char *old_target_prime_filename, *new_target_prime_filename;
	static const char *old_pending_code_filename, *new_pending_code_filename;
	size_t dum;

	B1_in = 0;
	thread_count_in = -1;
    if( argc < 2 ) /* no arguments? */
    {
      printf("Upper limit B1 required for Lucas chain generator!\n"
              "Example: LucasChainGen -B1 3e6\n");
      return -1;
    }

	/* get upper limit for target primes & number of threads */
    while ((argc > 1) && (argv[1][0] == '-'))
    {
    	if( strcmp( argv[1], "-B1") == 0)
    	{
    		B1_in = strtod (argv[2], &argv[2]);
    	}
    	else if( strcmp( argv[1], "-nT") == 0)
    	{
    		thread_count_in = strtod (argv[2], &argv[2]);
    	}
    	argv+=2;
	    argc-=2;
   }

    if( B1_in <= 0 ) /* bad syntax */
    {
      printf("ERROR: Upper limit B1 > 0 required for Lucas chain generator!\n"
              "Example: ./LucasChainGen -B1 3e6\n");
      return -1;
    }

    B1 = (u_int64_t)B1_in;
    printf("\nGenerator upper limit B1 = %lu\n\n", B1);

    if( thread_count_in > 0 )
    {
        if( (thread_count_in < 1) || (thread_count_in > 16) ) /* bad syntax */
        {
          printf("ERROR: number of threads must be an integer, 1 <= nT <= 16. Default number of threads is 4.\n"
        		  "Example: ./LucasChainGen -B1 3e6 -nT 8");
          return -1;
        }
        thread_count = (u_int8_t)thread_count_in;
    	printf("Number of threads set to %u\n\n", thread_count);
    }
    else
    {
    	thread_count = DEFAULT_THREAD_COUNT;
    	printf("Number of threads set to %u (default)\n\n", thread_count);
    }

    /* keep the compiler happy */
    old_tgt_p_list_read_file = (FILE *)NULL;
    old_pending_code_list_read_file = (FILE *)NULL;

	/* initialize pointers & arrays */
	dif_table = get_dif_table_ptr();
	sieve_space = get_sieve_space_ptr();;
	Fib = thread_mem[0].Fib;
	Luc = thread_mem[0].Luc;
	current_partial_length = &thread_mem[0].current_partial_length;
	c_list_start_index = &thread_mem[0].c_list_start_index;
	chain_count = thread_mem[0].chain_count;
	chain_count_max_dbls = thread_mem[0].chain_count_max_dbls;
	chain_max_dbl_count = thread_mem[0].chain_max_dbl_count;
	w_chain_length = &thread_mem[0].w_chain_length;
	tgt_p_count = &thread_mem[0].tgt_p_count;
/*	tgt_prime_list = thread_mem[0].tgt_prime_list; */
	chain_code_list = thread_mem[0].chain_code_list;
	chain_code_list_start_index = &thread_mem[0].chain_code_list_start_index;
	index_count_per_val = &thread_mem[0].index_count_per_val;
	tgt_prime_code_length = thread_mem[0].tgt_prime_code_length;

	if( thread_count > 1 )
	{
		/* create threads */
		for( i = 0; i < thread_count; i++)
		{
			work_request[i] = NULL;
			thrd_io[i].thrd_indx = (u_int8_t)i;

			rc = pthread_create(&tid[i], NULL, recursive_work, (void *)&thrd_io[i]);
			if (rc)
			{
				printf("Error:unable to create thread, %d", rc);
				return -1;
			}
		}
	}

	/* initialize the difference table and primes used for sieving */
	sieve_prime_count = sieve_init();

	/* for reference: pi(10^3) = 168, pi(10^4) = 1229, pi(10^5) = 9592
	pi(10^6) = 78498, pi(10^7) = 664579, pi(10^8) = 5761455,
	pi(10^9) = 50847534, pi(10^10) = 455052511, pi(2*10^10) = 882206716
	 */

	/* sieve and then count primes */
	standard_sieve( sieve_prime_count );
	/* start count with 17 */
	sieve_space_start_index = 8;
	dif_table_start_index = 1;
	p_count = prime_count( &sieve_space_start_index, &dif_table_start_index );
	total_p_count = p_count + 6; /* all multiples of 2, 3, 5, 7, 11, and 13 avoided */
	printf("pi(10^6) = %d, known value = 78498\n", total_p_count);

	/* clear the chain code list */
	for(i = 0;i < MAX_CODE_OR_PRIME_COUNT; i++)
		chain_code_list[i] = 0;

	truncating_for_B1 = _false_;
	gen_exit_flag = _false_;
	*current_partial_length = 2;
	*c_list_start_index = 0;

	init_working_chains();
	set_work_assignments();
	init_Fib_sequence();
	init_Luc_sequence();

	restart_flag = _false_; /* if restarting, set to _true_ and recompile */
/*	restart_flag = _true_; */ /* if restarting, set to _true_ and recompile */

	if( restart_flag == _false_ ) /* start from scratch */
	{
		indx = 8; /* target prime list starts at 17 = (2*8 + 1) */
		true_indx = 8;
		dif_index = 1;

		/* code_save_index = (n - 5) for the nth prime, assuming p(1) = 2, p(2) = 3, p(3) = 5, etc. */
		/* chain code list indexes for primes 11 and 13 are 0 and 1, respectively */
		code_save_index = 2; /* index for p(7) = 17 */
		*tgt_p_count = 0;
		max_code_length = 3;

		chain_code_list[0] = CHAIN_START_4_7;
		chain_code_list[1] = CHAIN_START_5_8_13;

		chain_code_file = fopen(file_1,"w");
		fwrite((int8_t *)chain_code_list, sizeof(u_int64_t), 2, chain_code_file);
		fclose(chain_code_file);
		printf("\nChain codes for 11 and 13 written to file %s\n", file_1);

		*chain_code_list_start_index = 2;
		old_smallest_unsaved_code_index = 2;
		smallest_unsaved_code_index = 2;

		min_test_length = 6; /* bypass small primes 5 7 11 and 13 */
		new_length_init = _true_;

		old_tgt_prime_list_count = 0;
		old_tgt_p_file_read_count = 0;
		old_pending_code_list_count = 0;
		old_pending_code_file_read_count = 0;
		new_tgt_prime_list_count = 0;
		new_pending_code_list_count = 0;
		chain_list_zero_count = 0;
		total_prime_chain_count = 0;
		total_chain_count_max_dbls = 0;
	}
	else /* resume from where we left off */
	{
		current_status_file = fopen(file_6,"r");
		dum = fread((int8_t *)chain_code_list_start_index, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&smallest_unsaved_code_index, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&old_smallest_unsaved_code_index, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&smallest_target_prime_next_list, sizeof(u_int64_t), 1, current_status_file);
		dum = fread((int8_t *)&old_tgt_prime_list_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&old_tgt_p_file_read_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&old_pending_code_list_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&old_pending_code_file_read_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&new_tgt_prime_list_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&new_pending_code_list_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&chain_list_zero_count, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&total_prime_chain_count, sizeof(u_int64_t), 1, current_status_file);
		dum = fread((int8_t *)&total_chain_count_max_dbls, sizeof(u_int64_t), 1, current_status_file);
		dum = fread((int8_t *)&indx, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&code_save_index, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&dif_index, sizeof(u_int32_t), 1, current_status_file);
		dum = fread((int8_t *)&true_indx, sizeof(u_int64_t), 1, current_status_file);
		dum = fread((int8_t *)&test_length_restart, sizeof(u_int8_t), 1, current_status_file);
		dum = fread((int8_t *)&new_length_init, sizeof(u_int8_t), 1, current_status_file);
		dum = fread((int8_t *)&max_code_length, sizeof(u_int8_t), 1, current_status_file);
		if(dum == 0)
			printf("ERROR: Size mismatch reading status file!\n");
		fclose(current_status_file);

		/* update sieve space to current interval */
		i64 = SIEVE_SPACE_SIZE - 1; /* maximum sieve space index */
		while(i64 < true_indx)
		{
			standard_sieve( sieve_prime_count );
			i64 += SIEVE_SPACE_SIZE;
		}

		min_test_length = test_length_restart;
	}

	/* set a maximum chain length.The value of B1 will determine
	 * what the final length will actually be */
	max_test_length = 50;

	for(test_length = min_test_length; test_length <= max_test_length; test_length++)
	{
		/* set file name pointers. Files are "ping-ponged" as test_length is incremented */
		if( (test_length & 1) != 0 )
		{
			old_target_prime_filename = file_3;
			new_target_prime_filename = file_5;
			old_pending_code_filename = file_2;
			new_pending_code_filename = file_4;
		}
		else
		{
			old_target_prime_filename = file_5;
			new_target_prime_filename = file_3;
			old_pending_code_filename = file_4;
			new_pending_code_filename = file_2;
		}

		/* open active input files and, if restarting, set file positions */
		if( old_tgt_prime_list_count > 0 )
		{
			old_tgt_p_list_read_file = fopen(old_target_prime_filename, "r");

			/* bypass any target primes which have already been processed */
			if( old_tgt_p_file_read_count > 0 )
				fseek( old_tgt_p_list_read_file, (long int)(old_tgt_p_file_read_count*sizeof(target_prime)), SEEK_SET );
		}

		if( old_pending_code_list_count > 0 )
		{
			old_pending_code_list_read_file = fopen(old_pending_code_filename, "r");

			/* bypass any chain codes which have already been processed */
			if( old_pending_code_file_read_count > 0 )
				fseek( old_pending_code_list_read_file, (long int)(old_pending_code_file_read_count*sizeof(u_int64_t)), SEEK_SET );
		}

		/* set maximum possible sieve index value for a chain of length == test_length */
        max_odd_val = Fib[test_length + 2];
        if( max_odd_val > B1 ) /* don't exceed B1 */
        {
        	truncating_for_B1 = _true_;
        	max_odd_val = B1;
        }
        if((max_odd_val & 1) == 0)
          max_odd_val -= 1;
        max_indx_value = (max_odd_val - 1)/2;

    	if( true_indx <= max_indx_value )
    	{
    		reached_last_prime = _false_;
    		last_save_index = 0xFFFFFFFF;
    	}

    	clock_start = cputime();

		final_interval = _false_;

		/* Note that longer chain lengths will require processing multiple subintervals to
		 * avoid overwhelming memory limits and array bounds */
		do
		{
			/* Note: chain_code_list_start_index is the true base index
			 * for the chain code array. It is used to adjust each
			 * target prime's save_index to within the chain code array bounds */

			if( new_length_init == _true_)
			{
				printf("\n\nStarting generation of chains of length = %u\n\n", test_length);
			}
			else /* current length has already started */
			{
				printf("\nResuming generation of chains of length = %u\n\n", test_length);
			}

			/* read in saved chain code array from previous test length */
			if( old_pending_code_list_count > 0 )
			{
				/* read records remaining in the pending chain code file, up to the array size */
				pending_list_read_count = old_pending_code_list_count - old_pending_code_file_read_count;
				if( pending_list_read_count > MAX_CODE_OR_PRIME_COUNT )
					pending_list_read_count = MAX_CODE_OR_PRIME_COUNT;
				dum = fread((int8_t *)chain_code_list, sizeof(u_int64_t), pending_list_read_count, old_pending_code_list_read_file);
				if(dum == 0)
					printf("ERROR: End-Of-File reading old pending list file!\n");
				old_pending_code_file_read_count += pending_list_read_count;
				if( old_pending_code_list_count == old_pending_code_file_read_count )
				{
					/* no more of the old list to process */
					old_pending_code_list_count = 0;
					old_pending_code_file_read_count = 0;
					fclose(old_pending_code_list_read_file);
					printf("info: old pending code file closed after interval start\n");
				}
			}
			else
				pending_list_read_count = 0;

			if( pending_list_read_count > 0 )
			{
				/* the old pending chain code file is processed after chain generation so that the
				 * records just read in above always start with a zero chain code, so there is
				 * always at least one target prime corresponding to chain_code_list[0] */
				if( chain_code_list[0] != 0 ) /* should never happen */
					printf("ERROR: chain_code_list[0] != 0 \n");

				/* count zeros in the chain code array just read in */
				pending_list_zero_count = 0;
				for( i = 0; i < (int32_t)pending_list_read_count; i++ )
					if( chain_code_list[i] == 0 )
						pending_list_zero_count++;

				/* # of records remaining in the old target prime list file */
				i = (int32_t)(old_tgt_prime_list_count - old_tgt_p_file_read_count);
				if( pending_list_zero_count > (u_int32_t)i ) /* sanity check - should never happen */
				{
					printf("ERROR: code array zero count %u > remaining records in old target prime file %d!\n",
							pending_list_zero_count, i);
					return -1;
				}

				/* read in target primes corresponding to the zeros in the code array so far */
				dum = fread((int8_t *)tgt_prime_list, sizeof(target_prime), pending_list_zero_count, old_tgt_p_list_read_file);
				if(dum == 0)
					printf("ERROR: End-Of-File reading old target prime list file!\n");
				old_tgt_p_file_read_count += pending_list_zero_count;
				if( old_tgt_prime_list_count == old_tgt_p_file_read_count )
				{
					/* no more of the old target prime list to process */
					old_tgt_prime_list_count = 0;
					old_tgt_p_file_read_count = 0;
					fclose(old_tgt_p_list_read_file);
				}

				*tgt_p_count = pending_list_zero_count;

				/* logic check - remove after testing */
				*chain_code_list_start_index = tgt_prime_list[0].save_index;
				for( i = 0; i < (int32_t)pending_list_zero_count; i++ )
				{
					k = tgt_prime_list[i].save_index - *chain_code_list_start_index;
					if( k >= pending_list_read_count )
						printf("ERROR: tgt_prime_list[%d].save_index >= current code array count\n", i);
					if( chain_code_list[k] != 0 )
						printf("ERROR: chain_code_list[%u] zero mismatch with tgt_prime_list[%d].save_index\n", k, i);
				}
				if( (old_pending_code_list_count == 0) && (old_tgt_prime_list_count != 0) )
					printf("ERROR: old pending list total zero count < old target prime total count\n");
				/* end logic check */
			}
			else
			{
				*tgt_p_count = 0;
				pending_list_zero_count = 0;
			}

			chain_code_array_space_remaining = MAX_CODE_OR_PRIME_COUNT - pending_list_read_count;
			if( chain_code_array_space_remaining > 0 )
			{
				 /* logic check - remove after testing */
				if( old_tgt_prime_list_count > 0 )
					printf("ERROR: %u unread records in old_tgt_p_list_read_file\n", (old_tgt_prime_list_count - old_tgt_p_file_read_count));
				/* end logic check */

				/* extract primes from the sieve until either the corresponding chain code array
				 * is full or the maximum prime <= min( B1, Fib[test_length + 2]) is reached */
				while( (chain_code_array_space_remaining > 0) && (true_indx <= max_indx_value) )
				{
					if( sieve_space[indx] != 0 )
					{
						tgt_prime_list[*tgt_p_count].prime = 2*true_indx + 1;
						tgt_prime_list[*tgt_p_count].save_index = code_save_index;
						*tgt_p_count += 1;
						chain_code_array_space_remaining--;
						code_save_index++;
					}
					indx += dif_table[dif_index];
					true_indx += dif_table[dif_index];
					dif_index++;
					if( dif_index == 5760 )
						dif_index = 0;
					if(indx >= SIEVE_SPACE_SIZE)
					{
						/* sieve the next interval */
						standard_sieve( sieve_prime_count );
						indx -= SIEVE_SPACE_SIZE;
					}
				}

				/* make sure that the sieve index is pointing to a prime for next time
				 * this fixes an edge case where the code array becomes full, leaving
				 * true_indx <= max_indx_value, but there are no more primes <= max_indx_value.
				 * We then get an interval with zero target primes. Which is bad. */
				while( sieve_space[indx] == 0 )
				{
					indx += dif_table[dif_index];
					true_indx += dif_table[dif_index];
					dif_index++;
					if( dif_index == 5760 )
						dif_index = 0;
					if(indx >= SIEVE_SPACE_SIZE)
					{
						/* sieve the next interval */
						standard_sieve( sieve_prime_count );
						indx -= SIEVE_SPACE_SIZE;
					}
				}

				if( !reached_last_prime && (true_indx > max_indx_value) )
				{
					last_save_index = tgt_prime_list[*tgt_p_count - 1].save_index;
					reached_last_prime = _true_;
					printf("last prime = %lu, last_save_index = %u\n", tgt_prime_list[*tgt_p_count - 1].prime, tgt_prime_list[*tgt_p_count - 1].save_index );
				}

				/* write zeros into the chain code array for the primes just generated */
				i = (int32_t)(*tgt_p_count - pending_list_zero_count);
				if( i > 0 )
				{
					k = pending_list_read_count + (u_int32_t)i;
					for(j = pending_list_read_count; j < k; j++)
						chain_code_list[ j ] = 0;
				}
			}

			/* "index_count_per_val" helps to give a rough estimate of where in the target prime list to start
			 * looking for a match between a final chain value and a prime in the target list */
			c_value_range = (double)(tgt_prime_list[*tgt_p_count - 1].prime - tgt_prime_list[0].prime);
			if( c_value_range > 0.0 )
				*index_count_per_val = ((double)(*tgt_p_count - 1))/c_value_range;
			else
				*index_count_per_val = 0.0;

			/* output status of current interval */
			chain_code_array_count = pending_list_read_count + (*tgt_p_count - pending_list_zero_count); /* total index count */
			printf("Info: this interval chain code array size = %u: # of target primes = %u, # of unsaved chain codes = %u\n",
					chain_code_array_count, *tgt_p_count, (pending_list_read_count - pending_list_zero_count));
			printf("Info: target prime range: smallest = %lu, largest = %lu\n",
					tgt_prime_list[0].prime, tgt_prime_list[*tgt_p_count - 1].prime );
			if( chain_code_array_count > MAX_CODE_OR_PRIME_COUNT ) /* sanity check - should never happen */
				printf("\nERROR: # of chain codes = %u > MAX_CODE_OR_PRIME_COUNT = %u\n", chain_code_array_count, MAX_CODE_OR_PRIME_COUNT);
			if( *tgt_p_count > MAX_CODE_OR_PRIME_COUNT ) /* sanity check - should never happen */
				printf("\nERROR: # of target primes = %u > MAX_CODE_OR_PRIME_COUNT = %u\n", *tgt_p_count, MAX_CODE_OR_PRIME_COUNT);

			/* initialize chain counts for the target primes */
			for(i = 0;i < (int32_t)(*tgt_p_count); i++)
			{
				chain_count[i] = 0;
				chain_max_dbl_count[i] = 0;
				chain_count_max_dbls[i] = 0;
				tgt_prime_code_length[i] = 0;
			}
#if 0
			if( tgt_prime_list[*tgt_p_count - 1].prime >= Luc[test_length] )
			{
				/* find any exception primes on the target list >= Luc[test_length]
				 * Note that 1-step exceptions for the chain [1,2,3] have the following symmetry:
				 * E(i) = E(L-3-i) for all i, 0 <= i < (L-3)/2 and L is the test length. Also,
				 * when L is odd with (L-3)/2 an integer, the middle term E((L-3)/2) is always composite,
				 * with E((L-3)/2) = 2*Fib[(L+1)/2]*Fib[(L+3)/2] */
				exception_list_1_step[0] = Luc[test_length];
				exception_count = 1;
				if( (test_length & 1) != 0 )
					j = (u_int32_t)((test_length - 3)/2);
				else
					j = (u_int32_t)((test_length - 2)/2);

				for(i=1;i < (int32_t)j;i++) /* note that min_test_length >= 6, so j >= 2 */
				{
					/* Starting with the length 2 chain [1,2,3], follow the fibonacci chain for i steps, (a(i) = Fib[i+2]),
					 *  take one step down (a(i + 1) = Fib[i+2] + Fib[i]),
					 * then the maximal continuation for (test_length - 3 - i) more steps */
					excp_1_val = Fib[i+4] + Fib[i+2];
					excp_1_val = Fib[test_length - 2 - i]*excp_1_val + Fib[test_length - 3 - i]*Fib[i+4];
					on_list_flag = 0;
					exception_index = 0;
					while(exception_index < exception_count)
					{
						if(excp_1_val == exception_list_1_step[exception_index])
						{
							on_list_flag = 1;
							break;
						}
						exception_index++;
					}
					if( (on_list_flag == 0) && ( (excp_1_val >= tgt_prime_list[0].prime) && (excp_1_val <= tgt_prime_list[*tgt_p_count - 1].prime) ) )
					{
						exception_list_1_step[exception_count] = excp_1_val;
						exception_count++;
					}
				}

				/* find & report prime exceptions >= Luc[test_length] */
				prime_exception_count = 0;
				for(i = 0;i < exception_count;i++)
				{
					temp_var = exception_list_1_step[i];
					if( ( (temp_var & 1) != 0 ) && (not_divisible_by_3(temp_var))
							&& (not_divisible_by_5(temp_var)) && ( temp_var % 7 ) && ( temp_var % 11 ) && ( temp_var % 13 ) )
					{
						j = 0;
						while( tgt_prime_list[j].prime < Luc[test_length] )
							j++;
						do
						{
							if( temp_var <= tgt_prime_list[j].prime)
							{
								if( temp_var == tgt_prime_list[j].prime)
								{
									prime_exception_count++;
									if( tgt_prime_list[j].prime == Luc[test_length] )
										printf("*** Lucas prime exception: Luc[%u] = %lu\n", test_length, Luc[test_length]);
									else
										printf("*** One-step prime exception > Luc[%u]: %lu\n", test_length, temp_var);
								}
								break;
							}
							j++;
						}
						while( j < *tgt_p_count );
					}
				}
			}

			if(tgt_prime_list[*tgt_p_count - 1].prime == Fib[test_length + 2])
			{
				printf("*** Fibonacci prime exception: Fib[%u] = %lu\n", (test_length + 2), Fib[test_length + 2] );
			}
#endif

			*w_chain_length = test_length;
			*chain_code_list_start_index = tgt_prime_list[0].save_index;

			/* sanity check - remove after testing */
			if( new_length_init == _true_)
			{
				if( smallest_unsaved_code_index != *chain_code_list_start_index )
					printf("ERROR: smallest_unsaved_code_index %u != tgt_prime_list[0].save_index %u\n",
							smallest_unsaved_code_index, tgt_prime_list[0].save_index);
			}
			new_length_init = _false_;

			printf("Processing chains...\n");
			/* this call will recursively generate any/all Lucas chains with length = test_length
				which can reach primes in the target list */


			if( (test_length <= 26) || (thread_count == 1) )
				generate_and_process_candidate_list();
			else
			{
				init_thread_memory();

				i = 0;
				do
				{
			        j = (u_int32_t)i;
			        for( k = 0; k < thread_count; k++ )
			        {
			        	if( work_request[k] != NULL)
			        	{
			        		*(work_request[k]) = i;
			        		work_request[k] = NULL;
			        		i++;
			        		if( i == TOTAL_WORK_COUNT )
			        			break;
			        	}
			        }

			        if( (u_int32_t)i != j )
			           	pthread_cond_broadcast(&my_cond);

			        if( i < TOTAL_WORK_COUNT )
			        	pthread_cond_wait(&my_cond, &my_cond_m);
				}
				while( i < TOTAL_WORK_COUNT );

				/* wait for all threads to finish */
		        while( num_waiting_threads < thread_count )
		            pthread_cond_wait(&my_cond, &my_cond_m);

				consolidate_results();
			}

			code_save_count = 0;
			i = 0;
			if( new_tgt_prime_list_count == 0 )
			{
				/* save chain codes to the final file until a target prime with no chain yet appears */

				/* find the first target prime with no chain yet */
				while( i < (int32_t)(*tgt_p_count) )
				{
					if(chain_count[i] == 0)
					{
						smallest_target_prime_next_list = tgt_prime_list[i].prime;
						break;
					}
					i++;
				}

				/* write chain codes to the final file */
				if( i < (int32_t)(*tgt_p_count) ) /* a target prime was found with chain_count == 0 */
				{
					code_save_count = tgt_prime_list[i].save_index - *chain_code_list_start_index;
					smallest_unsaved_code_index = tgt_prime_list[i].save_index;
				}
				else /* all target primes in the current interval received a chain code */
				{
					code_save_count = chain_code_array_count;
					smallest_unsaved_code_index += chain_code_array_count;
				}

				if( code_save_count > 0  )
				{
					chain_code_file = fopen(file_1,"a");
					fwrite((int8_t *)chain_code_list, sizeof(u_int64_t), code_save_count, chain_code_file);
					fclose(chain_code_file);
					printf("info: stored %u chain codes in file %s\n", code_save_count, file_1);
					/* sanity check: make sure all saved codes are nonzero */
					for( j = 0;j < code_save_count;j++ )
						if( chain_code_list[j] == 0 )
							printf("ERROR: chain code[%u] = 0 written to file %s!\n", j, file_1);
				}

				if( code_save_count < chain_code_array_count )
				{
					/* save the rest of the chain code array to the new pending chain code file */
					k = chain_code_array_count - code_save_count;
					new_pending_code_list_write_file = fopen(new_pending_code_filename, "w");
					fwrite((int8_t *)&chain_code_list[code_save_count], sizeof(u_int64_t), k, new_pending_code_list_write_file);
					printf("info: stored %u chain codes in file %s\n", k, new_pending_code_filename);
					new_pending_code_list_count = k;
					last_j = 0;

					if( old_pending_code_list_count > 0 )
					{
						/* read and save codes from the old pending code file until a zero code appears */
						k = 0;
						while( old_pending_code_file_read_count < old_pending_code_list_count )
						{
							dum = fread((int8_t *)&temp_var, sizeof(u_int64_t), 1, old_pending_code_list_read_file);
							if(dum == 0)
								printf("ERROR: Early EOF reading old pending list file!\n");
							if( temp_var != 0 )
							{
								fwrite((int8_t *)&temp_var, sizeof(u_int64_t), 1, new_pending_code_list_write_file);
								old_pending_code_file_read_count++;
								k++;
							}
							else
							{
								/* rewind one record */
								fseek(old_pending_code_list_read_file, -(long int)sizeof(u_int64_t), SEEK_CUR );
								break;
							}
						}

						if( k > 0 )
						{
							new_pending_code_list_count += k;
							printf("info: transferred %u chain codes from the old pending code file to file %s\n", k, new_pending_code_filename);
						}

						if( old_pending_code_file_read_count == old_pending_code_list_count )
						{
							/* no more of the old list to process */
							old_pending_code_list_count = 0;
							old_pending_code_file_read_count = 0;
							fclose(old_pending_code_list_read_file);
							printf("info: old pending code file closed after new prime list started\n");
						}
					}
					fclose(new_pending_code_list_write_file);

					/* save target primes with no chain yet to the new target prime file */
					new_tgt_p_list_write_file = fopen(new_target_prime_filename, "w");
					k = 0;
					for(j = (u_int32_t)i; j < *tgt_p_count; j++)
					{
						if( chain_count[j] == 0 )
						{
							fwrite((int8_t *)&tgt_prime_list[j].prime, sizeof(target_prime), 1, new_tgt_p_list_write_file);
							k++;
							last_j = j;
						}
					}
					fclose(new_tgt_p_list_write_file);
					new_tgt_prime_list_count = k;

					printf("info: saved %u primes in file %s, cumulative total = %u\n", k, new_target_prime_filename, new_tgt_prime_list_count);
					largest_target_prime_next_list = tgt_prime_list[last_j].prime;
					printf("This interval # of target primes receiving a chain code = %u\n", (*tgt_p_count - k));
				}
				else /* all target primes in the interval received a chain */
				{
					if( old_pending_code_list_count > 0 )
					{
						/* read and save codes from the old pending code file until a zero code appears */
						k = 0;
						chain_code_file = fopen(file_1,"a");
						while( old_pending_code_file_read_count < old_pending_code_list_count )
						{
							dum = fread((int8_t *)&temp_var, sizeof(u_int64_t), 1, old_pending_code_list_read_file);
							if(dum == 0)
								printf("ERROR: Early EOF reading old pending list file!\n");
							if( temp_var != 0 )
							{
								fwrite((int8_t *)&temp_var, sizeof(u_int64_t), 1, chain_code_file);
								old_pending_code_file_read_count++;
								k++;
							}
							else
							{
								/* rewind one record */
								fseek(old_pending_code_list_read_file, -(long int)sizeof(u_int64_t), SEEK_CUR );
								break;
							}
						}
						fclose(chain_code_file);

						if( k > 0 )
						{
							smallest_unsaved_code_index += k;
							printf("info: transferred %u chain codes from the old pending code file to file %s\n", k, file_1);
						}

						if( old_pending_code_file_read_count == old_pending_code_list_count )
						{
							/* no more of the old list to process */
							old_pending_code_list_count = 0;
							old_pending_code_file_read_count = 0;
							fclose(old_pending_code_list_read_file);
							printf("info: old pending code file closed after transfer to final chain code file\n");
						}
					}
				}
			}
			else /* new target prime list file has already started */
			{
				/* write the chain code array to the new pending code list file */
				new_pending_code_list_write_file = fopen(new_pending_code_filename, "a");
				fwrite((int8_t *)chain_code_list, sizeof(u_int64_t), chain_code_array_count, new_pending_code_list_write_file);
				new_pending_code_list_count += chain_code_array_count;
				printf("info: stored %u chain codes in file %s\n", chain_code_array_count, new_pending_code_filename);

				/* transfer any non-zero chain codes at the beginning of the unread portion of
				 * the old pending code file to the new pending code file */
				if( old_pending_code_list_count > 0 )
				{
					/* read and save codes from the old pending code file until a zero code appears */
					k = 0;
					while( old_pending_code_file_read_count < old_pending_code_list_count )
					{
						dum = fread((int8_t *)&temp_var, sizeof(u_int64_t), 1, old_pending_code_list_read_file);
						if(dum == 0)
							printf("ERROR: Early EOF reading old pending list file!\n");
						if( temp_var != 0 )
						{
							fwrite((int8_t *)&temp_var, sizeof(u_int64_t), 1, new_pending_code_list_write_file);
							old_pending_code_file_read_count++;
							new_pending_code_list_count++;
							k++;
						}
						else
						{
							/* rewind one record */
							fseek(old_pending_code_list_read_file, -(long int)sizeof(u_int64_t), SEEK_CUR );
							break;
						}
					}

					if( k > 0 )
						printf("info: transferred %u unsaved chain codes from old pending file to new pending file\n", k);

					if( old_pending_code_file_read_count == old_pending_code_list_count )
					{
						/* no more of the old list to process */
						old_pending_code_list_count = 0;
						old_pending_code_file_read_count = 0;
						fclose(old_pending_code_list_read_file);
						printf("info: old pending code file closed after transfer to new pending code file\n");
					}
				}
				fclose(new_pending_code_list_write_file);

				/* write target primes without a code yet to the new target prime file */
				new_tgt_p_list_write_file = fopen(new_target_prime_filename, "a");
				k = 0;
				for(j = 0; j < *tgt_p_count; j++)
				{
					if( chain_count[j] == 0 )
					{
						fwrite((int8_t *)&tgt_prime_list[j].prime, sizeof(target_prime), 1, new_tgt_p_list_write_file);
						k++;
						last_j = j;
					}
				}
				fclose(new_tgt_p_list_write_file);
				new_tgt_prime_list_count += k;

				if( k > 0 )
				{
					printf("Saved %u primes in file %s, cumulative total = %u\n", k, new_target_prime_filename, new_tgt_prime_list_count);
					largest_target_prime_next_list = tgt_prime_list[last_j].prime;
				}
				printf("This interval # of target primes receiving a chain code = %u\n", (*tgt_p_count - k));
			}

			/* count zeros in the chain code list to compare with # of primes in new tgt prime list*/
			for( k = code_save_count; k < chain_code_array_count; k++ )
				if( chain_code_list[k] == 0 )
					chain_list_zero_count++;

			if( chain_list_zero_count != new_tgt_prime_list_count )
				printf("\nERROR: next target prime count = %u does not match zero count in pending code list = %u\n",
						new_tgt_prime_list_count, chain_list_zero_count);

			interval_chain_count = 0;
			k = 0;
			for(i = 0;i < (int32_t)(*tgt_p_count); i++)
			{
				interval_chain_count += chain_count[i];
				k += chain_count_max_dbls[i];
			}

			printf("This interval # of chains leading to a target prime = %lu, # of max dbl chains = %u\n", interval_chain_count, k);

			total_prime_chain_count += interval_chain_count;
			total_chain_count_max_dbls += k;

/*			unique_chain_count = 0; */
			for(i = 0; i < (int32_t)(*tgt_p_count); i++)
			{
#if 0
				if( i < 10 ) /* print results for up to 10 primes at start of target list */
				{
					printf("target prime[%u] = %lu, total chain_count = %u, max # of doubles = %u, # of chains with max dbles = %u\n",
							i, tgt_prime_list[i].prime, chain_count[i], chain_max_dbl_count[i], chain_count_max_dbls[i]);
					printf("chain code for p: %016lX\n", chain_code_list[tgt_prime_list[i].save_index - *chain_code_list_start_index]);
				}
/* #endif */

				if( chain_count[i] == 1 )
				{
					printf("*** single unique chain for target prime[%u] = %lu, total chain_count = %u, max # of doubles = %u, # of chains with max dbles = %u\n",
							i, tgt_prime_list[i].prime, chain_count[i], chain_max_dbl_count[i], chain_count_max_dbls[i]);
					printf("*** chain code for p: %016lX\n", chain_code_list[tgt_prime_list[i].save_index - *chain_code_list_start_index]);
					unique_chain_count++;
				}

				if( (chain_count[i] > 1) && (chain_count_max_dbls[i] == 1) )
				{
					unique_chain_count++;
/* #if 0 */
					printf("*** unique chain for target prime[%u] = %lu, total chain_count = %u, max # of doubles = %u, # of chains with max dbles = %u\n",
							i, tgt_prime_list[i].prime, chain_count[i], chain_max_dbl_count[i], chain_count_max_dbls[i]);
					printf("*** chain code for p: %016lX\n", chain_code_list[tgt_prime_list[i].save_index - *chain_code_list_start_index]);
				}
#endif
				if( tgt_prime_code_length[i] > max_code_length )
				{
					max_code_length = tgt_prime_code_length[i];
					code_index = tgt_prime_list[i].save_index - *chain_code_list_start_index;
					printf("ALERT: new max code length = %u bits for prime = %lu, code = %016lX\n",
						max_code_length, tgt_prime_list[i].prime, chain_code_list[code_index]);
				}

				if( max_code_length > 64 )
					printf("BAD PROBLEM: maximum code length exceeds 64 bits. A program revision to handle longer codes is required.\n");
			}
#if 0
			if( unique_chain_count > 0 )
				printf("Total # of unique chains, including max double chains, = %u\n", unique_chain_count);
#endif

			if( (tgt_prime_list[0].save_index + MAX_CODE_OR_PRIME_COUNT) > last_save_index )
			{
				final_interval = _true_;
				printf("final_interval set to true\n");
			}

			if( final_interval )
			{
				printf("\nFor reference, Fib[%u] = %lu, Luc[%u] = %lu, Fib[%u] = %lu, Luc[%u] = %lu, Fib[%u] = %lu\n",
						test_length, Fib[test_length], test_length-1, Luc[test_length-1], test_length+1,
						Fib[test_length+1], test_length, Luc[test_length], test_length+2, Fib[test_length+2]);
				printf("For test_length = %u, total # of chains leading to a target prime = %lu, total # of max dbl chains = %lu\n",
						test_length, total_prime_chain_count, total_chain_count_max_dbls);

				if( chain_list_zero_count != new_tgt_prime_list_count )
					printf("\nERROR: target prime count = %u does not match zero count in pending code list = %u\n",
							new_tgt_prime_list_count, chain_list_zero_count);

				if(new_tgt_prime_list_count > 0)
				{
					if( B1 <= Fib[ test_length + 2 ] )
						printf("Total primes < B1 = %lu in next target list: %u. Smallest = %lu, largest = %lu\n",
								B1, new_tgt_prime_list_count, smallest_target_prime_next_list, largest_target_prime_next_list);
					else
						printf("Total primes < Fib[%u] in next target list: %u. Smallest = %lu, largest = %lu\n",
							(test_length + 2), new_tgt_prime_list_count, smallest_target_prime_next_list, largest_target_prime_next_list);
				}
				else
				{
					if( truncating_for_B1 )
						printf("No primes < B1 = %lu in next target prime list\n", B1);
					else
						printf("No primes < Fib[%u] in next target prime list\n", (test_length + 2));
				}

				if( new_pending_code_list_count > 0 )
					printf("Pending code list size = %u, total # of unsaved chain codes in list = %u\n",
							new_pending_code_list_count, (new_pending_code_list_count - chain_list_zero_count));
				else
					printf("Pending chain code list is empty.\n");

				printf("\n# of chain codes appended to file %s = %u\n", file_1, (smallest_unsaved_code_index - old_smallest_unsaved_code_index));
				printf("Cumulative total prime codes saved = %u\n", smallest_unsaved_code_index);

				clock_stop = cputime();
				printf("\nTime for test length = %u: %12.4f seconds\n\n", test_length, (double)(clock_stop - clock_start)/(double)1000000.0);

				test_length_restart = test_length + 1;
				old_tgt_prime_list_count = new_tgt_prime_list_count;
				old_pending_code_list_count = new_pending_code_list_count;
				new_tgt_prime_list_count = 0;
				new_pending_code_list_count = 0;
				chain_list_zero_count = 0;
				old_smallest_unsaved_code_index = smallest_unsaved_code_index;
				total_prime_chain_count = 0;
				total_chain_count_max_dbls = 0;
				new_length_init = _true_;

				if( truncating_for_B1 && (old_tgt_prime_list_count == 0) && (old_pending_code_list_count == 0) )
					gen_exit_flag = 1;
			}
			else
			{
				test_length_restart = test_length;
			}

			/* save current parameters for possible future restart */
			current_status_file = fopen(file_6,"w");
			fwrite((int8_t *)chain_code_list_start_index, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&smallest_unsaved_code_index, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&old_smallest_unsaved_code_index, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&smallest_target_prime_next_list, sizeof(u_int64_t), 1, current_status_file);
			fwrite((int8_t *)&old_tgt_prime_list_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&old_tgt_p_file_read_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&old_pending_code_list_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&old_pending_code_file_read_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&new_tgt_prime_list_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&new_pending_code_list_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&chain_list_zero_count, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&total_prime_chain_count, sizeof(u_int64_t), 1, current_status_file);
			fwrite((int8_t *)&total_chain_count_max_dbls, sizeof(u_int64_t), 1, current_status_file);
			fwrite((int8_t *)&indx, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&code_save_index, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&dif_index, sizeof(u_int32_t), 1, current_status_file);
			fwrite((int8_t *)&true_indx, sizeof(u_int64_t), 1, current_status_file);
			fwrite((int8_t *)&test_length_restart, sizeof(u_int8_t), 1, current_status_file);
			fwrite((int8_t *)&new_length_init, sizeof(u_int8_t), 1, current_status_file);
			fwrite((int8_t *)&max_code_length, sizeof(u_int8_t), 1, current_status_file);
			fclose(current_status_file);
		}
/*		while( true_indx <= max_indx_value); */
		while( !final_interval );
		if( gen_exit_flag )
		{
			dum = (size_t)system("rm -f Tgt_prime_list_1.dat");
			dum = (size_t)system("rm -f Tgt_prime_list_2.dat");
			dum = (size_t)system("rm -f Pending_Lchain_code_list_1.dat");
			dum = (size_t)system("rm -f Pending_Lchain_code_list_2.dat");
			dum = (size_t)system("rm -f current_status.dat");

			if( thread_count > 1 )
			{
				/* All threads will still be waiting, release them */
				all_work_done = 1;
				pthread_cond_broadcast(&my_cond);
			}
		    return EXIT_SUCCESS;
		}
	}

	return EXIT_SUCCESS;
}
