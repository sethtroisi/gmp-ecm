/*
 ============================================================================
 Name        : LucasChainGen.c
 Author      : P. B. McLaughlin
 Version     : 1.0
 Copyright   : Copyright 2023 by Philip McLaughlin. All rights reserved.
 Description : Minimal-length Lucas chain generator for prime numbers
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <time.h>
#include "LucasChainGen.h"

/* returns current clock count in microseconds */
u_int64_t cputime(void)
{
    return (u_int64_t)clock();
}

/* set up variables & variable arrays for multi-routine access */
chain_element *get_working_chain_ptr(void)
{
	static chain_element working_chain[MAX_WORKING_CHAIN_LENGTH];
	static int8_t init = 0;

	if(init == 0) /* populate the first three slots (always 1, 2, and 3) */
	{
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

		init = 1;
	}
	return working_chain;
}

u_int64_t *get_chain_values_ptr(void)
{
	static u_int64_t chain_values[MAX_WORKING_CHAIN_LENGTH];
	return chain_values;
}

chain_element *get_candidate_list_ptr(void)
{
	static chain_element candidate_list[MAX_CAND_LIST_COUNT];
	return candidate_list;
}

chain_element *get_raw_c_list_ptr(void)
{
	static chain_element raw_c_list[MAX_CANDIDATE_COUNT];
	return raw_c_list;
}

u_int64_t	*get_raw_c_counts_ptr(void)
{
	static u_int64_t raw_c_counts[MAX_CANDIDATE_COUNT];
	return raw_c_counts;
}

u_int8_t *get_check_result_ptr(void)
{
	static u_int8_t check_result[MAX_CANDIDATE_COUNT];
	return check_result;
}

u_int8_t *get_check_index_ptr(void)
{
	static u_int8_t check_index;
	return &check_index;
}



/* target primes for a length n chain, such that any prime in the list is <= Fib[n+2]
	and no chain of any shorter length leads to this prime */
target_prime *get_tgt_prime_list_ptr(void)
{
	static target_prime tgt_prime_list[MAX_CODE_OR_PRIME_COUNT];
	return tgt_prime_list;
}

/* # of primes in the target prime list */
u_int32_t *get_tgt_p_count_ptr(void)
{
	static u_int32_t tgt_p_count;
	return &tgt_p_count;
}

u_int64_t *get_chain_code_list_ptr(void)
{
	static u_int64_t chain_code_list[MAX_CODE_OR_PRIME_COUNT];
	return chain_code_list;
}

u_int32_t *get_chain_code_list_start_index_ptr(void)
{
	static u_int32_t chain_code_list_start_index;
	return &chain_code_list_start_index;
}

/* total # of different actual minimum length chains leading to a target prime */
u_int32_t *get_chain_count_ptr(void)
{
	static u_int32_t chain_count[MAX_CODE_OR_PRIME_COUNT];
	return chain_count;
}

/* the maximum # of doubled elements in any minimum length chain leading to a target prime */
u_int8_t *get_chain_max_dbl_count_ptr(void)
{
	static u_int8_t chain_max_dbl_count[MAX_CODE_OR_PRIME_COUNT];
	return chain_max_dbl_count;
}

/* # of different minimal length chains leading to a target prime
	with # of doubled elements = chain_max_dbl_count[] */
u_int16_t *get_chain_count_max_dbls_ptr(void)
{
	static u_int16_t chain_count_max_dbls[MAX_CODE_OR_PRIME_COUNT];
	return chain_count_max_dbls;
}

/* chain code bit counts for target primes */
u_int8_t *get_tgt_prime_code_length_ptr(void)
{
	static u_int8_t tgt_prime_code_length[MAX_CODE_OR_PRIME_COUNT];
	return tgt_prime_code_length;
}

u_int32_t	*get_code_length_problem_count_ptr(void)
{
	static u_int32_t code_length_problem_count;
	return &code_length_problem_count;
}

/* Fibonacci numbers */
u_int64_t *get_Fib_ptr(void)
{
	static u_int64_t Fib[FIB_LIMIT];
	static int8_t init = 0;

	if(init == 0)
	{
		int32_t i;
		Fib[0] = 0;
		Fib[1] = 1;
		for(i=2;i<FIB_LIMIT;i++)
		{
			Fib[i] = Fib[i-1] + Fib[i-2];
		}
		init = 1;
	}
	return Fib;
}

/* Lucas numbers */
u_int64_t *get_Luc_ptr(void)
{
	static u_int64_t Luc[FIB_LIMIT];
	static int8_t init = 0;

	if(init == 0)
	{
		int32_t i;
		Luc[0] = 2;
		Luc[1] = 1;
		for(i=2;i<FIB_LIMIT;i++)
		{
			Luc[i] = Luc[i-1] + Luc[i-2];
		}
		init = 1;
	}
	return Luc;
}

u_int8_t *get_w_chain_length_ptr(void)
{
	static u_int8_t w_chain_length;
	return &w_chain_length;
}

u_int8_t *get_current_partial_length_ptr(void)
{
	static u_int8_t current_partial_length;
	return &current_partial_length;
}

u_int16_t *get_c_list_start_index_ptr(void)
{
	static u_int16_t c_list_start_index;
	return &c_list_start_index;
}

u_int16_t *get_current_c_index_ptr(void)
{
	static u_int16_t current_c_index;
	return &current_c_index;
}

double *get_index_count_per_val_ptr(void)
{
	static double index_count_per_val;
	return &index_count_per_val;
}

u_int8_t *get_code_length_ptr(void)
{
	static u_int8_t code_length;
	return &code_length;
}

u_int8_t *get_max_code_length_ptr(void)
{
	static u_int8_t max_code_length;
	return &max_code_length;
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

#ifndef HAVE_ASM_GCD
/* simple gcd for small arguments */
u_int64_t gcd(u_int64_t arg1, u_int64_t arg2)
{
	u_int64_t max, min, rem;
	u_int32_t max32, min32, rem32;
	u_int16_t max16, min16, rem16;
	static u_int64_t lim32;
	static u_int32_t lim16, init = 0;

	if( init == 0 )
	{
		lim32 = 4294967296;
		lim16 = 65536;
		init = 1;
	}

	if(arg1 == arg2) return arg1;
	if(arg1 < arg2)
	{
		max = arg2;
		min = arg1;
	}
	else
	{
		max = arg1;
		min = arg2;
	}

	if( max >= lim32 )
	{
		do
		{
			rem = max % min;
			if( rem >= 2 )
			{
				if( min < lim16 )
				{
					max16 = (u_int16_t)min;
					min16 = (u_int16_t)rem;
					goto label_2;
				}
				if( min < lim32 )
				{
					max32 = (u_int32_t)min;
					min32 = (u_int32_t)rem;
					goto label_1;
				}
			}
			else
			{
				if(rem == 0)
					return min;
				else
					return 1;
			}
			max = min;
			min = rem;
		}
		while(1);
	}
	else if( max >= lim16 )
	{
		max32 = (u_int32_t)max;
		min32 = (u_int32_t)min;
		goto label_1;
	}
	else
	{
		max16 = (u_int16_t)max;
		min16 = (u_int16_t)min;
		goto label_2;
	}

label_1:
	do
	{
		rem32 = max32 % min32;
		if( rem32 >= 2 )
		{
			if( min32 < lim16 )
			{
				max16 = (u_int16_t)min32;
				min16 = (u_int16_t)rem32;
				goto label_2;
			}
		}
		else
		{
			if(rem32 == 0)
				return (u_int64_t)min32;
			else
				return 1;
		}
		max32 = min32;
		min32 = rem32;
	}
	while(1);

label_2:
	do
	{
		rem16 = max16 % min16;
		if( rem16 <= 1 )
		{
			if(rem16 == 0)
				return (u_int64_t)min16;
			else
				return 1;
		}
		max16 = min16;
		min16 = rem16;
	}
	while(1);
}
#endif

void copy_candidate_to_working_chain(void)
{
	static chain_element *working_chain, *candidate_list;
	static u_int8_t *current_partial_length;
	static u_int16_t *current_c_index;
	static u_int8_t init = 0;
	u_int16_t c_index;
	u_int8_t w_index;

	if(init == 0)
	{
		working_chain = get_working_chain_ptr();
		current_partial_length = get_current_partial_length_ptr();
		candidate_list = get_candidate_list_ptr();
		current_c_index = get_current_c_index_ptr();
		init = 1;
	}

	w_index = (u_int8_t)(*current_partial_length + 1);
	c_index = *current_c_index;

	working_chain[w_index].value           = candidate_list[c_index].value;
	working_chain[w_index].comp_offset_1   = candidate_list[c_index].comp_offset_1;
	working_chain[w_index].comp_offset_2   = candidate_list[c_index].comp_offset_2;
	working_chain[w_index].dif_offset      = candidate_list[c_index].dif_offset;
	working_chain[w_index].chain_dbl_count = candidate_list[c_index].chain_dbl_count;

	*current_partial_length = w_index;
}

u_int8_t extract_chain_values(void)
{
	static chain_element *working_chain;
	static u_int64_t *chain_values;
	static u_int8_t *current_partial_length;
	static u_int8_t init = 0;
	u_int8_t i, indx, double_count;

	if( init == 0 )
	{
		/* initialize array pointers */
		working_chain = get_working_chain_ptr();
		chain_values = get_chain_values_ptr();
		current_partial_length = get_current_partial_length_ptr();
		init = 1;
	}

	indx = *current_partial_length;
	i = 0;
	double_count = working_chain[ indx ].chain_dbl_count;
	do
	{
		chain_values[i] = working_chain[ indx ].value;
		i++;
		indx--;
	}
	while( indx > 0 );

	chain_values[i] = 1;
	return double_count;
}

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

	/* first 4 bits of code give the next three chain components */
	code_fragment = (u_int8_t)(chain_code & 0xF);
	chain_code >>= 4;
	switch( code_fragment )
	{
		case CHAIN_START_5_8_13:
		{
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

			start_frag_count[CHAIN_START_5_8_13]++;
			break;
		}
		case CHAIN_START_5_7_12:
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 2;
			Lchain[4].dif_offset = 1;

			Lchain[5].value = 12;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 3;

			start_frag_count[CHAIN_START_5_7_12]++;
			break;
		}
		case CHAIN_START_5_8_11:
		{
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

			start_frag_count[CHAIN_START_5_8_11]++;
			break;
		}
		case CHAIN_START_5_6_11:
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			Lchain[4].value = 6;
			Lchain[4].comp_offset_1 = 1;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 0;
			(*dbl_count)++;

			Lchain[5].value = 11;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 4;

			start_frag_count[CHAIN_START_5_6_11]++;
			break;
		}
		case CHAIN_START_4_7_11:
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

			Lchain[5].value = 11;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 2;

			start_frag_count[CHAIN_START_4_7_11]++;
			break;
		}
		case CHAIN_START_5_6_10:
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			Lchain[4].value = 6;
			Lchain[4].comp_offset_1 = 1;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 0;

			Lchain[5].value = 10;
			Lchain[5].comp_offset_1 = 1;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 0;
			*dbl_count += 2;

			start_frag_count[CHAIN_START_5_6_10]++;
			break;
		}
		case CHAIN_START_4_7_10:
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

			Lchain[5].value = 10;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 2;
			Lchain[5].dif_offset = 1;

			start_frag_count[CHAIN_START_4_7_10]++;
			break;
		}
		case CHAIN_START_5_6_9:
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			Lchain[4].value = 6;
			Lchain[4].comp_offset_1 = 1;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 0;
			(*dbl_count)++;

			Lchain[5].value = 9;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 2;
			Lchain[5].dif_offset = 2;

			start_frag_count[CHAIN_START_5_6_9]++;
			break;
		}
		case CHAIN_START_4_6_9:
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;

			Lchain[4].value = 6;
			Lchain[4].comp_offset_1 = 1;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 0;

			Lchain[5].value = 9;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 2;
			Lchain[5].dif_offset = 2;

			*dbl_count += 2;
			start_frag_count[CHAIN_START_4_6_9]++;
			break;
		}
		case CHAIN_START_4_5_9:
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;
			(*dbl_count)++;

			Lchain[4].value = 5;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 3;
			Lchain[4].dif_offset = 1;

			Lchain[5].value = 9;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 4;

			start_frag_count[CHAIN_START_4_5_9]++;
			break;
		}
		case CHAIN_START_4_7_8:
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 3;

			Lchain[5].value = 8;
			Lchain[5].comp_offset_1 = 1;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 0;

			*dbl_count += 2;
			start_frag_count[CHAIN_START_4_7_8]++;
			break;
		}
		case CHAIN_START_4_6_8:
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;

			Lchain[4].value = 6;
			Lchain[4].comp_offset_1 = 1;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 0;

			Lchain[5].value = 8;
			Lchain[5].comp_offset_1 = 1;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 0;

			*dbl_count += 3;
			start_frag_count[CHAIN_START_4_6_8]++;
			break;
		}
		case CHAIN_START_4_6_7:
		{
			Lchain[3].value = 4;
			Lchain[3].comp_offset_1 = 1;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 0;

			Lchain[4].value = 6;
			Lchain[4].comp_offset_1 = 1;
			Lchain[4].comp_offset_2 = 1;
			Lchain[4].dif_offset = 0;

			Lchain[5].value = 7;
			Lchain[5].comp_offset_1 = 1;
			Lchain[5].comp_offset_2 = 2;
			Lchain[5].dif_offset = 4;

			*dbl_count += 2;
			start_frag_count[CHAIN_START_4_6_7]++;
			break;
		}
		case CHAIN_START_5_8_10: /* extremely rare */
		{
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
			start_frag_count[CHAIN_START_5_8_10]++;
			break;
		}
		case CHAIN_START_5_7_9:  /* extremely rare */
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 2;
			Lchain[4].dif_offset = 1;

			Lchain[5].value = 9;
			Lchain[5].comp_offset_1 = 0;
			Lchain[5].comp_offset_2 = 3;
			Lchain[5].dif_offset = 1;

			start_frag_count[CHAIN_START_5_7_9]++;
			break;
		}
		case CHAIN_START_5_7_10:  /* extremely rare */
		{
			Lchain[3].value = 5;
			Lchain[3].comp_offset_1 = 0;
			Lchain[3].comp_offset_2 = 1;
			Lchain[3].dif_offset = 2;

			Lchain[4].value = 7;
			Lchain[4].comp_offset_1 = 0;
			Lchain[4].comp_offset_2 = 2;
			Lchain[4].dif_offset = 1;

			Lchain[5].value = 10;
			Lchain[5].comp_offset_1 = 1;
			Lchain[5].comp_offset_2 = 1;
			Lchain[5].dif_offset = 0;

			(*dbl_count)++;
			start_frag_count[CHAIN_START_5_7_10]++;
			break;
		}
		default:
		{
			printf("ERROR: bad chain code start value in generate_Lchain = %u\n", code_fragment);
			return 0;
		}
	}
	chain_length = 5;

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
					case 0: /* chain encoding did not require this step until chain length = 44 */
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
					default: /* cases 1, 2, & 3 reserved for any more rare steps which might occur */
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
0xkD, k > 0.  Reserved.
0xjE.  Type 5 with k = i – j - 3.
0xjF.  Type 6 with k = i – j - 2.
*/

u_int64_t encode_Lchain(void)
{
	static chain_element *working_chain, *raw_c_list;
	static u_int8_t *current_partial_length;
	static u_int8_t *check_index;
	static u_int8_t *code_length;
	static u_int8_t init = 0;

	u_int64_t chain_code;
	u_int64_t val3, val4, val5;
	u_int8_t i, j, len, index1, index2;
	u_int8_t step_type[MAX_WORKING_CHAIN_LENGTH];
	u_int8_t step_index2[MAX_WORKING_CHAIN_LENGTH];
	u_int8_t type_count, type_index, increment;

	if( init == 0 )
	{
		working_chain = get_working_chain_ptr();
		current_partial_length = get_current_partial_length_ptr();
		raw_c_list = get_raw_c_list_ptr();
		check_index = get_check_index_ptr();
		code_length = get_code_length_ptr();
		init = 1;
	}

	/* We assume that any application using Lucas chains for primes will
	 * handle the small primes 2, 3, 5, and 7 without using/requiring a chain code.
	 * Codes for primes 11 (CHAIN_START_4_7_11) and 13 (CHAIN_START_5_8_13)
	 * will be written directly to the beginning of the chain code save file */

	chain_code = 0;
	*code_length = 4;
	len = *current_partial_length + 1;

	if(len > 5)
	{
		/* set indexes for top element (raw_c_list[*check_index]) */
		index1 = raw_c_list[*check_index].comp_offset_1;
		index2 = raw_c_list[*check_index].comp_offset_2;

		/* ignore consecutive maximal steps at the top end (largest elements) of the chain */
		if( (index1 == 0) && (index2 == 1) )
		{
			do
			{
				len--;
				index1 = working_chain[len].comp_offset_1;
				index2 = working_chain[len].comp_offset_2;
			}
			while( (len > 5) && (index1 == 0) && (index2 == 1) );
		}

		if( len > 5 )
		{
			/* set step_type for each element with working chain index > 5 */
			type_count = 0;
			while( len > 5 )
			{
				if( (index1 == 0) && (index2 == 1) )
					step_type[type_count] = 1;
				else if( (index1 == 0) && (index2 == 2) )
					step_type[type_count] = 2;
				else if( (index1 == 1) && (index2 == 1) )
					step_type[type_count] = 3;
				else if( (index1 == 2) && (index2 == 2) )
					step_type[type_count] = 4;
				else if( (index1 == 0) && (index2 > 2) )
				{
					step_type[type_count] = 5;
					step_index2[type_count] = index2;
				}
				else if( (index1 == 1) && (index2 > 1) )
				{
					step_type[type_count] = 6;
					step_index2[type_count] = index2;
				}
				else
					printf("ERROR: unimplemented step_type: %u %u\n", index1, index2);

				type_count++;
				len--;
				index1 = working_chain[len].comp_offset_1;
				index2 = working_chain[len].comp_offset_2;
			}
			step_type[type_count] = 0;

			type_index = 0;
			do
			{
				switch( step_type[type_index] )
				{
					case 1:
						/* Note: during top-down encoding, a step of this type
						 * will always follow only step types 2 or 3, which both cover
						 * multiple subsequent maximal continuations (step type == 1) */
						printf("ERROR: isolated step type = 1.\n"); /* should never happen */
						return 0;

					case 2:
						j = 0;
						while( step_type[type_index + j + 1] == 1 )
							j++;

						increment = j + 1;

						i = 0;
						while( j >= 12 )
						{
							j -= 12;
							i++;
						}

						switch( j )
						{
							case 0:
								chain_code += 0x1;
								break;
							case 1:
								chain_code += 0x3;
								break;
							case 2:
								chain_code += 0x5;
								break;
							case 3:
								chain_code += 0x7;
								break;
							case 4:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x9;
								break;
							case 5:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x19;
								break;
							case 6:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x29;
								break;
							case 7:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x39;
								break;
							case 8:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0xB;
								break;
							case 9:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x1B;
								break;
							case 10:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x2B;
								break;
							case 11:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x3B;
								break;
							default: /* placeholder - should never get here */
								printf("ERROR: reached default #2\n");
								break;
						}
						switch( i ) /* step_type = 1, i*12 */
						{
							case 0:
								break;
							case 1:
								chain_code <<= 6;
								*code_length += 6;
								/* chain_code += 0x0;  12 maximal elements  */
								break;
							case 2:
								chain_code <<= 6;
								*code_length += 6;
								chain_code += 0x10; /* 24 maximal elements */
								break;
							case 3:
								chain_code <<= 6;
								*code_length += 6;
								chain_code += 0x20; /* 36 maximal elements */
								break;
							/* Note: the code fragment 0x30 is used for step type = 4, since a string of 48
							   maximal elements will never appear until chain lengths are well into the 50's */
							default: /* placeholder - should never get here */
								printf("ERROR: reached default #3\n");
								break;
						}
						break;

					case 3:
						j = 0;
						while( step_type[type_index + j + 1] == 1 )
							j++;

						increment = j + 1;

						i = 0;
						while( j >= 12 )
						{
							j -= 12;
							i++;
						}

						switch( j )
						{
							case 0:
								chain_code += 0x2;
								break;
							case 1:
								chain_code += 0x4;
								break;
							case 2:
								chain_code += 0x6;
								break;
							case 3:
								chain_code += 0x8;
								break;
							case 4:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0xA;
								break;
							case 5:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x1A;
								break;
							case 6:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x2A;
								break;
							case 7:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x3A;
								break;
							case 8:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0xC;
								break;
							case 9:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x1C;
								break;
							case 10:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x2C;
								break;
							case 11:
								chain_code <<= 2;
								*code_length += 2;
								chain_code += 0x3C;
								break;
							default: /* placeholder - should never get here */
								printf("ERROR: reached default #4\n");
								break;
						}
						switch( i ) /* step_type = 1, i*12 */
						{
						case 0:
							break;
						case 1:
							chain_code <<= 6;
							*code_length += 6;
							/* chain_code += 0x0;  12 maximal elements */
							break;
						case 2:
							chain_code <<= 6;
							*code_length += 6;
							chain_code += 0x10; /* 24 maximal elements */
							break;
						case 3:
							chain_code <<= 6;
							*code_length += 6;
							chain_code += 0x20; /* 36 maximal elements */
							break;
						default: /* placeholder - should never get here */
							printf("ERROR: reached default #5\n");
							break;
						/* Note: the code fragment 0x30 is used for step type = 4, since a string of 48
						   maximal elements will never appear until chain lengths are well into the 50's */
						}
						break;

					case 4:
						increment = 1;
						chain_code <<= 2;
						*code_length += 2;
						chain_code += 0x30;
						break;

					case 5:
						i = step_index2[type_index] - 3;
						chain_code <<= 2;
						*code_length += 2;
						increment = 1;
						switch(i)
						{
							case 0:
							case 1:
							case 2:
							case 3:
								chain_code += (u_int64_t)(i*16 + 0xE);
								break;
							case 4:
								chain_code += 0xD; /* not needed until chain length = 44 */
								break;
							default:
								printf("ERROR: no code for step, index1 = 0, index2 = %u\n", (i+3));
								return 0;
						}
						break;

					case 6:
						i = step_index2[type_index] - 2;
						if( i > 3 )
						{
							printf("ERROR: no code for step, index1 = 1, index2 = %u\n", (i+2));
							return 0;
						}
						chain_code <<= 2;
						*code_length += 2;
						chain_code += (u_int64_t)(i*16 + 0xF);
						increment = 1;
						break;

					default:
						printf("ERROR: bad step_type = %u\n", step_type[type_index]);
						return 0;
				} /* end switch */

				chain_code <<= 4;
				*code_length += 4;
				type_index += increment;
			}
			while( type_index < type_count );

		} /* end if len > 5 after top max continuation elements removed */
	} /* end if len > 5 */

	/*  special codes for the first three chain elements after 3 ([1,2,3] are implicit and not coded)
	 *  Easy inspection shows that there are 23 possible chain extensions which have total length == 5,
	 *  specifically [4,5,6], [4,5,7], [4,5,8], [4,5,9], [4,6,7], [4,6,8], [4,6,9], [4,6,10], [4,7,8],
	 *  [4,7,10], [4,7,11], [5,6,7], [5,6,8], [5,6,9], [5,6,10], [5,6,11], [5,7,8], [5,7,9], [5,7,10],
	 *  [5,7,12], [5,8,10], [5,8,11], and [5,8,13]. Of these, [4,6,10] only leads to composite integers.
	 *  [4,5,6] and [4,5,7] have not (for chain lengths <= 45) appeared in any minimal length chain, and
	 *  four of the remaining 20 ([4,5,8], [5,6,7], [5,6,8], and [5,7,8]
	 *  do appear in minimum length chains for primes, but never (for chain lengths <= 45) in a minimal
	 *  length chain with the maximum number of doubled elements. This suggests that a chain containing
	 *  one of these four will have a substitute sequence with an extra doubled element, but we
	 *  won't try to prove this. */
	switch( *current_partial_length )
	{
		case 2:
			val3 = raw_c_list[*check_index].value;
			val4 = val3 + 3;
			val5 = val4 + val3;
			break;

		case 3:
			val3 = working_chain[3].value;
			val4 = raw_c_list[*check_index].value;
			val5 = val4 + val3;
			break;

		case 4:
			val3 = working_chain[3].value;
			val4 = working_chain[4].value;
			val5 = raw_c_list[*check_index].value;
			break;

		default:
			val3 = working_chain[3].value;
			val4 = working_chain[4].value;
			val5 = working_chain[5].value;
	}

	switch( val3 )
	{
		case 4:
		{
			switch( val4 )
			{
				case 5:
				{
					switch( val5 )
					{
						case 9:
							chain_code += CHAIN_START_4_5_9;
							return chain_code;
						default:
							goto report_error;
					}
				}
				case 6:
				{
					switch( val5 )
					{
						case 7:
							chain_code += CHAIN_START_4_6_7;
							return chain_code;
						case 8:
							chain_code += CHAIN_START_4_6_8;
							return chain_code;
						case 9:
							chain_code += CHAIN_START_4_6_9;
							return chain_code;
						default:
							goto report_error;
					}
				}
				case 7:
				{
					switch( val5 )
					{
						case 8:
							chain_code += CHAIN_START_4_7_8;
							return chain_code;
						case 10:
							chain_code += CHAIN_START_4_7_10;
							return chain_code;
						case 11:
							chain_code += CHAIN_START_4_7_11;
							return chain_code;
						default:
							goto report_error;
					}
				}
				default:
					goto report_error;
			}
		} /* end case val3 == 4 */
		case 5:
		{
			switch( val4 )
			{
				case 6:
				{
					switch( val5 )
					{
						case 9:
							chain_code += CHAIN_START_5_6_9;
							return chain_code;
						case 10:
							chain_code += CHAIN_START_5_6_10;
							return chain_code;
						case 11:
							chain_code += CHAIN_START_5_6_11;
							return chain_code;
						default:
							goto report_error;
					}
				}
				case 7:
				{
					switch( val5 )
					{
						case 9:
							chain_code += CHAIN_START_5_7_9;
							return chain_code;
						case 10:
							chain_code += CHAIN_START_5_7_10;
							return chain_code;
						case 12:
							chain_code += CHAIN_START_5_7_12;
							return chain_code;
						default:
							goto report_error;
					}
				}
				case 8:
				{
					switch( val5 )
					{
						case 10:
							chain_code += CHAIN_START_5_8_10;
							return chain_code;
						case 11:
							chain_code += CHAIN_START_5_8_11;
							return chain_code;
						case 13:
							chain_code += CHAIN_START_5_8_13;
							return chain_code;
						default:
							goto report_error;
					}
				}
				default:
					goto report_error;
			}
		}
		default:
			goto report_error;
	}

	/* report unimplemented start sequence */
	report_error:
		printf("ERROR: unimplemented chain start sequence: %lu %lu %lu\n", val3, val4, val5);
		return 0;
}

void check_candidate(void)
{
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
	static chain_element *working_chain;
	static chain_element *raw_c_list;
	static target_prime *tgt_prime_list;
	static u_int8_t *check_result, *check_index;
	static u_int64_t *chain_code_list;
	static u_int32_t *chain_code_list_start_index;
	static u_int32_t *chain_count, *tgt_p_count;
	static u_int64_t *Fib;
	static u_int64_t *chain_values;
	static u_int16_t *chain_count_max_dbls;
	static u_int8_t *chain_max_dbl_count, *current_partial_length;
	static u_int8_t *code_length;
	static u_int8_t *tgt_prime_code_length;
	static u_int32_t *code_length_problem_count;
	static u_int8_t *w_chain_length, init = 0;
	static double *index_count_per_val;
	double approx_start_index, del;
	u_int64_t temp_var;
	u_int64_t gcd_c_p;
	u_int64_t c, max_val, p_val, compare_val;
	u_int32_t k, k_old, code_index;
	u_int64_t beta, upper_limit_1_step = 0, a1, a2 = 0;
	u_int8_t kk, steps_to_go, c_length, element_count, doubles_count;
	u_int8_t compare_val_in_chain, next_c_count, max_c_flag;
	u_int64_t next_step_c_list[MAX_CANDIDATE_COUNT];
	u_int64_t cand, dif;
	u_int8_t i, j, c_flag, ii;

	if( init == 0 )
	{
		working_chain = get_working_chain_ptr();
		chain_values = get_chain_values_ptr();
		raw_c_list = get_raw_c_list_ptr();
		check_index = get_check_index_ptr();
		check_result = get_check_result_ptr();
		tgt_prime_list = get_tgt_prime_list_ptr();
		chain_code_list = get_chain_code_list_ptr();
		chain_code_list_start_index = get_chain_code_list_start_index_ptr();
		chain_count = get_chain_count_ptr();
		chain_count_max_dbls = get_chain_count_max_dbls_ptr();
		chain_max_dbl_count = get_chain_max_dbl_count_ptr();
		tgt_p_count = get_tgt_p_count_ptr();
		current_partial_length = get_current_partial_length_ptr();
		code_length = get_code_length_ptr();
		tgt_prime_code_length = get_tgt_prime_code_length_ptr();
		code_length_problem_count = get_code_length_problem_count_ptr();
		w_chain_length = get_w_chain_length_ptr();
		Fib = get_Fib_ptr();
		index_count_per_val = get_index_count_per_val_ptr();
		init = 1;
	}

	c = raw_c_list[*check_index].value;
	steps_to_go = (u_int8_t)(*w_chain_length - *current_partial_length - 1);
	max_val = Fib[steps_to_go + 1]*c + Fib[steps_to_go]*chain_values[0];
	/* note that max_val == c when steps_to_go == 0 */

	if( (max_val < tgt_prime_list[0].prime) || (c > tgt_prime_list[*tgt_p_count - 1].prime) )
	{
		check_result[*check_index] = 0;
		return;
	}

	gcd_c_p = gcd( c, chain_values[0] );

	if( steps_to_go > 0 ) /* c is not the last (largest) element */
	{
		/* check for chains which never lead to primes */
		if( gcd_c_p > 1 )
		{
			/* check if candidate gives a chain of the form [...,x,...,2x,3x],
			 * with all of the elements between x and 2x less than 3x/2 */
			p_val = chain_values[0];
			if((working_chain[*current_partial_length].dif_offset == 0) /* doubled element = 2x */
				&& (c == 3*p_val/2) /* candidate = 3x */
				&& (2*chain_values[1] < 3*p_val/2 )) /* largest element below 2x is < 3x/2 */
			{
				/* check if 2x - e is an element in the chain for any element e such that x < e < 2x  */
				c_length = *current_partial_length;
				element_count = 2; /* find the index where chain_values[index] == p_val/2 */
				while( chain_values[element_count] > p_val/2 )
					element_count++;

				/* sanity check */
				if( 2*chain_values[element_count] != p_val )
					printf("ERROR: no chain element == p_val/2!\n");

				k = 1;
				compare_val_in_chain = 0;
				do
				{
					compare_val = p_val - chain_values[k];
					kk = c_length;
					do
					{
						if(compare_val == chain_values[kk] )
						{
							compare_val_in_chain = 1; /* a candidate to extend the chain exists which is NOT a multiple of p_val/2 */
							break;
						}
						kk--;
					}
					while( kk > element_count );
					if(compare_val_in_chain != 0) /* the 2x3x theorem does not apply */
						break;
					k++;
				}
				while( chain_values[k] > p_val/2 );

				if(compare_val_in_chain == 0) /* the 2x3x theorem applies; we can ignore the candidate */
				{
					check_result[*check_index] = 0;
					return;
				}

			} /* end check if 2x3x chain theorem applies */

			/* check if c, chain_values[0], and chain_values[1] share a nontrivial divisor */
			kk = 1;
			do
			{
				gcd_c_p = gcd( gcd_c_p, chain_values[kk] );
				kk++;
			}
			while(gcd_c_p > 1);
			kk--;

			if(kk >= 2)
			{
				 /* c, chain_values[0],..., chain_values[kk - 1]
				    share a nontrivial common divisor */
				compare_val = chain_values[kk]; /* gcd(compare_val, gcd_c_p) = 1 */

				/* If 2*compare_val > c, then it is a valid candidate to extend the chain including c.
					If 2*compare_val == c, then c + compare_val is a valid candidate, since
					c - compare_val = compare_val is in the chain */
				if(c <= 2*compare_val) /* "3-strikes" theorem does not apply */
					goto bypass;

				/* for each smaller element divisible by gcd_c_p, check that either
				   (element + compare_val) < c OR element > 2*compare_val (or both)
				   so that there are NO candidates > c involving compare_val, or
				   any chain element < compare_val. It is then an easy inductive proof
				   to show that all chain elements > c are divisible by gcd_c_p
				   in any chain containing [chain_values[...], c], and hence never prime */
				for(k=0;k < kk;k++)
				{
					if( ((chain_values[k] + compare_val) > c ) && (chain_values[k] <= 2*compare_val) ) /* might be a candidate */
					{
						 /* if chain_values[k] > 2*compare_val, then (chain_values[k] - compare_val) > compare_val.
						    Since the only elements in the chain greater than compare_val are divisible by gcd_c_p,
						    (chain_values[k] - compare_val) is never a chain member */
						goto bypass; /* theorem does not apply */
						/* note that chain_values[k] + compare_val != c, since equality would force
						   compare_val to be divisible by gcd_c_p, a contradiction */
					}
				}

				/* divisibility theorem applies; toss candidate out */
				check_result[*check_index] = 0;
				return;
			}
		} /* end if gcd_c_p > 1 */

bypass:
		/* check if the one-step upper limit lemma applies */

		/* get a list of all candidates for the next chain including c
		   which are greater than or equal to 2*chain_values[0] */
		next_step_c_list[0] = 2*chain_values[0];
		next_c_count = 1;

		j = 0;
		do
		{
			cand = c + chain_values[j];
			if( cand <= next_step_c_list[0] )
				break;

			dif = c - chain_values[j];
			if( dif > chain_values[0] )
				break;

			/* check if dif is an element in the chain */
			k = 0;
			do
			{
				if( dif >= chain_values[k] )
				{
					if( dif == chain_values[k] )
					{
						/* make sure that cand is not already in the candidate list */
						c_flag = 1;
						ii = 0;
						while( ii < next_c_count )
						{
							if( cand == next_step_c_list[ii] )
							{
								c_flag = 0;
								break;
							}
							ii++;
						}
						if( c_flag )
						{
							/* add cand to the list */
							next_step_c_list[next_c_count] = cand;
							next_c_count++;
						}
					}
					break;
				}
				k++;
			}
			while( 1 );
			j++;
		}
		while( j <= (*current_partial_length) );

		/* bubble sort the list to get largest to smallest */
		for( i = (next_c_count - 1); i > 0; i--)
		{
			for(j=0;j<i;j++)
			{
				if(next_step_c_list[j] < next_step_c_list[j+1])
				{
					dif = next_step_c_list[j];
					next_step_c_list[j] = next_step_c_list[j+1];
					next_step_c_list[j+1] = dif;
				}
			}
		}

		/* find beta (the largest candidate strictly less than the maximal element) */
		if( next_step_c_list[0] == ( c + chain_values[0] ) )
		{
			max_c_flag = 1; /* maximum chain continuation exists */
			beta = next_step_c_list[1];
		}
		else /* next_step_c_list[0] < (c + chain_values[0]) */
		{
			max_c_flag = 0; /* maximal chain continuation does not exist */
			beta = next_step_c_list[0];
		}

		/* determine the upper limit value that the current working chain can reach
		 *  in steps_to_go by choosing, at one of the steps, the largest candidate
		 *  less than (a(i) + a(i-1)), and the maximal candidate otherwise */
		a1 = Fib[steps_to_go]*beta + Fib[steps_to_go - 1]*c;
		if( (steps_to_go <= 2) || (max_c_flag == 0) )
			upper_limit_1_step = a1;
		 /* steps_to_go >= 3 && max_c_flag == 1 from here*/
		else if(steps_to_go == 3)
			a2 = 2*c + 3*chain_values[0];
		else if( (2*c <= 3*chain_values[0]) ) /* steps_to_go >= 4 */
			a2 = Fib[steps_to_go]*c + Fib[steps_to_go + 1]*chain_values[0];
		else
			a2 = (Fib[steps_to_go] + 2*Fib[steps_to_go - 3])*c + (Fib[steps_to_go - 1] + 2*Fib[steps_to_go - 4])*chain_values[0];

		if( (steps_to_go >= 3) && (max_c_flag == 1) )
		{
			if(a1 >= a2)
				upper_limit_1_step = a1;
			else
				upper_limit_1_step = a2;
		}

		if(upper_limit_1_step >= tgt_prime_list[0].prime)
		{
			check_result[*check_index] = 1;
			return;
		}
		else
		{
			if( max_val > tgt_prime_list[*tgt_p_count - 1].prime )
			{
				check_result[*check_index] = 0; /* an upcoming interval will catch this one */
				return;
			}


			if( (max_c_flag == 1) && ((max_val & 1) != 0) )
			{
				/* check if max_val is on the target list */
				if( (gcd_c_p == 1) && (gcd( max_val, 15015 ) == 1) )
				{
					doubles_count = raw_c_list[*check_index].chain_dbl_count;

					/* the code blitzes through the shorter lengths so fast
					 * that we don't need to adjust the start point for k */
					if( *w_chain_length > 15)
					{
						/* find approximate index to start comparison loop
						 	note that max_val >= tgt_prime_list[0].prime here, so approx_start_index >= 0 */
						approx_start_index = (double)( max_val - tgt_prime_list[0].prime )*( *index_count_per_val );

						k = (u_int32_t)approx_start_index;
						/* Note: since  tgt_prime_list[0].prime <= max_val <= Fib[test_length + 2],
						 * and *index_count_per_val = (*tgt_p_count - 1) / (Fib[test_length + 2] - tgt_prime_list[0].prime)
						 * in double precision floating point, then  0 <= k <= (*tgt_p_count - 1) */

						/* if max_val < tgt_prime_list[k].prime,
						 *  reduce k until max_val >= tgt_prime_list[k].prime */
						while( max_val < tgt_prime_list[k].prime ) /* if true, k must be > 0 */
						{
							del = ( (double)( tgt_prime_list[k].prime - max_val ) )*( *index_count_per_val ) + 1.0;
							approx_start_index -= del;
							if( approx_start_index < 0 ) approx_start_index = 0;
							k = (u_int32_t)approx_start_index;
						}

						del = ( (double)( max_val - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;
						k_old = k + (u_int32_t)del;

						while( (k_old < *tgt_p_count) && (max_val > tgt_prime_list[k_old].prime) )
						{
							k = k_old;
							del = ( (double)( max_val - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;
							k_old = k + (u_int32_t)del;
						}
					}
					else
						k = 0;

					while(k < *tgt_p_count)
					{
						if(max_val <= tgt_prime_list[k].prime)
						{
							if(max_val == tgt_prime_list[k].prime)
							{
								chain_count[k]++;
								if( doubles_count >= chain_max_dbl_count[k] )
								{
									if(doubles_count > chain_max_dbl_count[k])
									{
										/* best chain so far */
										chain_max_dbl_count[k] = doubles_count;
										chain_count_max_dbls[k] = 1;
										/* Encode & save this chain. Code will be overwritten
											if/when a better chain with more doubled elements is found */
										code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;
										chain_code_list[code_index] = encode_Lchain();

										if( (tgt_prime_code_length[k]) > 64 && (*code_length <= 64) )
										{
											printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",
													k, tgt_prime_list[k].prime);
											*code_length_problem_count -= 1;
										}

										tgt_prime_code_length[k] = *code_length;
										if(*code_length > 64)
										{
											printf("PROBLEM: chain code length > 64 bits for tgt_prime_list[%u] = %lu\n",
													k, tgt_prime_list[k].prime);
											*code_length_problem_count += 1;
										}
									}
									else
									{
										chain_count_max_dbls[k]++;
										if( tgt_prime_code_length[k] > 64 )
										{
											/* given this choice, see if a code exists with <= 64 bits */
											temp_var = encode_Lchain();
											if( *code_length <= 64 )
											{
												code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;
												chain_code_list[code_index] = temp_var;
												tgt_prime_code_length[k] = *code_length;
												printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",
														k, tgt_prime_list[k].prime);
												*code_length_problem_count -= 1;
											}
										}
									}
								}
							}
							break;
						}
						k++;
					}
				}
			}
			check_result[*check_index] = 0; /* this candidate has been resolved */
			return;
		}
	} /* end if steps_to_go > 0 */
	else /* steps_to_go == 0 */
	{
		if(  (gcd_c_p == 1) && (gcd( c, 15015 ) == 1) )
		{
			doubles_count = raw_c_list[*check_index].chain_dbl_count;

			/* the code blitzes through the shorter lengths so fast
			 * that we don't need to adjust the start point for k */
			if( *w_chain_length > 15)
			{
				/* find approximate index to start comparison loop
				 	note that c >= tgt_prime_list[0].prime here, so approx_start_index >= 0 */
				approx_start_index = (double)( c - tgt_prime_list[0].prime )*( *index_count_per_val );

				k = (u_int32_t)approx_start_index;

				/* if c < tgt_prime_list[k].prime,
				 *  reduce k until c >= tgt_prime_list[k].prime */
				while( c < tgt_prime_list[k].prime  ) /* if true, k must be > 0 */
				{
					del = ( (double)( tgt_prime_list[k].prime - c ) )*( *index_count_per_val ) + 1.0;
					approx_start_index -= del;
					if( approx_start_index < 0 ) approx_start_index = 0;
					k = (u_int32_t)approx_start_index;
				}

				del = ( (double)( c - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;
				k_old = k + (u_int32_t)del;

				while( (k_old < *tgt_p_count) && (c > tgt_prime_list[k_old].prime) )
				{
					k = k_old;
					del = ( (double)( c - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;
					k_old = k + (u_int32_t)del;
				}
			}
			else
				k = 0;

			while(k < *tgt_p_count)
			{
				if(c <= tgt_prime_list[k].prime)
				{
					if(c == tgt_prime_list[k].prime)
					{
						chain_count[k]++;
						if( doubles_count >= chain_max_dbl_count[k] )
						{
							if(doubles_count > chain_max_dbl_count[k])
							{
								/* best chain so far */
								chain_max_dbl_count[k] = doubles_count;
								chain_count_max_dbls[k] = 1;

								/* Encode & save this chain. Code will be overwritten
									if/when a better chain with more doubled elements is found */
								code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;
								chain_code_list[code_index] = encode_Lchain();

								if( (tgt_prime_code_length[k]) > 64 && (*code_length <= 64) )
								{
									printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",
											k, tgt_prime_list[k].prime);
									*code_length_problem_count -= 1;
								}

								tgt_prime_code_length[k] = *code_length;

								if(*code_length > 64)
								{
									printf("PROBLEM: chain code length > 64 bits for tgt_prime_list[%u] = %lu\n",
											k, tgt_prime_list[k].prime);
									*code_length_problem_count -= 1;
								}
							}
							else
							{
								chain_count_max_dbls[k]++;
								if( tgt_prime_code_length[k] > 64 )
								{
									/* given this choice, see if a code exists with <= 64 bits */
									temp_var = encode_Lchain();
									if( *code_length <= 64 )
									{
										code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;
										chain_code_list[code_index] = temp_var;
										tgt_prime_code_length[k] = *code_length;
										printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",
												k, tgt_prime_list[k].prime);
										*code_length_problem_count -= 1;
									}
								}
							}
						}
					}
					break;
				}
				k++;
			}
		}
		check_result[*check_index] = 0; /* this candidate has been resolved */
		return;
	} /* end steps_to_go == 0 */
}

/* find all candidates to extend the working chain that pass
 * the filters in check_candidate(). If steps_to_go == 1, valid candidates
 * will be compared to the list of target primes and all successes will be recorded */
u_int16_t gen_candidate_list(void)
{
	static chain_element *candidate_list;
	static chain_element *raw_c_list;
	static u_int64_t *chain_values;
	static u_int16_t *c_list_start_index;
	static u_int8_t *current_partial_length;
	static u_int8_t *w_chain_length, init = 0;
	static u_int8_t *check_result, *check_index;
	static u_int64_t *raw_c_counts;

	u_int64_t dif, c;
	u_int16_t c_index, c_count, ii;
	u_int8_t i, j, k, c_flag, doubles_count, steps_to_go, rcl_index;

	if( init == 0 )
	{
		/* initialize pointers */
		candidate_list = get_candidate_list_ptr();
		raw_c_list = get_raw_c_list_ptr();
		check_result = get_check_result_ptr();
		check_index = get_check_index_ptr();
		c_list_start_index = get_c_list_start_index_ptr();
		chain_values = get_chain_values_ptr();
		current_partial_length = get_current_partial_length_ptr();
		w_chain_length = get_w_chain_length_ptr();
		raw_c_counts = get_raw_c_counts_ptr();
		for( i = 0; i < 16; i++)
		{
			check_result[i] = 0;
			raw_c_counts[i] = 0;
		}

		init = 1;
	}

	steps_to_go = *w_chain_length - *current_partial_length;
	rcl_index = 0;

	/* extract the chain into the chain_values array. Return value
	 * is the # of doubled elements in the working chain */
	doubles_count = extract_chain_values();

	/* First, find all possible doubled candidates. It is easy to prove that
		there is always at least one, namely 2*chain_values[1] */
	if( steps_to_go > 1 ) /* if == 1, doubled candidates cannot be prime */
	{
		i = 1;
		while( (c = 2*chain_values[i]) > chain_values[0] )
		{
			raw_c_list[rcl_index].value = c;
			raw_c_list[rcl_index].comp_offset_1 = i;
			raw_c_list[rcl_index].comp_offset_2 = i;
			raw_c_list[rcl_index].dif_offset = 0; /* indicates a doubled element */
			raw_c_list[rcl_index].chain_dbl_count = doubles_count + 1;
			*check_index = rcl_index;
			check_candidate(); /* results written to check_result[] array */
			rcl_index++;
			i++;
		}
	}

	i = 0;
	j = 1;
	do
	{
		do
		{
			c = chain_values[i] + chain_values[j];
			if( c <= chain_values[0] )
				break;

			dif = chain_values[i] - chain_values[j];
			if( dif > chain_values[i + 1] )
				break;

			/* check if dif is an element in the chain */
			k = i + 1;
			do
			{
				if( dif >= chain_values[k] )
				{
					if( dif == chain_values[k] )
					{
						/* make sure that c is not already in the candidate list */
						c_flag = 1;
						ii = 0;
						while( ii < rcl_index )
						{
							if( c == raw_c_list[ii].value )
							{
								c_flag = 0;
								break;
							}
							ii++;
						}
						if( c_flag && ( (steps_to_go > 1) || (c & 0x1)) ) /* don't list if steps_to_go == 1 && c is even */
						{
							/* populate the next available slot. Done this way in case
							  the full chain, or a maximally extended partial chain,
							  requires encoding, which is done in the check_candidate() routine */
							raw_c_list[rcl_index].value = c;
							raw_c_list[rcl_index].comp_offset_1 = i;
							raw_c_list[rcl_index].comp_offset_2 = j;
							raw_c_list[rcl_index].dif_offset = k;
							raw_c_list[rcl_index].chain_dbl_count = doubles_count;
							*check_index = rcl_index;
							check_candidate(); /* results written to check_result[] array */
							rcl_index++;
							if( rcl_index == MAX_CANDIDATE_COUNT )
								printf("\nALERT: rcl_index reached maximum of %d - increase maximum\n\n", (int)MAX_CANDIDATE_COUNT );

						}
					}
					break;
				}
				k++;
			}
			while( 1 );
			j++;
		}
		while( j <= *current_partial_length );
		i++;
		j = i+1;
	}
	while( (chain_values[i] + chain_values[j]) > chain_values[0] );

/*	if( steps_to_go > 1 ) */
		raw_c_counts[rcl_index]++;
	c_index = *c_list_start_index;
	c_count = 0;
	for(i = 0; i < rcl_index; i++)
	{
		if( check_result[i] )
		{
			check_result[i] = 0;
			candidate_list[c_index].value = raw_c_list[i].value;
			candidate_list[c_index].comp_offset_1 = raw_c_list[i].comp_offset_1;
			candidate_list[c_index].comp_offset_2 = raw_c_list[i].comp_offset_2;
			candidate_list[c_index].dif_offset = raw_c_list[i].dif_offset;
			candidate_list[c_index].chain_dbl_count = raw_c_list[i].chain_dbl_count;
			c_index++;
			c_count++;
			if(c_index == MAX_CAND_LIST_COUNT)
			{
				printf("Candidate array is full. Increase MAX_CAND_LIST_COUNT\n");
				return c_count;
			}
		}
	}
	return c_count; /* total # of valid candidates to extend the current partial working chain */
}


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
{
	static u_int16_t *c_list_start_index, *current_c_index;
	static u_int16_t c_count[MAX_WORKING_CHAIN_LENGTH];
	static u_int16_t c_indx[MAX_WORKING_CHAIN_LENGTH];
	static u_int8_t *current_partial_length;
	static u_int8_t init = 0, r_level = 0;

	if( init == 0 )
	{
		/* initialize common variable & variable array pointers */
		current_partial_length = get_current_partial_length_ptr();
		c_list_start_index = get_c_list_start_index_ptr();
		current_c_index = get_current_c_index_ptr();
		init = 1;
	}

	/* find all candidates to be the next element of the working chain */
	/* final-step-in-the-chain candidates are processed in the following call */
	c_count[r_level] = gen_candidate_list();

	if(c_count[r_level] == 0) /* no candidates to process. Final step always returns zero */
		return;

	c_indx[r_level] = 0;
	do
	{
		*current_c_index = *c_list_start_index + c_indx[r_level];
		/* working_chain[*current_partial_length + 1] = candidate_list[*current_c_index] */
		/* copy candidate to working chain (increments *current_partial_length) */
		copy_candidate_to_working_chain();
		*c_list_start_index += c_count[r_level];
		r_level++;
		generate_and_process_candidate_list(); /* recursive call */
		r_level--;
		*c_list_start_index -= c_count[r_level];
		(*current_partial_length)--;

		c_indx[r_level]++;
	}
	while(c_indx[r_level] < c_count[r_level]);
}


int32_t main( int argc, char *argv[])
{
	u_int8_t *dif_table;
	u_int8_t *sieve_space;
	u_int32_t dif_table_start_index, sieve_prime_count, p_count, total_p_count;
	u_int32_t sieve_space_start_index;
	u_int32_t indx, dif_index;
	u_int64_t true_indx, max_indx_value, max_odd_val, i64;

	target_prime *tgt_prime_list;
	u_int64_t *chain_code_list, clock_start, clock_stop, temp_var;
	u_int32_t *chain_code_list_start_index;
	u_int32_t code_save_index;
	u_int32_t j, last_j;
	u_int32_t *tgt_p_count, *chain_count;
	u_int64_t *Fib, *Luc;
/*	u_int64_t a, b, c; */
	u_int64_t excp_1_val, exception_list_1_step[40];
	int32_t i;
	u_int32_t k;
	u_int16_t exception_count, exception_index, prime_exception_count;
	u_int16_t *c_list_start_index;
	u_int8_t *current_partial_length, *w_chain_length;
	u_int16_t *chain_count_max_dbls;
	u_int8_t *chain_max_dbl_count;
	u_int64_t total_prime_chain_count, interval_chain_count;
	u_int64_t total_chain_count_max_dbls;
	u_int8_t test_length, min_test_length, max_test_length, test_length_restart;
	u_int8_t on_list_flag, restart_flag, new_length_init;
    u_int8_t truncating_for_B1, gen_exit_flag;
	u_int32_t code_save_count;
	double c_value_range, *index_count_per_val;
	double B1_in;
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
	u_int8_t *max_code_length;
	u_int8_t *tgt_prime_code_length;
	u_int32_t *code_length_problem_count, code_index;
	u_int32_t unique_chain_count;
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
	static u_int64_t *raw_c_counts;
	size_t dum;


	B1_in = 0;
    if( argc < 2 ) /* no arguments? */
    {
      printf("Upper limit B1 required for Lucas chain generator!\n"
              "Example: LucasChainGen -B1 3e6\n");
      return -1;
    }

	/* get upper limit for target primes */
	while ((argc > 1) && (argv[1][0] == '-'))
    {
      if (strcmp (argv[1], "-B1") == 0)
        B1_in = strtod (argv[2], &argv[2]);
      argv+=2;
      argc-=2;
    }

    if( B1_in <= 0 ) /* bad syntax */
    {
      printf("Upper limit B1 > 0 required for Lucas chain generator!\n"
              "Example: LucasChainGen -B1 3e6\n");
      return -1;
    }

    B1 = (u_int64_t)B1_in;
    printf("\nGenerator upper limit B1 = %lu\n\n", B1);

    /* test gcd */
/*    b = 2570;
    a = 19937;
    c = gcd(a, b);
    printf("gcd(%lu, %lu) = %lu\n", a, b, c);
    return 0; */

    /* keep the compiler happy */
    old_tgt_p_list_read_file = (FILE *)NULL;
    old_pending_code_list_read_file = (FILE *)NULL;

    /* generate Fibonacci & Lucas numbers */
	Fib = get_Fib_ptr();
	Luc = get_Luc_ptr();

	/* initialize pointers & arrays */
	dif_table = get_dif_table_ptr();
	sieve_space = get_sieve_space_ptr();;
	current_partial_length = get_current_partial_length_ptr();
	c_list_start_index = get_c_list_start_index_ptr();
	chain_count = get_chain_count_ptr();
	chain_count_max_dbls = get_chain_count_max_dbls_ptr();
	chain_max_dbl_count = get_chain_max_dbl_count_ptr();
	w_chain_length = get_w_chain_length_ptr();
	tgt_p_count = get_tgt_p_count_ptr();
	tgt_prime_list = get_tgt_prime_list_ptr();
	chain_code_list = get_chain_code_list_ptr();
	chain_code_list_start_index = get_chain_code_list_start_index_ptr();
	index_count_per_val = get_index_count_per_val_ptr();
	max_code_length = get_max_code_length_ptr();
	tgt_prime_code_length = get_tgt_prime_code_length_ptr();
	code_length_problem_count = get_code_length_problem_count_ptr();
	raw_c_counts = get_raw_c_counts_ptr();

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

	*code_length_problem_count = 0;

	truncating_for_B1 = _false_;
	gen_exit_flag = _false_;
	*current_partial_length = 2;
	*c_list_start_index = 0;

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
		*max_code_length = 4;

		chain_code_list[0] = CHAIN_START_4_7_11;
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
		dum = fread((int8_t *)max_code_length, sizeof(u_int8_t), 1, current_status_file);
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

		clock_start = cputime();

		/* Note that longer chain lengths will require processing multiple subintervals to
		 * avoid overwhelming memory limits and array bounds */
		do
		{
			/* Note: chain_code_list_start_index is the true base index
			 * for the chain code array. It is used to adjust each
			 * target prime's save_index to within the chain code array bounds */

			if( new_length_init == _true_)
			{
				printf("\n\nStarting test_length = %u\n\n", test_length);
			}
			else /* current test length has already started */
			{
				printf("\nResuming test_length = %u\n\n", test_length);
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
					if( ( (exception_list_1_step[i] & 1) != 0 ) && (gcd( exception_list_1_step[i], 15015 ) == 1) )
					{
						j = 0;
						while( tgt_prime_list[j].prime < Luc[test_length] )
							j++;
						do
						{
							if(exception_list_1_step[i] <= tgt_prime_list[j].prime)
							{
								if(exception_list_1_step[i] == tgt_prime_list[j].prime)
								{
									prime_exception_count++;
									if( tgt_prime_list[j].prime == Luc[test_length] )
										printf("*** Lucas prime exception: Luc[%u] = %lu\n", test_length, Luc[test_length]);
									else
										printf("*** One-step prime exception > Luc[%u]: %lu\n", test_length, exception_list_1_step[i]);
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
			generate_and_process_candidate_list();


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
					if( truncating_for_B1 ) /* we are done; set exit flag */
						gen_exit_flag = _true_;
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

			unique_chain_count = 0;
			for(i = 0; i < (int32_t)(*tgt_p_count); i++)
			{
#if 0
				if( i < 10 ) /* print results for up to 10 primes at start of target list */
				{
					printf("target prime[%u] = %lu, total chain_count = %u, max # of doubles = %u, # of chains with max dbles = %u\n",
							i, tgt_prime_list[i].prime, chain_count[i], chain_max_dbl_count[i], chain_count_max_dbls[i]);
					printf("chain code for p: %016lX\n", chain_code_list[tgt_prime_list[i].save_index - *chain_code_list_start_index]);
				}
#endif

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
#if 0
					printf("*** unique chain for target prime[%u] = %lu, total chain_count = %u, max # of doubles = %u, # of chains with max dbles = %u\n",
							i, tgt_prime_list[i].prime, chain_count[i], chain_max_dbl_count[i], chain_count_max_dbls[i]);
					printf("*** chain code for p: %016lX\n", chain_code_list[tgt_prime_list[i].save_index - *chain_code_list_start_index]);
#endif
				}

				if( tgt_prime_code_length[i] > *max_code_length )
				{
					*max_code_length = tgt_prime_code_length[i];
					code_index = tgt_prime_list[i].save_index - *chain_code_list_start_index;
					printf("ALERT: new max code length = %u bits for prime = %lu, code = %016lX\n",
						*max_code_length, tgt_prime_list[i].prime, chain_code_list[code_index]);
				}

				if( *code_length_problem_count > 0 )
					printf("BAD PROBLEM: maximum code length exceeds 64 bits. A program revision to handle longer codes is required.\n");
			}

			if( unique_chain_count > 0 )
				printf("Total # of unique chains, including max double chains, = %u\n", unique_chain_count);

			if( true_indx > max_indx_value ) /* current test length is complete */
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
			fwrite((int8_t *)max_code_length, sizeof(u_int8_t), 1, current_status_file);
			fclose(current_status_file);
		}
		while( true_indx <= max_indx_value );
		if( gen_exit_flag )
		{
			dum = (size_t)system("rm -f Tgt_prime_list_1.dat");
			dum = (size_t)system("rm -f Tgt_prime_list_2.dat");
			dum = (size_t)system("rm -f Pending_Lchain_code_list_1.dat");
			dum = (size_t)system("rm -f Pending_Lchain_code_list_2.dat");
			dum = (size_t)system("rm -f current_status.dat");

			for( i = 0; i < MAX_CANDIDATE_COUNT; i++ )
				printf("raw candidate count[%d] = %lu\n", i, raw_c_counts[i]);

			return EXIT_SUCCESS;
		}
	}

	return EXIT_SUCCESS;
}
