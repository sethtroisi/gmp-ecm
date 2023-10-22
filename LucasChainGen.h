/*
 * LucasChainGen.h
 *
 *  Created on: Sep 9, 2022
 *      Author: Philip McLaughlin
 */

#ifndef LUCASCHAINGEN_H_
#define LUCASCHAINGEN_H_

#define FIB_LIMIT 56				/* maximum # of Fibonacci numbers */
#define MAX_WORKING_CHAIN_LENGTH 64	/* maximum # of chain elements in the current working chain */
#define MAX_CODE_OR_PRIME_COUNT 6000000	/* maximum size for both chain code & target prime arrays */
#define MAX_CAND_LIST_COUNT 500		/* maximum total # of candidates in the recursive list */
#define MAX_CANDIDATE_COUNT 20		/* maximum number of candidates to extend any given chain */

/* sieve parameters */
#define SBL 400						/* sieve bootstrap limit: all primes less than this limit
 	 	 	 	 	 	 	 	 	 	 are used to generate the sieve prime list */
#define SPL 160000				/* sieve prime list limit: all primes p such that 17 <= p < SPL will be sieve primes */
#define MAX_SIEVE_LIMIT 25600000000	/* integer limit above which composite integers will appear after sieving */
#define NEWSIEVE_SIZE 80000		/* byte sieve size to generate all sieve primes (SPL/2) */
#define SIEVE_SPACE_SIZE 500000		/* byte sieve size; (NOTE: must be >= NEWSIEVE_SIZE */
#define SIEVE_PRIME_COUNT 14700		/* array size estimate based on pi(SPL) approx = SPL/(ln(SPL) - 1.1) */
 	 	 	 	 	 	 	 	 	/* The true number of sieve primes is 14677 from 17 to 159979 */

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
/* 4-bit chain start sequence codes (16 observed out of 22 possible for chain lengths <= 40)
 * If # of observed start sequences is ever > 16, chain encoding will need revision */
#define CHAIN_START_5_7_10 0x3 /* never appears for encoding until test_length == 40 */
#define CHAIN_START_4_5_9  0x2
#define CHAIN_START_4_6_7  0x4
#define CHAIN_START_4_6_8  0x6
#define CHAIN_START_4_6_9  0x8
#define CHAIN_START_4_7_8  0xA
#define CHAIN_START_4_7_10 0xC
#define CHAIN_START_4_7_11 0xE
#define CHAIN_START_5_7_9  0x1 /* never appears for encoding until test_length == 38 */
#define CHAIN_START_5_8_10 0x0 /* never appears for encoding until test_length == 36 */
#define CHAIN_START_5_6_9  0x5
#define CHAIN_START_5_6_10 0x7
#define CHAIN_START_5_6_11 0x9
#define CHAIN_START_5_7_12 0xB
#define CHAIN_START_5_8_11 0xD
#define CHAIN_START_5_8_13 0xF

#define _true_ (1)
#define _false_ (0)

typedef struct
{
	u_int32_t	prime;
	u_int64_t	sieve_space_start_index;
	u_int32_t	dif_table_start_index;

} sieve_params;

/* when a target prime is the n-th prime, save_index = n - 7. For example,
 *  for prime = 17, save_index = 0; for 19 save -index = 1, etc. */
typedef struct
{
	u_int64_t	prime;
	u_int32_t	save_index;

} target_prime;

typedef struct
{
	u_int64_t	value;			/* integer value of this chain element */
	u_int64_t	gcd;			/* gcd( value, parent_value ), where "parent" is the next smaller chain element */
	u_int8_t	comp_offset_1;	/* larger summand (summand_1) component index counting back from parent (parent = 0) */
	u_int8_t	comp_offset_2;	/* smaller summand (summand_2) component index counting back from parent */
	u_int8_t	dif_offset;		/* component index of (summand_1 - summand_2) counting back from parent */
								/* note: dif_offset = 0 will indicate this is a doubled element */
	u_int8_t	chain_dbl_count;	/* total # of doubled elements in the working chain up to and including this element */

} chain_element;

/* prototypes */

u_int64_t	cputime(void);
chain_element	*get_working_chain_ptr(void);
u_int8_t	*get_current_partial_length_ptr(void); /* current # of elements in the working chain */
u_int64_t	*get_chain_values_ptr(void);
chain_element	*get_candidate_list_ptr(void);
u_int16_t	*get_c_list_start_index_ptr(void); /* next available slot in the candidate list */
u_int16_t	*get_current_c_index_ptr(void);
target_prime	*get_tgt_prime_list_ptr(void);
u_int32_t	*get_tgt_p_count_ptr(void);
u_int64_t	*get_chain_code_list_ptr(void);
u_int32_t	*get_chain_code_list_start_index_ptr(void);
u_int32_t	*get_chain_count_ptr(void);
u_int8_t	*get_chain_max_dbl_count_ptr(void);
u_int16_t	*get_chain_count_max_dbls_ptr(void);
u_int64_t	*get_Fib_ptr(void);
u_int64_t	*get_Luc_ptr(void);
u_int8_t	*get_w_chain_length_ptr(void);
u_int64_t	*get_next_step_c_list_ptr(void);
double		*get_index_count_per_val_ptr(void);
u_int8_t	*get_code_length_ptr(void);
u_int8_t	*get_max_code_length_ptr(void);
u_int8_t	*get_tgt_prime_code_length_ptr(void);
u_int32_t	*get_code_length_problem_count_ptr(void);

u_int8_t	*get_dif_table_ptr(void);
u_int8_t	*get_sieve_space_ptr(void);
sieve_params	*get_sieve_primes_ptr(void);
u_int32_t	sieve_init(void);
void		standard_sieve(u_int32_t);
u_int32_t	prime_count( u_int32_t *, u_int32_t *);
u_int8_t	extract_chain_values(void);
u_int16_t	gen_candidate_list(void);
u_int8_t	gen_next_step_candidates(void);
u_int64_t	gcd(u_int64_t, u_int64_t);
void		copy_candidate_to_working_chain(void);
u_int8_t	check_candidate(void);
u_int64_t	encode_Lchain(void);
u_int8_t	generate_Lchain( u_int64_t, u_int64_t, chain_element *, u_int8_t *, u_int8_t *, u_int32_t * );
void		max_continuation( chain_element *, u_int8_t *, u_int8_t );

void		generate_and_process_candidate_list(void);

#endif /* LUCASCHAINGEN_H_ */
