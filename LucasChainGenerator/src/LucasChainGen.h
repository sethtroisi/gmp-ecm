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
#define MAX_CANDIDATE_COUNT 24		/* maximum number of candidates to extend any given chain */
#define MAX_THREADS 16
#define TOTAL_WORK_COUNT 50			/* number of work assignments */
#define DEFAULT_THREAD_COUNT 4

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

/* 3-bit chain start sequences. Covers all possible Lucas chains */
#define CHAIN_START_5_8_13 0x7
#define CHAIN_START_5_8_11 0x6
#define CHAIN_START_5_8_10 0x5
#define CHAIN_START_5_7    0x4
#define CHAIN_START_5_6    0x3
#define CHAIN_START_4_7    0x2
#define CHAIN_START_4_5    0x1
#define CHAIN_START_4_6    0x0 /* precludes a completely zero code */

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
	u_int8_t	comp_offset_1;	/* larger summand (summand_1) component index counting back from parent (parent = 0) */
	u_int8_t	comp_offset_2;	/* smaller summand (summand_2) component index counting back from parent */
	u_int8_t	dif_offset;		/* component index of (summand_1 - summand_2) counting back from parent */
								/* note: dif_offset = 0 will indicate this is a doubled element */
	u_int8_t	chain_dbl_count;	/* total # of doubled elements in the working chain up to and including this element */

} chain_element;

/* variables used by multiple routines, one copy for each thread */
typedef struct
{
/*	target_prime tgt_prime_list[MAX_CODE_OR_PRIME_COUNT]; */
	u_int64_t chain_code_list[MAX_CODE_OR_PRIME_COUNT];
	u_int32_t chain_count[MAX_CODE_OR_PRIME_COUNT];
	u_int8_t chain_max_dbl_count[MAX_CODE_OR_PRIME_COUNT];
	u_int16_t chain_count_max_dbls[MAX_CODE_OR_PRIME_COUNT];
	u_int8_t tgt_prime_code_length[MAX_CODE_OR_PRIME_COUNT];
	chain_element working_chain[MAX_WORKING_CHAIN_LENGTH];
	u_int64_t chain_values[MAX_WORKING_CHAIN_LENGTH];
	chain_element candidate_list[MAX_CAND_LIST_COUNT];
	u_int64_t Fib[FIB_LIMIT];
	u_int64_t Luc[FIB_LIMIT];
	chain_element raw_c_list[MAX_CANDIDATE_COUNT];
	u_int8_t check_result[MAX_CANDIDATE_COUNT];
	u_int8_t check_index;
	u_int32_t chain_code_list_start_index;
	u_int32_t tgt_p_count;
	u_int8_t w_chain_length;
	u_int8_t current_partial_length; /* current # of elements in the working chain */
	u_int16_t c_list_start_index; /* next available slot in the candidate list */
	u_int16_t current_c_index;
	double index_count_per_val;
	u_int8_t code_length;

} mem_struct;

typedef struct
{
	chain_element working_chain[15];
	u_int8_t current_partial_length; /* starting # of elements in the working chain */

} work_struct;

typedef struct
{
	u_int8_t thrd_indx;

} thread_io_struct;

/* prototypes */

u_int8_t	*get_dif_table_ptr(void);
u_int8_t	*get_sieve_space_ptr(void);
sieve_params	*get_sieve_primes_ptr(void);
u_int32_t	sieve_init(void);
void		standard_sieve(u_int32_t);
u_int32_t	prime_count( u_int32_t *, u_int32_t *);

u_int64_t	cputime(void);
void		init_working_chains(void);
void 		set_work_assignments(void);
void		init_Fib_sequence(void);
void		init_Luc_sequence(void);
void		init_thread_memory(void);
void		consolidate_results(void);
void		*recursive_work(void *);

/* subroutines requiring templates */
void copy_work_assignment_to_thread( u_int8_t );
void copy_work_assignment_to_thread_01( u_int8_t );
void copy_work_assignment_to_thread_02( u_int8_t );
void copy_work_assignment_to_thread_03( u_int8_t );
void copy_work_assignment_to_thread_04( u_int8_t );
void copy_work_assignment_to_thread_05( u_int8_t );
void copy_work_assignment_to_thread_06( u_int8_t );
void copy_work_assignment_to_thread_07( u_int8_t );
void copy_work_assignment_to_thread_08( u_int8_t );
void copy_work_assignment_to_thread_09( u_int8_t );
void copy_work_assignment_to_thread_10( u_int8_t );
void copy_work_assignment_to_thread_11( u_int8_t );
void copy_work_assignment_to_thread_12( u_int8_t );
void copy_work_assignment_to_thread_13( u_int8_t );
void copy_work_assignment_to_thread_14( u_int8_t );
void copy_work_assignment_to_thread_15( u_int8_t );

u_int64_t	encode_Lchain(void);
u_int64_t	encode_Lchain_01(void);
u_int64_t	encode_Lchain_02(void);
u_int64_t	encode_Lchain_03(void);
u_int64_t	encode_Lchain_04(void);
u_int64_t	encode_Lchain_05(void);
u_int64_t	encode_Lchain_06(void);
u_int64_t	encode_Lchain_07(void);
u_int64_t	encode_Lchain_08(void);
u_int64_t	encode_Lchain_09(void);
u_int64_t	encode_Lchain_10(void);
u_int64_t	encode_Lchain_11(void);
u_int64_t	encode_Lchain_12(void);
u_int64_t	encode_Lchain_13(void);
u_int64_t	encode_Lchain_14(void);
u_int64_t	encode_Lchain_15(void);

u_int8_t	not_divisible_by_3( u_int64_t );
u_int8_t	not_divisible_by_3_01( u_int64_t );
u_int8_t	not_divisible_by_3_02( u_int64_t );
u_int8_t	not_divisible_by_3_03( u_int64_t );
u_int8_t	not_divisible_by_3_04( u_int64_t );
u_int8_t	not_divisible_by_3_05( u_int64_t );
u_int8_t	not_divisible_by_3_06( u_int64_t );
u_int8_t	not_divisible_by_3_07( u_int64_t );
u_int8_t	not_divisible_by_3_08( u_int64_t );
u_int8_t	not_divisible_by_3_09( u_int64_t );
u_int8_t	not_divisible_by_3_10( u_int64_t );
u_int8_t	not_divisible_by_3_11( u_int64_t );
u_int8_t	not_divisible_by_3_12( u_int64_t );
u_int8_t	not_divisible_by_3_13( u_int64_t );
u_int8_t	not_divisible_by_3_14( u_int64_t );
u_int8_t	not_divisible_by_3_15( u_int64_t );

u_int8_t	not_divisible_by_5( u_int64_t );
u_int8_t	not_divisible_by_5_01( u_int64_t );
u_int8_t	not_divisible_by_5_02( u_int64_t );
u_int8_t	not_divisible_by_5_03( u_int64_t );
u_int8_t	not_divisible_by_5_04( u_int64_t );
u_int8_t	not_divisible_by_5_05( u_int64_t );
u_int8_t	not_divisible_by_5_06( u_int64_t );
u_int8_t	not_divisible_by_5_07( u_int64_t );
u_int8_t	not_divisible_by_5_08( u_int64_t );
u_int8_t	not_divisible_by_5_09( u_int64_t );
u_int8_t	not_divisible_by_5_10( u_int64_t );
u_int8_t	not_divisible_by_5_11( u_int64_t );
u_int8_t	not_divisible_by_5_12( u_int64_t );
u_int8_t	not_divisible_by_5_13( u_int64_t );
u_int8_t	not_divisible_by_5_14( u_int64_t );
u_int8_t	not_divisible_by_5_15( u_int64_t );

u_int8_t	extract_chain_values(void);
u_int8_t	extract_chain_values_01(void);
u_int8_t	extract_chain_values_02(void);
u_int8_t	extract_chain_values_03(void);
u_int8_t	extract_chain_values_04(void);
u_int8_t	extract_chain_values_05(void);
u_int8_t	extract_chain_values_06(void);
u_int8_t	extract_chain_values_07(void);
u_int8_t	extract_chain_values_08(void);
u_int8_t	extract_chain_values_09(void);
u_int8_t	extract_chain_values_10(void);
u_int8_t	extract_chain_values_11(void);
u_int8_t	extract_chain_values_12(void);
u_int8_t	extract_chain_values_13(void);
u_int8_t	extract_chain_values_14(void);
u_int8_t	extract_chain_values_15(void);

void		copy_candidate_to_working_chain(void);
void		copy_candidate_to_working_chain_01(void);
void		copy_candidate_to_working_chain_02(void);
void		copy_candidate_to_working_chain_03(void);
void		copy_candidate_to_working_chain_04(void);
void		copy_candidate_to_working_chain_05(void);
void		copy_candidate_to_working_chain_06(void);
void		copy_candidate_to_working_chain_07(void);
void		copy_candidate_to_working_chain_08(void);
void		copy_candidate_to_working_chain_09(void);
void		copy_candidate_to_working_chain_10(void);
void		copy_candidate_to_working_chain_11(void);
void		copy_candidate_to_working_chain_12(void);
void		copy_candidate_to_working_chain_13(void);
void		copy_candidate_to_working_chain_14(void);
void		copy_candidate_to_working_chain_15(void);

u_int16_t	gen_candidate_list(void);
u_int16_t	gen_candidate_list_01(void);
u_int16_t	gen_candidate_list_02(void);
u_int16_t	gen_candidate_list_03(void);
u_int16_t	gen_candidate_list_04(void);
u_int16_t	gen_candidate_list_05(void);
u_int16_t	gen_candidate_list_06(void);
u_int16_t	gen_candidate_list_07(void);
u_int16_t	gen_candidate_list_08(void);
u_int16_t	gen_candidate_list_09(void);
u_int16_t	gen_candidate_list_10(void);
u_int16_t	gen_candidate_list_11(void);
u_int16_t	gen_candidate_list_12(void);
u_int16_t	gen_candidate_list_13(void);
u_int16_t	gen_candidate_list_14(void);
u_int16_t	gen_candidate_list_15(void);

void		check_candidate(void);
void		check_candidate_01(void);
void		check_candidate_02(void);
void		check_candidate_03(void);
void		check_candidate_04(void);
void		check_candidate_05(void);
void		check_candidate_06(void);
void		check_candidate_07(void);
void		check_candidate_08(void);
void		check_candidate_09(void);
void		check_candidate_10(void);
void		check_candidate_11(void);
void		check_candidate_12(void);
void		check_candidate_13(void);
void		check_candidate_14(void);
void		check_candidate_15(void);

/* u_int8_t	generate_Lchain( u_int64_t, u_int64_t, chain_element *, u_int8_t *, u_int8_t *, u_int32_t * ); */
/* void		max_continuation( chain_element *, u_int8_t *, u_int8_t ); */

void		generate_and_process_candidate_list(void);
void		generate_and_process_candidate_list_01(void);
void		generate_and_process_candidate_list_02(void);
void		generate_and_process_candidate_list_03(void);
void		generate_and_process_candidate_list_04(void);
void		generate_and_process_candidate_list_05(void);
void		generate_and_process_candidate_list_06(void);
void		generate_and_process_candidate_list_07(void);
void		generate_and_process_candidate_list_08(void);
void		generate_and_process_candidate_list_09(void);
void		generate_and_process_candidate_list_10(void);
void		generate_and_process_candidate_list_11(void);
void		generate_and_process_candidate_list_12(void);
void		generate_and_process_candidate_list_13(void);
void		generate_and_process_candidate_list_14(void);
void		generate_and_process_candidate_list_15(void);

#endif /* LUCASCHAINGEN_H_ */
