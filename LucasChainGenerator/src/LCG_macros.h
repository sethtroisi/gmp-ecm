
#define CHECK_CAN_TEMPLATE_1(thrd) \
{\
	static chain_element *working_chain;\
	static chain_element *raw_c_list;\
/*	static target_prime *tgt_prime_list; */\
	static u_int8_t *check_result, *check_index;\
	static u_int64_t *chain_code_list;\
	static u_int32_t *chain_code_list_start_index;\
	static u_int32_t *chain_count, *tgt_p_count;\
	static u_int64_t *Fib;\
	static u_int64_t *chain_values;\
	static u_int16_t *chain_count_max_dbls;\
	static u_int8_t *chain_max_dbl_count, *current_partial_length;\
	static u_int8_t *code_length;\
	static u_int8_t *tgt_prime_code_length;\
	static u_int8_t *w_chain_length, thrd_indx = (thrd), init = 0;\
	static double *index_count_per_val;\
	double approx_start_index, del;\
	u_int64_t temp_var;\
	u_int64_t gcd_c_p = 1;\
	u_int64_t c, max_val, p_val, compare_val;\
	u_int32_t k, k_old, code_index;\
	u_int64_t beta, upper_limit_1_step = 0, a1, a2 = 0;\
	u_int8_t kk, steps_to_go, c_length, element_count, doubles_count;\
	u_int8_t compare_val_in_chain, next_c_count, max_c_flag;\
	u_int64_t next_step_c_list[MAX_CANDIDATE_COUNT];\
	u_int64_t cand, dif;\
	u_int8_t i, j, c_flag, ii;\
\
	if( init == 0 )\
	{\
		working_chain = thread_mem[thrd_indx].working_chain;\
		chain_values = thread_mem[thrd_indx].chain_values;\
		raw_c_list = thread_mem[thrd_indx].raw_c_list;\
		check_index = &thread_mem[thrd_indx].check_index;\
		check_result = thread_mem[thrd_indx].check_result;\
/*		tgt_prime_list = thread_mem[thrd_indx].tgt_prime_list; */\
		chain_code_list = thread_mem[thrd_indx].chain_code_list;\
		chain_code_list_start_index = &thread_mem[thrd_indx].chain_code_list_start_index;\
		chain_count = thread_mem[thrd_indx].chain_count;\
		chain_count_max_dbls = thread_mem[thrd_indx].chain_count_max_dbls;\
		chain_max_dbl_count = thread_mem[thrd_indx].chain_max_dbl_count;\
		tgt_p_count = &thread_mem[thrd_indx].tgt_p_count;\
		current_partial_length = &thread_mem[thrd_indx].current_partial_length;\
		code_length = &thread_mem[thrd_indx].code_length;\
		tgt_prime_code_length = thread_mem[thrd_indx].tgt_prime_code_length;\
		w_chain_length = &thread_mem[thrd_indx].w_chain_length;\
		Fib = thread_mem[thrd_indx].Fib;\
		index_count_per_val = &thread_mem[thrd_indx].index_count_per_val;\
		init = 1;\
	}\
\
	c = raw_c_list[*check_index].value;\
	steps_to_go = (u_int8_t)(*w_chain_length - *current_partial_length - 1);\
	max_val = Fib[steps_to_go + 1]*c + Fib[steps_to_go]*chain_values[0];\
	/* note that max_val == c when steps_to_go == 0 */\
\
	if( (max_val < tgt_prime_list[0].prime) || (c > tgt_prime_list[*tgt_p_count - 1].prime) )\
	{\
		check_result[*check_index] = 0;\
		return;\
	}\
\
	if( steps_to_go > 0 ) /* c is not the last (largest) element */\
	{\
		/* check if candidate gives a chain of the form [...,x,...,2x,3x],\
		 * with all of the elements between x and 2x less than 3x/2 */\
		if(working_chain[*current_partial_length].dif_offset == 0) /* doubled element = 2x */\
		{\
			p_val = chain_values[0];\
			if( (c == 3*p_val/2) /* candidate = 3x */\
				&& (2*chain_values[1] < 3*p_val/2 )) /* largest element below 2x is < 3x/2 */\
			{\
				/* check if 2x - e is an element in the chain for any element e such that x < e < 2x  */\
				c_length = *current_partial_length;\
				element_count = 2; /* find the index where chain_values[index] == p_val/2 */\
				while( chain_values[element_count] > p_val/2 )\
					element_count++;\
\
				/* sanity check */\
				if( 2*chain_values[element_count] != p_val )\
					printf("ERROR: no chain element == p_val/2!\n");\
\
				k = 1;\
				compare_val_in_chain = 0;\
				do\
				{\
					compare_val = p_val - chain_values[k];\
					kk = c_length;\
					do\
					{\
						if(compare_val == chain_values[kk] )\
						{\
							compare_val_in_chain = 1; /* a candidate to extend the chain exists which is NOT a multiple of p_val/2 */\
							break;\
						}\
						kk--;\
					}\
					while( kk > element_count );\
					if(compare_val_in_chain != 0) /* the 2x3x theorem does not apply */\
						break;\
					k++;\
				}\
				while( chain_values[k] > p_val/2 );\
\
				if(compare_val_in_chain == 0) /* the 2x3x theorem applies; we can ignore the candidate */\
				{\
					check_result[*check_index] = 0;\
					return;\
				}\
			}\
		} /* end check if 2x3x chain theorem applies */\
\
		/* The following only checks for a string of even elements */\
		if( !(c & 0x1) && !(chain_values[0] & 0x1 ) )\
			gcd_c_p = 2;\
		else\
			gcd_c_p = 1;\
\
		/* check for chains which never lead to primes */\
		if( gcd_c_p > 1 )\
		{\
			/* check if c, chain_values[0], and chain_values[1] are all even */\
			kk = 1;\
			do\
			{\
				if( chain_values[kk] & 0x1 )\
					break;\
				kk++;\
			}\
			while(1);\
\
			if(kk >= 2)\
			{\
				 /* c, chain_values[0],..., chain_values[kk - 1]\
				    share a nontrivial common divisor */\
				compare_val = chain_values[kk]; /* gcd(compare_val, gcd_c_p) = 1 */\
\
				/* If 2*compare_val > c, then it is a valid candidate to extend the chain including c.\
					If 2*compare_val == c, then c + compare_val is a valid candidate, since\
					c - compare_val = compare_val is in the chain */\
				if(c <= 2*compare_val) /* "3-strikes" theorem does not apply */\
					goto bypass;\
\
				/* for each smaller element divisible by gcd_c_p, check that either\
				   (element + compare_val) < c OR element > 2*compare_val (or both)\
				   so that there are NO candidates > c involving compare_val, or\
				   any chain element < compare_val. It is then an easy inductive proof\
				   to show that all chain elements > c are divisible by gcd_c_p\
				   in any chain containing [chain_values[...], c], and hence never prime */\
				for(k=0;k < kk;k++)\
				{\
					if( ((chain_values[k] + compare_val) > c ) && (chain_values[k] <= 2*compare_val) ) /* might be a candidate */\
					{\
						 /* if chain_values[k] > 2*compare_val, then (chain_values[k] - compare_val) > compare_val.\
						    Since the only elements in the chain greater than compare_val are divisible by gcd_c_p,\
						    (chain_values[k] - compare_val) is never a chain member */\
						goto bypass; /* theorem does not apply */\
						/* note that chain_values[k] + compare_val != c, since equality would force\
						   compare_val to be divisible by gcd_c_p, a contradiction */\
					}\
				}\
\
				/* divisibility theorem applies; toss candidate out */\
				check_result[*check_index] = 0;\
				return;\
			}\
		} /* end if gcd_c_p > 1 */\
\
bypass:\
		/* check if the one-step upper limit lemma applies */\
\
		/* get a list of all candidates for the next chain including c\
		   which are greater than or equal to 2*chain_values[0] */\
		next_step_c_list[0] = 2*chain_values[0];\
		next_c_count = 1;\
\
		j = 0;\
		do\
		{\
			cand = c + chain_values[j];\
			if( cand <= next_step_c_list[0] )\
				break;\
\
			dif = c - chain_values[j];\
			if( dif > chain_values[0] )\
				break;\
\
			/* check if dif is an element in the chain */\
			k = 0;\
			do\
			{\
				if( dif >= chain_values[k] )\
				{\
					if( dif == chain_values[k] )\
					{\
						/* make sure that cand is not already in the candidate list */\
						c_flag = 1;\
						ii = 0;\
						while( ii < next_c_count )\
						{\
							if( cand == next_step_c_list[ii] )\
							{\
								c_flag = 0;\
								break;\
							}\
							ii++;\
						}\
						if( c_flag )\
						{\
							/* add cand to the list */\
							next_step_c_list[next_c_count] = cand;\
							next_c_count++;\
						}\
					}\
					break;\
				}\
				k++;\
			}\
			while( 1 );\
			j++;\
		}\
		while( j <= (*current_partial_length) );\
\
		/* bubble sort the list to get largest to smallest */\
		for( i = (next_c_count - 1); i > 0; i--)\
		{\
			for(j=0;j<i;j++)\
			{\
				if(next_step_c_list[j] < next_step_c_list[j+1])\
				{\
					dif = next_step_c_list[j];\
					next_step_c_list[j] = next_step_c_list[j+1];\
					next_step_c_list[j+1] = dif;\
				}\
			}\
		}\
\
		/* find beta (the largest candidate strictly less than the maximal element) */\
		if( next_step_c_list[0] == ( c + chain_values[0] ) )\
		{\
			max_c_flag = 1; /* maximum chain continuation exists */\
			beta = next_step_c_list[1];\
		}\
		else /* next_step_c_list[0] < (c + chain_values[0]) */\
		{\
			max_c_flag = 0; /* maximal chain continuation does not exist */\
			beta = next_step_c_list[0];\
		}\
\
		/* determine the upper limit value that the current working chain can reach\
		 *  in steps_to_go by choosing, at one of the steps, the largest candidate\
		 *  less than (a(i) + a(i-1)), and the maximal candidate otherwise */\
		a1 = Fib[steps_to_go]*beta + Fib[steps_to_go - 1]*c;\
		if( (steps_to_go <= 2) || (max_c_flag == 0) )\
			upper_limit_1_step = a1;\
		 /* steps_to_go >= 3 && max_c_flag == 1 from here*/\
		else if(steps_to_go == 3)\
			a2 = 2*c + 3*chain_values[0];\
		else if( (2*c <= 3*chain_values[0]) ) /* steps_to_go >= 4 */\
			a2 = Fib[steps_to_go]*c + Fib[steps_to_go + 1]*chain_values[0];\
		else\
			a2 = (Fib[steps_to_go] + 2*Fib[steps_to_go - 3])*c + (Fib[steps_to_go - 1] + 2*Fib[steps_to_go - 4])*chain_values[0];\
\
		if( (steps_to_go >= 3) && (max_c_flag == 1) )\
		{\
			if(a1 >= a2)\
				upper_limit_1_step = a1;\
			else\
				upper_limit_1_step = a2;\
		}\
\
		if(upper_limit_1_step >= tgt_prime_list[0].prime)\
		{\
			check_result[*check_index] = 1;\
			return;\
		}\
		else\
		{\
			if( max_val > tgt_prime_list[*tgt_p_count - 1].prime )\
			{\
				check_result[*check_index] = 0; /* an upcoming interval will catch this one */\
				return;\
			}\


#define CHECK_CAN_TEMPLATE_4 \
				{\
					doubles_count = raw_c_list[*check_index].chain_dbl_count;\
\
					/* the code blitzes through the shorter lengths so fast\
					 * that we don't need to adjust the start point for k */\
					if( *w_chain_length > 15)\
					{\
						/* find approximate index to start comparison loop\
						 	note that max_val >= tgt_prime_list[0].prime here, so approx_start_index >= 0 */\
						approx_start_index = (double)( max_val - tgt_prime_list[0].prime )*( *index_count_per_val );\
\
						k = (u_int32_t)approx_start_index;\
						/* Note: since  tgt_prime_list[0].prime <= max_val <= Fib[test_length + 2],\
						 * and *index_count_per_val = (*tgt_p_count - 1) / (Fib[test_length + 2] - tgt_prime_list[0].prime)\
						 * in double precision floating point, then  0 <= k <= (*tgt_p_count - 1) */\
\
						/* if max_val < tgt_prime_list[k].prime,\
						 *  reduce k until max_val >= tgt_prime_list[k].prime */\
						while( max_val < tgt_prime_list[k].prime ) /* if true, k must be > 0 */\
						{\
							del = ( (double)( tgt_prime_list[k].prime - max_val ) )*( *index_count_per_val ) + 1.0;\
							approx_start_index -= del;\
							if( approx_start_index < 0 ) approx_start_index = 0;\
							k = (u_int32_t)approx_start_index;\
						}\
\
						del = ( (double)( max_val - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;\
						k_old = k + (u_int32_t)del;\
\
						while( (k_old < *tgt_p_count) && (max_val > tgt_prime_list[k_old].prime) )\
						{\
							k = k_old;\
							del = ( (double)( max_val - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;\
							k_old = k + (u_int32_t)del;\
						}\
					}\
					else\
						k = 0;\
\
					while(k < *tgt_p_count)\
					{\
						if(max_val <= tgt_prime_list[k].prime)\
						{\
							if(max_val == tgt_prime_list[k].prime)\
							{\
								chain_count[k]++;\
								if( doubles_count >= chain_max_dbl_count[k] )\
								{\
									if(doubles_count > chain_max_dbl_count[k])\
									{\
										/* best chain so far */\
										chain_max_dbl_count[k] = doubles_count;\
										chain_count_max_dbls[k] = 1;\
										/* Encode & save this chain. Code will be overwritten\
											if/when a better chain with more doubled elements is found */\
										code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;

#define CHECK_CAN_TEMPLATE_5 \
\
/*										if( (tgt_prime_code_length[k]) > 64 && (*code_length <= 64) )\
										{\
											printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",\
													k, tgt_prime_list[k].prime);\
										} */\
\
										tgt_prime_code_length[k] = *code_length;\
/*										if(*code_length > 64)\
										{\
											printf("PROBLEM: chain code length > 64 bits for tgt_prime_list[%u] = %lu\n",\
													k, tgt_prime_list[k].prime);\
										} */\
									}\
									else\
									{\
										chain_count_max_dbls[k]++;\
										if( tgt_prime_code_length[k] > 64 )\
										{\
											/* given this choice, see if a code exists with <= 64 bits */

#define CHECK_CAN_TEMPLATE_6 \
											if( *code_length <= 64 )\
											{\
												code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;\
												chain_code_list[code_index] = temp_var;\
												tgt_prime_code_length[k] = *code_length;\
/*												printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",\
														k, tgt_prime_list[k].prime); */\
											}\
										}\
									}\
								}\
							}\
							break;\
						}\
						k++;\
					}\
				}\
			}\
			check_result[*check_index] = 0; /* this candidate has been resolved */\
			return;\
		}\
	} /* end if steps_to_go > 0 */\
	else /* steps_to_go == 0 */\
	{\
		if( ( c % 7 ) && ( c % 11 ) && ( c % 13 ) )\
		{\
			doubles_count = raw_c_list[*check_index].chain_dbl_count;\
\
			/* the code blitzes through the shorter lengths so fast\
			 * that we don't need to adjust the start point for k */\
			if( *w_chain_length > 15)\
			{\
				/* find approximate index to start comparison loop\
				 	note that c >= tgt_prime_list[0].prime here, so approx_start_index >= 0 */\
				approx_start_index = (double)( c - tgt_prime_list[0].prime )*( *index_count_per_val );\
\
				k = (u_int32_t)approx_start_index;\
\
				/* if c < tgt_prime_list[k].prime,\
				 *  reduce k until c >= tgt_prime_list[k].prime */\
				while( c < tgt_prime_list[k].prime  ) /* if true, k must be > 0 */\
				{\
					del = ( (double)( tgt_prime_list[k].prime - c ) )*( *index_count_per_val ) + 1.0;\
					approx_start_index -= del;\
					if( approx_start_index < 0 ) approx_start_index = 0;\
					k = (u_int32_t)approx_start_index;\
				}\
\
				del = ( (double)( c - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;\
				k_old = k + (u_int32_t)del;\
\
				while( (k_old < *tgt_p_count) && (c > tgt_prime_list[k_old].prime) )\
				{\
					k = k_old;\
					del = ( (double)( c - tgt_prime_list[k].prime ) )*( *index_count_per_val )*0.5 + 1.0;\
					k_old = k + (u_int32_t)del;\
				}\
			}\
			else\
				k = 0;\
\
			while(k < *tgt_p_count)\
			{\
				if(c <= tgt_prime_list[k].prime)\
				{\
					if(c == tgt_prime_list[k].prime)\
					{\
						chain_count[k]++;\
						if( doubles_count >= chain_max_dbl_count[k] )\
						{\
							if(doubles_count > chain_max_dbl_count[k])\
							{\
								/* best chain so far */\
								chain_max_dbl_count[k] = doubles_count;\
								chain_count_max_dbls[k] = 1;\
\
								/* Encode & save this chain. Code will be overwritten\\
									if/when a better chain with more doubled elements is found */\
								code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;

#define CHECK_CAN_TEMPLATE_8 \
/*								if( (tgt_prime_code_length[k]) > 64 && (*code_length <= 64) )\
								{\
									printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",\
											k, tgt_prime_list[k].prime);\
								} */\
\
								tgt_prime_code_length[k] = *code_length;\
\
/*								if(*code_length > 64)\
								{\
									printf("PROBLEM: chain code length > 64 bits for tgt_prime_list[%u] = %lu\n",\
											k, tgt_prime_list[k].prime);\
								} */\
							}\
							else\
							{\
								chain_count_max_dbls[k]++;\
								if( tgt_prime_code_length[k] > 64 )\
								{\
									/* given this choice, see if a code exists with <= 64 bits */

#define CHECK_CAN_TEMPLATE_9 \
									if( *code_length <= 64 )\
									{\
										code_index = tgt_prime_list[k].save_index - *chain_code_list_start_index;\
										chain_code_list[code_index] = temp_var;\
										tgt_prime_code_length[k] = *code_length;\
/*										printf("Resolved problem: new chain code length <= 64 bits for tgt_prime_list[%u] = %lu\n",\
												k, tgt_prime_list[k].prime); */\
									}\
								}\
							}\
						}\
					}\
					break;\
				}\
				k++;\
			}\
		}\
		check_result[*check_index] = 0; /* this candidate has been resolved */\
		return;\
	} /* end steps_to_go == 0 */\
}





#define GEN_AND_PROCESS_C_LIST_1(thrd) \
{\
	static u_int16_t *c_list_start_index, *current_c_index;\
	static u_int16_t c_count[MAX_WORKING_CHAIN_LENGTH];\
	static u_int16_t c_indx[MAX_WORKING_CHAIN_LENGTH];\
	static u_int8_t *current_partial_length;\
	static u_int8_t thrd_indx = (thrd), init = 0, r_level = 0;\
\
	if( init == 0 )\
	{\
		/* initialize common variable & variable array pointers */\
		current_partial_length = &thread_mem[thrd_indx].current_partial_length;\
		c_list_start_index = &thread_mem[thrd_indx].c_list_start_index;\
		current_c_index = &thread_mem[thrd_indx].current_c_index;\
		init = 1;\
	}\
\
	/* find all candidates to be the next element of the working chain */\
	/* final-step-in-the-chain candidates are processed in the following call */

#define GEN_AND_PROCESS_C_LIST_2 \
\
	if(c_count[r_level] == 0) /* no candidates to process. Final step always returns zero */\
		return;\
\
	c_indx[r_level] = 0;\
	do\
	{\
		*current_c_index = *c_list_start_index + c_indx[r_level];\
		/* working_chain[*current_partial_length + 1] = candidate_list[*current_c_index] */\
		/* copy candidate to working chain (increments *current_partial_length) */

#define GEN_AND_PROCESS_C_LIST_3 \
		*c_list_start_index += c_count[r_level];\
		r_level++;

#define GEN_AND_PROCESS_C_LIST_4 \
		r_level--;\
		*c_list_start_index -= c_count[r_level];\
		(*current_partial_length)--;\
\
		c_indx[r_level]++;\
	}\
	while(c_indx[r_level] < c_count[r_level]);\
}


/* u_int8_t extract_chain_values(void) */
#define EXTRACT_CHAIN_VALUES(thrd) {\
	static chain_element *working_chain;\
	static u_int64_t *chain_values;\
	static u_int8_t *current_partial_length;\
	static u_int8_t thrd_indx = (thrd), init = 0;\
	u_int8_t i, indx, double_count;\
\
	if( init == 0 )\
	{\
		/* initialize array pointers */\
		working_chain = thread_mem[thrd_indx].working_chain;\
		chain_values = thread_mem[thrd_indx].chain_values;\
		current_partial_length = &thread_mem[thrd_indx].current_partial_length;\
		init = 1;\
	}\
\
	indx = *current_partial_length;\
	i = 0;\
	double_count = working_chain[ indx ].chain_dbl_count;\
	do\
	{\
		chain_values[i] = working_chain[ indx ].value;\
		i++;\
		indx--;\
	}\
	while( indx > 0 );\
\
	chain_values[i] = 1;\
	return double_count;\
}

#define COPY_C_TO_W_CHAIN(thrd) \
{\
	static chain_element *working_chain, *candidate_list;\
	static u_int8_t *current_partial_length;\
	static u_int16_t *current_c_index;\
	static u_int8_t thrd_indx = (thrd), init = 0;\
	u_int16_t c_index;\
	u_int8_t w_index;\
\
	if(init == 0)\
	{\
		working_chain = thread_mem[thrd_indx].working_chain;\
		current_partial_length = &thread_mem[thrd_indx].current_partial_length;\
		candidate_list = thread_mem[thrd_indx].candidate_list;\
		current_c_index = &thread_mem[thrd_indx].current_c_index;\
		init = 1;\
	}\
\
	w_index = (u_int8_t)(*current_partial_length + 1);\
	c_index = *current_c_index;\
\
	working_chain[w_index].value           = candidate_list[c_index].value;\
	working_chain[w_index].comp_offset_1   = candidate_list[c_index].comp_offset_1;\
	working_chain[w_index].comp_offset_2   = candidate_list[c_index].comp_offset_2;\
	working_chain[w_index].dif_offset      = candidate_list[c_index].dif_offset;\
	working_chain[w_index].chain_dbl_count = candidate_list[c_index].chain_dbl_count;\
\
	*current_partial_length = w_index;\
}

#define GEN_C_LIST_1(thrd) \
{\
	static chain_element *candidate_list;\
	static chain_element *raw_c_list;\
	static u_int64_t *chain_values;\
	static u_int16_t *c_list_start_index;\
	static u_int8_t *current_partial_length;\
	static u_int8_t *w_chain_length, thrd_indx = (thrd), init = 0;\
	static u_int8_t *check_result, *check_index;\
\
	u_int64_t dif, c;\
	u_int16_t c_index, c_count, ii;\
	u_int8_t i, j, k, c_flag, doubles_count, steps_to_go, rcl_index;\
\
	if( init == 0 )\
	{\
		/* initialize pointers */\
		candidate_list = thread_mem[thrd_indx].candidate_list;\
		raw_c_list = thread_mem[thrd_indx].raw_c_list;\
		check_result = thread_mem[thrd_indx].check_result;\
		check_index = &thread_mem[thrd_indx].check_index;\
		c_list_start_index = &thread_mem[thrd_indx].c_list_start_index;\
		chain_values = thread_mem[thrd_indx].chain_values;\
		current_partial_length = &thread_mem[thrd_indx].current_partial_length;\
		w_chain_length = &thread_mem[thrd_indx].w_chain_length;\
		for( i = 0; i < MAX_CANDIDATE_COUNT; i++)\
			check_result[i] = 0;\
		init = 1;\
	}\
\
	steps_to_go = *w_chain_length - *current_partial_length;\
	rcl_index = 0;\
\
	/* extract the chain into the chain_values array. Return value\
	 * is the # of doubled elements in the working chain */

#define GEN_C_LIST_2 \
\
	/* First, find all possible doubled candidates. It is easy to prove that\
		there is always at least one, namely 2*chain_values[1] */\
	if( steps_to_go > 1 ) /* if == 1, doubled candidates cannot be prime */\
	{\
		i = 1;\
		while( (c = 2*chain_values[i]) > chain_values[0] )\
		{\
			raw_c_list[rcl_index].value = c;\
			raw_c_list[rcl_index].comp_offset_1 = i;\
			raw_c_list[rcl_index].comp_offset_2 = i;\
			raw_c_list[rcl_index].dif_offset = 0; /* indicates a doubled element */\
			raw_c_list[rcl_index].chain_dbl_count = doubles_count + 1;\
			*check_index = rcl_index;

#define GEN_C_LIST_3 \
			rcl_index++;\
			i++;\
		}\
	}\
\
	i = 0;\
	j = 1;\
	do\
	{\
		do\
		{\
			c = chain_values[i] + chain_values[j];\
			if( c <= chain_values[0] )\
				break;\
\
			dif = chain_values[i] - chain_values[j];\
			if( dif > chain_values[i + 1] )\
				break;\
\
			/* check if dif is an element in the chain */\
			k = i + 1;\
			do\
			{\
				if( dif >= chain_values[k] )\
				{\
					if( dif == chain_values[k] )\
					{\
						/* make sure that c is not already in the candidate list */\
						c_flag = 1;\
						ii = 0;\
						while( ii < rcl_index )\
						{\
							if( c == raw_c_list[ii].value )\
							{\
								c_flag = 0;\
								break;\
							}\
							ii++;\
						} /* don't list if steps_to_go == 1 && c is divisible by 2, 3, or 5 */

#define GEN_C_LIST_3A \
						{\
							/* populate the next available slot. Done this way in case\
							  the full chain, or a maximally extended partial chain,\
							  requires encoding, which is done in the check_candidate() routine */\
							raw_c_list[rcl_index].value = c;\
							raw_c_list[rcl_index].comp_offset_1 = i;\
							raw_c_list[rcl_index].comp_offset_2 = j;\
							raw_c_list[rcl_index].dif_offset = k;\
							raw_c_list[rcl_index].chain_dbl_count = doubles_count;\
							*check_index = rcl_index;

#define GEN_C_LIST_4 \
							rcl_index++;\
							if( rcl_index == MAX_CANDIDATE_COUNT )\
								printf("\nALERT: rcl_index reached maximum of %d - increase maximum\n\n", (int)MAX_CANDIDATE_COUNT );\
						}\
					}\
					break;\
				}\
				k++;\
			}\
			while( 1 );\
			j++;\
		}\
		while( j <= *current_partial_length );\
		i++;\
		j = i+1;\
	}\
	while( (chain_values[i] + chain_values[j]) > chain_values[0] );\
\
	c_index = *c_list_start_index;\
	c_count = 0;\
	for(i = 0; i < rcl_index; i++)\
	{\
		if( check_result[i] )\
		{\
			check_result[i] = 0;\
			candidate_list[c_index].value = raw_c_list[i].value;\
			candidate_list[c_index].comp_offset_1 = raw_c_list[i].comp_offset_1;\
			candidate_list[c_index].comp_offset_2 = raw_c_list[i].comp_offset_2;\
			candidate_list[c_index].dif_offset = raw_c_list[i].dif_offset;\
			candidate_list[c_index].chain_dbl_count = raw_c_list[i].chain_dbl_count;\
			c_index++;\
			c_count++;\
			if(c_index == MAX_CAND_LIST_COUNT)\
			{\
				printf("Candidate array is full. Increase MAX_CAND_LIST_COUNT\n");\
				return c_count;\
			}\
		}\
	}\
	return c_count; /* total # of valid candidates to extend the current partial working chain */\
}

#define NOT_DIVISIBLE_3 \
{\
	u_int64_t q, r;\
\
	q = arg*0x55555555;\
	r = q & 0xFFFFFFFF;\
	q >>= 32;\
	if( 0xFFFFFFFF - (q + r) )\
		return 1;\
	else\
		return 0;\
}

#define NOT_DIVISIBLE_5 \
{\
	u_int64_t q, r;\
\
	q = arg*0x33333333;\
	r = q & 0xFFFFFFFF;\
	q >>= 32;\
	if( 0xFFFFFFFF - (q + r) )\
		return 1;\
	else\
		return 0;\
}


#define COPY_WORK_TO_THREAD(thrd) \
{\
	chain_element *work, *thrd_wrk;\
	u_int8_t cpl, i, thrd_indx = (thrd);\
\
	work = work_assignment[wrk_indx].working_chain;\
	cpl = work_assignment[wrk_indx].current_partial_length;\
	thrd_wrk = thread_mem[thrd_indx].working_chain;\
\
	thread_mem[thrd_indx].current_partial_length = cpl;\
	thread_mem[thrd_indx].c_list_start_index = 0;\
\
	for( i = 3; i <= cpl; i++)\
	{\
		thrd_wrk[i].value = work[i].value;\
		thrd_wrk[i].comp_offset_1 = work[i].comp_offset_1;\
		thrd_wrk[i].comp_offset_2 = work[i].comp_offset_2;\
		thrd_wrk[i].dif_offset = work[i].dif_offset;\
		thrd_wrk[i].chain_dbl_count = work[i].chain_dbl_count;\
	}\
}

/* u_int64_t encode_Lchain(void) */
#define ENCODE_LCHAIN_TEMPLATE(thrd) \
{\
	static chain_element *working_chain, *raw_c_list;\
	static u_int8_t *current_partial_length;\
	static u_int8_t *check_index;\
	static u_int8_t *code_length;\
	static u_int8_t thrd_indx = (thrd), init = 0;\
\
	u_int64_t chain_code;\
	u_int64_t val3, val4, val5;\
	u_int8_t i, j, len, index1, index2;\
	u_int8_t step_type[MAX_WORKING_CHAIN_LENGTH];\
	u_int8_t step_index2[MAX_WORKING_CHAIN_LENGTH];\
	u_int8_t type_count, type_index, increment, lower_limit;\
\
	if( init == 0 )\
	{\
		working_chain = thread_mem[thrd_indx].working_chain;\
		current_partial_length = &thread_mem[thrd_indx].current_partial_length;\
		raw_c_list = thread_mem[thrd_indx].raw_c_list;\
		check_index = &thread_mem[thrd_indx].check_index;\
		code_length = &thread_mem[thrd_indx].code_length;\
		init = 1;\
	}\
\
	/* We assume that any application using Lucas chains for primes will\
	 * handle the small primes 2, 3, 5, and 7 without using/requiring a chain code.\
	 * Codes for primes 11 (CHAIN_START_4_7_11) and 13 (CHAIN_START_5_8_13)\
	 * will be written directly to the beginning of the chain code save file */\
\
	chain_code = 0;\
	*code_length = 0;\
	len = *current_partial_length + 1;\
\
	/* check for chain start 1 2 3 5 8 */\
	if( len >= 5 )\
	{\
		if( (working_chain[3].value == 5) && (working_chain[4].value == 8) )\
			lower_limit = 5;\
		else\
			lower_limit = 4;\
	}\
	else if(len == 4)\
	{\
		if( (working_chain[3].value == 5) && (raw_c_list[*check_index].value == 8) )\
			lower_limit = 5;\
		else\
			lower_limit = 4;\
	}\
	else if( len == 3 )\
	{\
		if( raw_c_list[*check_index].value == 5 )\
			lower_limit = 5;\
		else\
			lower_limit = 4;\
	}\
	else\
		printf("ERROR: encode_lchain called for chain length (%u) less than 3!\n", len);\
\
	if( len > lower_limit )\
	{\
		/* set indexes for top element (raw_c_list[*check_index]) */\
		index1 = raw_c_list[*check_index].comp_offset_1;\
		index2 = raw_c_list[*check_index].comp_offset_2;\
\
		/* ignore consecutive maximal steps at the top end (largest elements) of the chain */\
		if( (index1 == 0) && (index2 == 1) )\
		{\
			do\
			{\
				len--;\
				index1 = working_chain[len].comp_offset_1;\
				index2 = working_chain[len].comp_offset_2;\
			}\
			while( (len > lower_limit) && (index1 == 0) && (index2 == 1) );\
		}\
\
		if( len > lower_limit )\
		{\
			/* set step_type for each element with working chain index > 5 */\
			type_count = 0;\
			while( len > lower_limit )\
			{\
				if( (index1 == 0) && (index2 == 1) )\
					step_type[type_count] = 1;\
				else if( (index1 == 0) && (index2 == 2) )\
					step_type[type_count] = 2;\
				else if( (index1 == 1) && (index2 == 1) )\
					step_type[type_count] = 3;\
				else if( (index1 == 2) && (index2 == 2) )\
					step_type[type_count] = 4;\
				else if( (index1 == 0) && (index2 > 2) )\
				{\
					step_type[type_count] = 5;\
					step_index2[type_count] = index2;\
				}\
				else if( (index1 == 1) && (index2 > 1) )\
				{\
					step_type[type_count] = 6;\
					step_index2[type_count] = index2;\
				}\
				else if( (index1 == 2) && (index2 == 3) )\
					step_type[type_count] = 7;\
				else if( (index1 == 2) && (index2 == 4) )\
					step_type[type_count] = 8;\
				else\
					printf("ERROR: unimplemented step_type: %u %u\n", index1, index2);\
\
				type_count++;\
				len--;\
				index1 = working_chain[len].comp_offset_1;\
				index2 = working_chain[len].comp_offset_2;\
			}\
			step_type[type_count] = 0;\
\
			type_index = 0;\
			do\
			{\
				*code_length += 4;\
				chain_code <<= 4;\
				switch( step_type[type_index] )\
				{\
					case 1:\
						/* Note: during top-down encoding, a step of this type\
						 * will always follow only step types 2 or 3, which both cover\
						 * multiple subsequent maximal continuations (step type == 1) */\
						printf("ERROR: isolated step type = 1.\n"); /* should never happen */\
						return 0;\
\
					case 2:\
						j = 0;\
						while( step_type[type_index + j + 1] == 1 )\
							j++;\
\
						increment = j + 1;\
\
						i = 0;\
						while( j >= 12 )\
						{\
							j -= 12;\
							i++;\
						}\
\
						switch( j )\
						{\
							case 0:\
								chain_code += 0x1;\
								break;\
							case 1:\
								chain_code += 0x3;\
								break;\
							case 2:\
								chain_code += 0x5;\
								break;\
							case 3:\
								chain_code += 0x7;\
								break;\
							case 4:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x9;\
								break;\
							case 5:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x19;\
								break;\
							case 6:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x29;\
								break;\
							case 7:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x39;\
								break;\
							case 8:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0xB;\
								break;\
							case 9:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x1B;\
								break;\
							case 10:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x2B;\
								break;\
							case 11:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x3B;\
								break;\
							default: /* placeholder - should never get here */\
								printf("ERROR: reached default #2\n");\
								break;\
						}\
						switch( i ) /* step_type = 1, i*12 */\
						{\
							case 0:\
								break;\
							case 1:\
								chain_code <<= 6;\
								*code_length += 6;\
								/* chain_code += 0x0;  12 maximal elements  */\
								break;\
							case 2:\
								chain_code <<= 6;\
								*code_length += 6;\
								chain_code += 0x10; /* 24 maximal elements */\
								break;\
							case 3:\
								chain_code <<= 6;\
								*code_length += 6;\
								chain_code += 0x20; /* 36 maximal elements */\
								break;\
							/* Note: the code fragment 0x30 is used for step type = 4, since a string of 48\
							   maximal elements will never appear until chain lengths are well into the 50's */\
							default: /* placeholder - should never get here */\
								printf("ERROR: reached default #3\n");\
								break;\
						}\
						break;\
\
					case 3:\
						j = 0;\
						while( step_type[type_index + j + 1] == 1 )\
							j++;\
\
						increment = j + 1;\
\
						i = 0;\
						while( j >= 12 )\
						{\
							j -= 12;\
							i++;\
						}\
\
						switch( j )\
						{\
							case 0:\
								chain_code += 0x2;\
								break;\
							case 1:\
								chain_code += 0x4;\
								break;\
							case 2:\
								chain_code += 0x6;\
								break;\
							case 3:\
								chain_code += 0x8;\
								break;\
							case 4:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0xA;\
								break;\
							case 5:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x1A;\
								break;\
							case 6:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x2A;\
								break;\
							case 7:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x3A;\
								break;\
							case 8:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0xC;\
								break;\
							case 9:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x1C;\
								break;\
							case 10:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x2C;\
								break;\
							case 11:\
								chain_code <<= 2;\
								*code_length += 2;\
								chain_code += 0x3C;\
								break;\
							default: /* placeholder - should never get here */\
								printf("ERROR: reached default #4\n");\
								break;\
						}\
						switch( i ) /* step_type = 1, i*12 */\
						{\
						case 0:\
							break;\
						case 1:\
							chain_code <<= 6;\
							*code_length += 6;\
							/* chain_code += 0x0;  12 maximal elements */\
							break;\
						case 2:\
							chain_code <<= 6;\
							*code_length += 6;\
							chain_code += 0x10; /* 24 maximal elements */\
							break;\
						case 3:\
							chain_code <<= 6;\
							*code_length += 6;\
							chain_code += 0x20; /* 36 maximal elements */\
							break;\
						default: /* placeholder - should never get here */\
							printf("ERROR: reached default #5\n");\
							break;\
						/* Note: the code fragment 0x30 is used for step type = 4, since a string of 48\
						   maximal elements will never appear until chain lengths are well into the 50's */\
						}\
						break;\
\
					case 4:\
						increment = 1;\
						chain_code <<= 2;\
						*code_length += 2;\
						chain_code += 0x30;\
						break;\
\
					case 5:\
						i = step_index2[type_index] - 3;\
						chain_code <<= 2;\
						*code_length += 2;\
						increment = 1;\
						switch(i)\
						{\
							case 0:\
							case 1:\
							case 2:\
							case 3:\
								chain_code += (u_int64_t)(i*16 + 0xE);\
								break;\
							case 4:\
								chain_code += 0x0D;\
								break;\
							case 5:\
								chain_code += 0x1D;\
								break;\
							default:\
								printf("ERROR: no code for step, index1 = 0, index2 = %u\n", (i+3));\
								return 0;\
						}\
						break;\
\
					case 6:\
						i = step_index2[type_index] - 2;\
						if( i > 3 )\
						{\
							printf("ERROR: no code for step, index1 = 1, index2 = %u\n", (i+2));\
							return 0;\
						}\
						chain_code <<= 2;\
						*code_length += 2;\
						chain_code += (u_int64_t)(i*16 + 0xF);\
						increment = 1;\
						break;\
\
					case 7: /* index1 == 2 & index2 == 3 */\
						chain_code <<= 2;\
						*code_length += 2;\
						chain_code += 0x2D;\
						increment = 1;\
						break;\
\
					case 8: /* index1 == 2 & index2 == 4 */\
						chain_code <<= 2;\
						*code_length += 2;\
						chain_code += 0x3D;\
						increment = 1;\
						break;\
\
					default:\
						printf("ERROR: bad step_type = %u\n", step_type[type_index]);\
						return 0;\
				} /* end switch */\
\
				type_index += increment;\
			}\
			while( type_index < type_count );\
\
		} /* end if len > lower_limit after top max continuation elements removed */\
	} /* end if len > lower_limit */\
\
	*code_length += 3;\
	chain_code <<= 3;\
\
/*  special 3-bit codes for the first two or three chain elements after 3 ([1,2,3] are implicit and not coded) */\
\
	switch( *current_partial_length )\
	{\
		case 2:\
			val3 = raw_c_list[*check_index].value;\
			val4 = val3 + 3;\
			val5 = val4 + val3;\
			break;\
\
		case 3:\
			val3 = working_chain[3].value;\
			val4 = raw_c_list[*check_index].value;\
			val5 = val4 + val3;\
			break;\
\
		case 4:\
			val3 = working_chain[3].value;\
			val4 = working_chain[4].value;\
			val5 = raw_c_list[*check_index].value;\
			break;\
\
		default:\
			val3 = working_chain[3].value;\
			val4 = working_chain[4].value;\
			val5 = working_chain[5].value;\
	}\
\
	switch( val3 )\
	{\
		case 4:\
		{\
			switch( val4 )\
			{\
				case 5:\
				{\
					chain_code += CHAIN_START_4_5;\
					return chain_code;\
				}\
				case 6:\
				{\
					chain_code += CHAIN_START_4_6;\
					return chain_code;\
				}\
				case 7:\
				{\
					chain_code += CHAIN_START_4_7;\
					return chain_code;\
				}\
				default:\
					goto report_error;\
			}\
		}\
		case 5:\
		{\
			switch( val4 )\
			{\
				case 6:\
				{\
					chain_code += CHAIN_START_5_6;\
					return chain_code;\
				}\
				case 7:\
				{\
					chain_code += CHAIN_START_5_7;\
					return chain_code;\
				}\
				case 8:\
				{\
					switch( val5 )\
					{\
						case 10:\
							chain_code += CHAIN_START_5_8_10;\
							return chain_code;\
						case 11:\
							chain_code += CHAIN_START_5_8_11;\
							return chain_code;\
						case 13:\
							chain_code += CHAIN_START_5_8_13;\
							return chain_code;\
						default:\
							goto report_error;\
					}\
				}\
				default:\
					goto report_error;\
			}\
		}\
		default:\
			goto report_error;\
	}\
\
	/* report unimplemented start sequence */\
	report_error:\
		printf("ERROR: unimplemented chain start sequence: %lu %lu %lu\n", val3, val4, val5);\
		return 0;\
}

