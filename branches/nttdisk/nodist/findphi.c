#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define MAX_P_FACTORS 10
#define MAX_POWERS 4
#define P_SEARCH_DONE (uint64_t)(-1)

/* structure for building arithmetic progressions */

typedef struct {
	uint32_t p;
	uint32_t max_power;
	uint32_t powers[MAX_POWERS];
	uint64_t cofactor_max;
} aprog_t;

typedef struct {
	aprog_t *aprogs;
	uint32_t num_aprogs;
	uint32_t num_aprogs_alloc;
} aprog_list_t;

/* structures for finding arithmetic progressions via
   explicit enumeration */

typedef struct {
	uint32_t num_factors;
	uint32_t curr_factor[MAX_P_FACTORS + 1];
	uint32_t curr_power[MAX_P_FACTORS + 1];
	uint64_t curr_phi[MAX_P_FACTORS + 1];
	uint64_t curr_prod[MAX_P_FACTORS + 1];
} p_enum_t;

/* a factory for building arithmetic progressions */

typedef struct {
	uint64_t p_min, p_max;
	aprog_list_t aprog_data;
	p_enum_t p_enum;
} sieve_fb_t;


typedef struct {
	double score;
	uint64_t p;
} phi_score_t;

typedef struct {
	uint32_t num_entries;
	uint32_t max_entries;
	phi_score_t *entries;
} phi_score_heap_t;

/*------------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void
heapify(phi_score_t *h, uint32_t index, uint32_t size) {

	uint32_t c;
	phi_score_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].score > h[c+1].score)
			c++;

		if (h[index].score > h[c].score) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].score > h[c].score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void
make_heap(phi_score_t *h, uint32_t size) {

	int32_t i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32_t)i, size);
}

static void
save_phi_score(phi_score_heap_t *heap, uint64_t p, double score)
{
	phi_score_t *h = heap->entries;

	if (heap->num_entries <= heap->max_entries - 1) {
		phi_score_t *s = h + heap->num_entries++;
		s->p = p;
		s->score = score;
		if (heap->num_entries == heap->max_entries)
			make_heap(h, heap->max_entries);
	}
	else if (h->score < score) {
		h->p = p;
		h->score = score;
		heapify(h, 0, heap->max_entries);
	}
}

/*------------------------------------------------------------------------*/
static void
sieve_add_aprog(sieve_fb_t *s, uint32_t p) 
{
	uint32_t i;
	uint32_t power, power_limit;
	aprog_list_t *list = &s->aprog_data;
	aprog_t *a;

	if (list->num_aprogs == list->num_aprogs_alloc) {
		list->num_aprogs_alloc *= 2;
		list->aprogs = (aprog_t *)realloc(list->aprogs,
						list->num_aprogs_alloc *
						sizeof(aprog_t));
	}

	a = list->aprogs + list->num_aprogs;
	a->p = p;

	list->num_aprogs++;

	power = p;
	power_limit = (uint32_t)(-1) / p;
	i = 1;
	if (p <= 19) {
		while (i < MAX_POWERS && power < power_limit) {
			power *= p;
			i++;
		}
	}

	a->max_power = i;
	a->powers[0] = p;
	power = p;
	for (i = 1; i < a->max_power; i++) {
		power *= p;
		a->powers[i] = power;
	}
}

/*------------------------------------------------------------------------*/
const uint32_t factors[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,
				  47,53,59,61,67,71,73,79,83,89,97,101,
				  103,107,109,113,127,131,137,139,149,
				  151,157,163,167,173,179,181,191,193,
				  197,199,211,223,227,229,233,239,241,251};
void
sieve_fb_init(sieve_fb_t *s, uint32_t factor_min, uint32_t factor_max)
{
	uint32_t i, p;
	aprog_list_t *aprog = &s->aprog_data;

	memset(s, 0, sizeof(sieve_fb_t));

	aprog->num_aprogs_alloc = 100;
	aprog->aprogs = (aprog_t *)malloc(aprog->num_aprogs_alloc *
						sizeof(aprog_t));

	for (i = 0; i < sizeof(factors) / sizeof(uint32_t); i++) {
		p = factors[i];

		if (p > factor_max)
			break;

		sieve_add_aprog(s, p);
	}
}

/*------------------------------------------------------------------------*/
void 
sieve_fb_reset(sieve_fb_t *s, uint64_t p_min, uint64_t p_max)
{
	uint32_t i;
	aprog_t *aprogs = s->aprog_data.aprogs;
	uint32_t num_aprogs = s->aprog_data.num_aprogs;

	s->p_min = p_min;
	s->p_max = p_max;

	/* set up for finding arithmetic progressions by
	   enumerating combinations of small factors p_i
	   from the aprog list */

	if (num_aprogs > 0) {
		p_enum_t *p_enum = &s->p_enum; 

		/* clear the enum state */

		p_enum->curr_factor[0] = 0;
		p_enum->num_factors = 0;
		p_enum->curr_power[0] = 0;
		p_enum->curr_prod[0] = 1;
		p_enum->curr_phi[0] = 1;

		for (i = 0; i < num_aprogs; i++) {
			aprog_t *a = aprogs + i;

			a->cofactor_max = p_max / a->p;
		}
	}
}

/*------------------------------------------------------------------------*/
static uint64_t
get_next_enum(sieve_fb_t *s)
{
	p_enum_t *p_enum = &s->p_enum;
	uint32_t *curr_factor = p_enum->curr_factor;
	uint32_t *curr_power = p_enum->curr_power;
	uint64_t *curr_phi = p_enum->curr_phi;
	uint64_t *curr_prod = p_enum->curr_prod;

	while (1) {

		uint32_t i = p_enum->num_factors;
		uint32_t power_up = (i && curr_factor[i] == curr_factor[i - 1]);
		aprog_t *a = NULL;

		if (curr_factor[i] < s->aprog_data.num_aprogs)
			a = s->aprog_data.aprogs + curr_factor[i];

		if (a != NULL && curr_prod[i] <= a->cofactor_max
		    && !(power_up && ++curr_power[i - 1] >= a->max_power)
		    && (power_up || i < MAX_P_FACTORS)) {

			uint64_t p = curr_prod[i] * a->p;
                        uint64_t phi = curr_phi[i];

			if (power_up) {
				phi *= a->p;
                        }
                        else {
				phi *= a->p - 1;
				p_enum->num_factors = ++i;
			}

			curr_prod[i] = p;
			curr_phi[i] = phi;
			curr_factor[i] = 0;
			curr_power[i] = 0;

			if (p >= s->p_min)
				return p;
		}
		else if (i) {

			i--;
			curr_factor[i]++;
			curr_power[i] = 0;
			p_enum->num_factors = i;
		}
		else {
			return P_SEARCH_DONE;
		}
	}
}

/*------------------------------------------------------------------------*/
int main(int argc, char **argv) {

        uint32_t i;
	sieve_fb_t s;
        phi_score_heap_t heap;

        if (argc != 5) {
                printf("usage: %s <start> <end> <fb_limit> <heap_size>\n",
				argv[0]);
		return 0;
	}

	sieve_fb_init(&s, 2, atol(argv[3]));

        memset(&heap, 0, sizeof(heap));
        heap.max_entries = atol(argv[4]);
        heap.entries = (phi_score_t *)malloc(heap.max_entries *
                                              sizeof(phi_score_t));

	sieve_fb_reset(&s, strtoull(argv[1], NULL, 0), 
			strtoull(argv[2], NULL, 0));

	while (1) {
		uint64_t p = get_next_enum(&s);
		if (p == P_SEARCH_DONE)
			break;
                save_phi_score(&heap, p, (double)p / 
                                 s.p_enum.curr_phi[s.p_enum.num_factors]);
	}

        for (i = 0; i < heap.num_entries; i++)
		printf("%" PRIu64 " %lf\n", heap.entries[i].p,
		      			heap.entries[i].score);

	return 0;
}
