#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

#define MAX_P_FACTORS 10
#define MAX_POWERS 4
#define P_SEARCH_DONE (uint64_t)(-1)

/* structure for building arithmetic progressions */

typedef struct {
	uint32_t p;
	uint32_t max_power; /* The maximum power of p we precompute */
	uint32_t powers[MAX_POWERS]; /* powers[i] = p^(i+1), 0 <= i < max_power */
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
	uint64_t P, phiP;
} P_phiP_t;

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
sieve_add_aprog(aprog_list_t *list, uint32_t p) 
{
	uint32_t i;
	uint32_t power, power_limit;
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
		while (i < MAX_POWERS && power <= power_limit) {
			power *= p;
			i++;
		}
	}

	a->max_power = i;
	power = p;
	for (i = 0; i < a->max_power; i++) {
		a->powers[i] = power;
		power *= p;
	}
}

/*------------------------------------------------------------------------*/
const uint32_t factors[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,
				  47,53,59,61,67,71,73,79,83,89,97,101,
				  103,107,109,113,127,131,137,139,149,
				  151,157,163,167,173,179,181,191,193,
				  197,199,211,223,227,229,233,239,241,251};

/* Init aprog data for each prime p in factors[] with p <= factor_max */
void
aprog_init(aprog_list_t *aprog, uint32_t factor_min, uint32_t factor_max)
{
	uint32_t i, p;
	const uint32_t num_aprogs = sizeof(factors) / sizeof(uint32_t);

	aprog->num_aprogs_alloc = num_aprogs;
	aprog->aprogs = (aprog_t *)malloc(num_aprogs * sizeof(aprog_t));

	for (i = 0; i < num_aprogs; i++) {
		p = factors[i];

		if (p > factor_max)
			break;

		sieve_add_aprog(aprog, p);
	}
}

/*------------------------------------------------------------------------*/
void
sieve_fb_init(sieve_fb_t *s, uint32_t factor_min, uint32_t factor_max)
{
	memset(s, 0, sizeof(sieve_fb_t));

	aprog_init (&s->aprog_data, factor_min, factor_max);
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

double score_fudge (uint64_t p)
{
        /* lim inf phi(p) / p * log(log(p)) = e^-gamma. We compute 
           p/phi(p) * e^-gamma/log(log(p + 5)), so good p values will have a 
           value > 1. The +5 is another fudge term to avoid lots of bad 
           small p surviving */
        const double expgamma = 1.7810724179901979852365041031071795492;
        return expgamma * log(log((double)(p + 10)));
}

/*------------------------------------------------------------------------*/

uint64_t 
eulerphi_u64 (uint64_t n)
{
  uint64_t phi = 1UL, p;

  for (p = 2; p * p <= n; p += 2)
    {
      if (n % p == 0)
	{
	  phi *= p - 1;
	  n /= p;
	  while (n % p == 0)
	    {
	      phi *= p;
	      n /= p;
	    }
	}

      if (p == 2)
	p--;
    }

  /* now n is prime or 1 */

  return (n == 1) ? phi : phi * (n - 1);
}

/* Sort by increasing phiP, ties by increasing P */
int 
phiP_cmp (const void *a, const void* b) {
        const uint64_t aP = ((P_phiP_t *)a)->P;
        const uint64_t bP = ((P_phiP_t *)b)->P;
        const uint64_t aphiP = ((P_phiP_t *)a)->phiP;
        const uint64_t bphiP = ((P_phiP_t *)b)->phiP;
        
        return (aphiP > bphiP) ? 1 : (aphiP < bphiP) ? -1 : (aP > bP) ? 1 : (aP < bP) ? -1: 0;
}

int 
P_cmp (const void *a, const void* b) {
        const uint64_t aP = ((P_phiP_t *)a)->P;
        const uint64_t bP = ((P_phiP_t *)b)->P;
        
        return (aP > bP) ? 1 : (aP < bP) ? -1: 0;
}


void
remove_dupes (phi_score_heap_t *heap)
{
        P_phiP_t *PphiP;
        uint32_t i, written;
  
        PphiP = (P_phiP_t *) malloc (heap->num_entries * sizeof(P_phiP_t));

        for (i = 0; i < heap->num_entries; i++) {
                PphiP[i].P = heap->entries[i].p;
                PphiP[i].phiP = eulerphi_u64(heap->entries[i].p);
        }

        qsort (PphiP, heap->num_entries, sizeof(P_phiP_t), phiP_cmp);        

        for (i = 1, written = 0; i < heap->num_entries; i++) {
                if (PphiP[i - 1].phiP > PphiP[i].phiP ||
                        (PphiP[i - 1].phiP == PphiP[i].phiP && PphiP[i - 1].P > PphiP[i].P))
                    abort();
        }
        
        for (i = 1, written = 0; i < heap->num_entries; i++) {
                if (PphiP[i - 1].phiP != PphiP[i].phiP)
                  PphiP[written++] = PphiP[i - 1];
        }
        PphiP[written++] = PphiP[i - 1];

        qsort (PphiP, written, sizeof(P_phiP_t), P_cmp);        

        for (i = 0; i < written; i++) {
                const uint64_t P = PphiP[i].P;
                const uint64_t phiP = PphiP[i].phiP;
                const double ratio = (double) P / (double) PphiP[i].phiP;
		printf("%" PRIu64 " %" PRIu64 " %.10f %.10f\n", P, phiP, ratio, ratio / score_fudge(P));
        }
}

/*------------------------------------------------------------------------*/
int main(int argc, char **argv) {

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
                                 (double) s.p_enum.curr_phi[s.p_enum.num_factors]
                                 / score_fudge(p));
	}

	remove_dupes (&heap);

	return 0;
}
