/* Memory allocation used during tests.

Copyright 2001, 2002 Free Software Foundation, Inc.

This file was copied from the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>  /* for abort */
#include <gmp.h>
#include "ecm.h"

void *__gmp_default_allocate (size_t);
void *__gmp_default_reallocate (void *, size_t, size_t);
void __gmp_default_free (void *, size_t);

/* Each block allocated is a separate malloc, for the benefit of a redzoning
   malloc debugger during development or when bug hunting.

   Sizes passed when reallocating or freeing are checked (the default
   routines don't care about these).

   Memory leaks are checked by requiring that all blocks have been freed
   when tests_memory_end() is called.  Test programs must be sure to have
   "clear"s for all temporary variables used.  */

#define NAME_LEN 8

struct header {
  void           *ptr;
  size_t         size;
  char           name[NAME_LEN];
  unsigned int   line;
  struct header  *next;
};

struct header  *tests_memory_list = NULL;

static unsigned long nr_realloc = 0, nr_realloc_move = 0;
static char cur_name[NAME_LEN];
static unsigned int cur_line;
static unsigned long cur_mem, peak_mem;

/* Return a pointer to a pointer to the found block (so it can be updated
   when unlinking). */
static struct header **
tests_memory_find (void *ptr)
{
  struct header  **hp;

  for (hp = &tests_memory_list; *hp != NULL; hp = &((*hp)->next))
    if ((*hp)->ptr == ptr)
      return hp;

  return NULL;
}

#if 0
static int
tests_memory_valid (void *ptr)
{
  return (tests_memory_find (ptr) != NULL);
}
#endif

static void *
tests_allocate (size_t size)
{
  struct header  *h;
  int i;

  if (size == 0)
    {
      printf ("tests_allocate(): attempt to allocate 0 bytes\n");
      abort ();
    }

  if (cur_name[0] == 0)
    cur_name[1] = 0; /* Set breakpoint here to find untagged allocs */

  h = (struct header *) __gmp_default_allocate (sizeof (*h));
  h->next = tests_memory_list;
  tests_memory_list = h;

  h->size = size;
  h->ptr = (struct header*) __gmp_default_allocate (size);
  for (i = 0; i < NAME_LEN; i++)
    h->name[i] = cur_name[i];
  h->line = cur_line;
  cur_mem += size;
  if (cur_mem > peak_mem)
    peak_mem = cur_mem;
  return h->ptr;
}

static void *
tests_reallocate (void *ptr, size_t old_size, size_t new_size)
{
  struct header  **hp, *h;

  if (new_size == 0)
    {
      printf ("tests_reallocate(): attempt to reallocate 0x%lX to 0 bytes\n",
              (unsigned long) ptr);
      abort ();
    }

  hp = tests_memory_find (ptr);
  if (hp == NULL)
    {
      printf ("tests_reallocate(): attempt to reallocate bad pointer 0x%lX\n",
              (unsigned long) ptr);
      abort ();
    }
  h = *hp;

  if (h->size != old_size)
    {
      printf ("tests_reallocate(): bad old size %zd, should be %zd\n",
              old_size, h->size);
      abort ();
    }

  if (h->size > cur_mem)
    {
      printf ("tests_reallocate(): h->size = %zd but cur_mem = %lu\n",
              h->size, cur_mem);
      abort();
    }

  cur_mem = cur_mem - h->size + new_size;
  if (cur_mem > peak_mem)
    peak_mem = cur_mem;

#if 0
  printf ("Reallocating %p, first allocated in %s, line %d, from %d to %d\n",
          ptr, h->name, h->line, h->size, new_size);
  if (new_size <= h->size)
    printf ("Unnecessary realloc!\n");
#endif
  nr_realloc++;
  h->size = new_size;
  h->ptr = (struct header*) __gmp_default_reallocate (ptr, old_size, new_size);
  if (h->ptr != ptr)
    nr_realloc_move++;
  return h->ptr;
}

static struct header **
tests_free_find (void *ptr)
{
  struct header  **hp = tests_memory_find (ptr);
  if (hp == NULL)
    {
      printf ("tests_free(): attempt to free bad pointer 0x%lX\n",
              (unsigned long) ptr);
      abort ();
    }
  return hp;
}

static void
tests_free_nosize (void *ptr)
{
  struct header  **hp = tests_free_find (ptr);
  struct header  *h = *hp;

  if (h->size > cur_mem)
    {
      printf ("tests_free_nosize(): h->size = %zd but cur_mem = %lu\n",
              h->size, cur_mem);
      abort();
    }

  cur_mem -= h->size;

  *hp = h->next;  /* unlink */

  __gmp_default_free (ptr, h->size);
  __gmp_default_free (h, sizeof (*h));
}

void
tests_free (void *ptr, size_t size)
{
  struct header  **hp = tests_free_find (ptr);
  struct header  *h = *hp;

  if (h->size != size)
    {
      printf ("tests_free(): bad size %zd, should be %zd\n", size, h->size);
      abort ();
    }

  tests_free_nosize (ptr);
}

void
tests_memory_start (void)
{
  mp_set_memory_functions (tests_allocate, tests_reallocate, tests_free);
  cur_name[0] = 0;
  cur_line = 0;
  cur_mem = 0L;
  peak_mem = 0L;
}

void
tests_memory_reset (void)
{
  mp_set_memory_functions (__gmp_default_allocate, __gmp_default_reallocate,
                           __gmp_default_free);
}

void
tests_memory_end (void)
{
  if (tests_memory_list != NULL)
    {
      struct header  *h;
      unsigned  count;

      printf ("tests_memory_end(): not all memory freed\n");

      count = 0;
      for (h = tests_memory_list; h != NULL; h = h->next)
        {
	  count++;
	  printf ("Memory at %p, allocated by %s, line %d\n", 
	          h->ptr, h->name, h->line);
	}

      printf ("    %u block(s) remaining\n", count);
      abort ();
    }
  
  if (cur_mem != 0)
    {
      printf ("tests_memory_end(): cur_mem = %lu but list of allocated "
              "memory empty\n", cur_mem);
      abort ();
    }
  
  printf ("%lu reallocates, %lu reallocates with move, peak_mem = %lu\n", 
          nr_realloc, nr_realloc_move, peak_mem);
}

void
tests_memory_status (void)
{
  unsigned count = 0, size = 0;

  if (tests_memory_list != NULL)
    {
      struct header  *h;

      for (h = tests_memory_list; h != NULL; h = h->next)
        {
          count++;
          size += h->size;
        }

    }

  if (size != cur_mem)
    {
      printf ("tests_memory_status(): size = %d but cur_mem = %lu", 
              size, cur_mem);
      abort();
    }

  printf ("    %u blocks remaining, total size %u\n", count, size);
}

void tests_memory_set_location (char *name, unsigned int line)
{
  unsigned int i;

  for (i = 0; i < NAME_LEN; i++)
    cur_name[i] = name[i];
  cur_line = line;
}

