
#include "libntt.h"

#define CACHE_LINE_SIZE 64

void *
sp_aligned_malloc (size_t len)
{
  void *ptr, *aligned_ptr;
  size_t addr;

  ptr = malloc (len + CACHE_LINE_SIZE);
  if (ptr == NULL)
    return NULL;

  addr = (size_t)ptr;				
  addr = CACHE_LINE_SIZE - (addr % CACHE_LINE_SIZE);
  aligned_ptr = (void *)((char *)ptr + addr);

  *( (void **)aligned_ptr - 1 ) = ptr;
  return aligned_ptr;
}

void
sp_aligned_free (void *newptr) 
{
  void *ptr;

  if (newptr == NULL) 
    return;
  ptr = *( (void **)newptr - 1 );
  free (ptr);
}
