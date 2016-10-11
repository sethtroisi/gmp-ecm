#define _GNU_SOURCE
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
/* For EAGAIN etc. */
#include <errno.h>
#endif
#include "listz_handle.h"

/* #define TRACE_ITER yes */

/* Init a listz_handle_t to store up to len residues (modulo m). 
   If filename != NULL, uses disk storage, otherwise memory.
   Returns NULL if something goes wrong (i.e., if a memory allocation
   or opening a file fails) */

listz_handle_t 
listz_handle_init (const char *filename, const uint64_t len, const mpz_t m)
{
  listz_handle_t F;
  void *buf;

  F = malloc (sizeof (_listz_handle_t));
  if (F == NULL)
    return NULL;
  
  /* Find out how many file_word_t's m has */
  buf = (file_word_t *) mpz_export (NULL, &F->words, -1, sizeof(file_word_t), 
                                   -1, 0, m);
  if (buf == NULL)
    {
      free (F);
      return NULL;
    }
  free(buf);

  F->len = len;
  if (filename == NULL)
    {
      F->storage = 0; /* Memory storage */
      F->data.mem = init_list2 (len, mpz_sizeinbase(m, 2));
      if (F->data.mem == NULL)
        {
          free (F);
          F = NULL;
        }
    } else {
      F->storage = 1; /* Disk storage */
      F->filename = (char *) malloc ((strlen (filename) + 1) * sizeof(char));
      if (F->filename == NULL)
        {
          free (F);
          return NULL;
        }
      strcpy (F->filename, filename);
      F->data.file = fopen (F->filename, "rb+");
      if (F->data.file == NULL)
        F->data.file = fopen (F->filename, "wb+");
      if (F->data.file == NULL)
        {
          free (F->filename);
          free (F);
          return NULL;
        }
#ifdef HAVE_FALLOCATE
      fallocate (fileno(F->data.file), 0, (off_t) 0, 
                 F->words * sizeof(file_word_t) * len);
#endif
#if defined(HAVE_SETVBUF) && defined(HAVE_AIO_READ)
      /* Set to unbuffered mode as we use aio_*() functions for reading
         in the background */
      setvbuf (F->data.file, NULL, _IONBF, 0);
#endif
    }

  return F;
}


listz_handle_t 
listz_handle_from_listz (const listz_t l, const uint64_t len, const mpz_t m)
{
  listz_handle_t F;
  void *buf;

  F = malloc (sizeof (_listz_handle_t));
  if (F == NULL)
    return NULL;
  
  /* Find out how many file_word_t's m has */
  buf = (file_word_t *) mpz_export (NULL, &F->words, -1, sizeof(file_word_t), 
                                   -1, 0, m);
  if (buf == NULL)
    {
      free (F);
      return NULL;
    }
  free(buf);

  F->len = len;
  F->storage = 0; /* Memory storage */
  F->data.mem = l;

  return F;
}


void 
listz_handle_clear (listz_handle_t F)
{
  if (F->storage == 0)
    {
      clear_list (F->data.mem, F->len);
      F->data.mem = NULL;
    }
  else
    {
      fclose (F->data.file);
      F->data.file = NULL;
      remove (F->filename);
      free (F->filename);
    }
  free (F);
}


static inline void 
export_residue (file_word_t *buf, const size_t bufsize, const mpz_t r)
{
  size_t nr;

  /* Export r to buf */
  mpz_export (buf, &nr, -1, sizeof(file_word_t), 0, 0, r);
  ASSERT_ALWAYS (nr <= bufsize);

  /* Pad buf with zeroes */
  for ( ; nr < bufsize; nr++)
    buf[nr] = 0;
}


static inline int 
listz_handle_seek_entry (listz_handle_t F, const uint64_t index)
{
  int64_t foffset;

  ASSERT (F->storage == 1);
  ASSERT (index <= INT64_MAX);
  ASSERT (index <= INT64_MAX / sizeof(file_word_t) / F->words);
  foffset = (int64_t) index * F->words * sizeof(file_word_t);
  return aux_fseek64 (F->data.file, foffset, SEEK_SET);
}


/* Output a polynomial of degree len-1, or a monic polynomial of degree len.
   In either case, len is the number of coefficients read from "l".
   If symmetric == 1, then the polynomial is printed as a reciprocal Laurent
   polynomial where the coefficients stored in l (and perhaps the leading 
   monomial) are in standard basis. */
   
void
listz_handle_output_poly (const listz_handle_t l, const uint64_t len, 
                          const int monic, const int symmetric, 
                          const char *prefix, const char *suffix, 
                          const int verbosity)
{
  uint64_t i;
  mpz_t m;
  listz_iterator_t *iter;

  if (!test_verbose(verbosity))
    return;
  
  if (prefix != NULL)
    outputf (verbosity, prefix);

  if (len == 0)
    {
      if (monic)
	outputf (verbosity, "1\n");
      else
	outputf (verbosity, "0\n");
      return;
    }

  mpz_init (m);
  iter = listz_iterator_init (l, 0);
  for (i = 0; i < len; i++)
    {
      const uint64_t deg = len - ((monic != 0) ? 0U : 1U);
      const char *plus = (i < deg) ? " + " : "";
      listz_iterator_read (iter, m);
      if (symmetric)
        outputf (verbosity, "Mod(%Zd,N) * (x^%" PRIu64 " + x^-%" PRIu64 ")%s", 
                 m, i, i, plus);
      else
        outputf (verbosity, "Mod(%Zd,N) * x^%" PRIu64 "%s", m, i, plus);
    }
  if (monic)
    {
      if (symmetric)
	outputf (verbosity, "(x^%" PRIu64 " + x^-%" PRIu64 ")", len, len);
      else
	outputf (verbosity, "x^%" PRIu64, len);
    }
  listz_iterator_clear (iter);
  if (suffix != NULL)
    outputf (verbosity, suffix);
  mpz_clear (m);
}


/* Iterator functions for sequential access to elements of a
   list_handle_t. */

static void 
listz_iterator_fetch (listz_iterator_t *iter, const uint64_t offset)
{
  ASSERT (iter->handle->storage == 1);
  iter->offset = offset;
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
  iter->cb.aio_offset = (off_t) offset * iter->handle->words * sizeof(file_word_t);
  iter->cb.aio_buf = iter->buf[iter->active_buffer];
  iter->cb.aio_nbytes = iter->bufsize * iter->handle->words * sizeof(file_word_t);
  {
    int r = aio_read (&iter->cb);
    if (r != 0)
      {
        fprintf (stderr, "%s(): aio_read() returned %d\n", __func__, r);
        abort ();
      }
  }
#else /* ifdef HAVE_AIO_READ */
#ifdef _OPENMP
#pragma omp critical 
#endif
  {
    listz_handle_seek_entry (iter->handle, iter->offset);
    iter->valid = fread (iter->buf, iter->handle->words * sizeof(file_word_t), 
                         iter->bufsize, iter->handle->data.file);
  }
#endif /* ifdef HAVE_AIO_READ else */
}


#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
static size_t 
listz_iterator_suspend (struct aiocb * const cb)
{
  int r;
  ssize_t s;
  
  do {
    const struct aiocb * aiocb_list[1] = {cb};
    r = aio_suspend (aiocb_list, 1, NULL);
  } while (r == EAGAIN || r == EINTR);
  if (r != 0)
    {
      fprintf (stderr, "%s(): aio_suspend() returned error, errno = %d\n",
               __func__, errno);
      abort ();
    }

  s = aio_return (cb);
  if (s < 0)
    {
      fprintf (stderr, "%s(): aio_return() returned error code %ld\n", 
               __func__, (long int) s);
      abort();
    }
  return (size_t) s;
}
#endif  


static void 
listz_iterator_flush (listz_iterator_t *iter)
{

  ASSERT (iter->handle->storage == 1);
  if (iter->writeptr == 0)
    return;

#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
  {
    size_t nbytes, written;
    int r;

    iter->cb.aio_offset = 
        (off_t) iter->offset * iter->handle->words * sizeof(file_word_t);
    iter->cb.aio_buf = iter->buf[iter->active_buffer];
    nbytes = iter->writeptr * iter->handle->words * sizeof(file_word_t);
    iter->cb.aio_nbytes = nbytes;
    r = aio_write (&iter->cb);
    if (r != 0)
      {
        fprintf (stderr, "%s(): aio_write() returned error, errno = %d\n", 
                 __func__, errno);
        abort ();
      }
    written = listz_iterator_suspend (&iter->cb);
    ASSERT_ALWAYS (written == nbytes);
  }
#else
#ifdef _OPENMP
#pragma omp critical 
#endif
  {
    size_t written;
    listz_handle_seek_entry (iter->handle, iter->offset);
    written = fwrite (iter->buf, sizeof(file_word_t) * iter->handle->words, 
                      iter->writeptr, iter->handle->data.file);
    ASSERT_ALWAYS (written == iter->writeptr);
  }
#endif
  iter->writeptr = 0;
}


listz_iterator_t *  
listz_iterator_init2 (listz_handle_t h, const uint64_t firstres, 
                      const size_t nr_buffered)
{
  listz_iterator_t *iter;
  
  iter = (listz_iterator_t *) malloc (sizeof(listz_iterator_t));
  if (iter == NULL)
    return NULL;
  iter->handle = h;

  if (iter->handle->storage == 0)
    iter->readptr = iter->writeptr = (size_t) firstres;
  else
    {
      iter->offset = firstres;
      iter->readptr = iter->writeptr = iter->valid = 0;
      iter->bufsize = nr_buffered;
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
      iter->buf[0] = malloc (iter->bufsize * iter->handle->words * 
                             sizeof(file_word_t));
      iter->buf[1] = malloc (iter->bufsize * iter->handle->words * 
                             sizeof(file_word_t));
      if (iter->buf[0] == NULL || iter->buf[1] == NULL) 
        {
          free (iter->buf[0]);
          free (iter->buf[1]);
          free (iter);
          return NULL;
        }
#ifdef _OPENMP
#pragma omp critical 
#endif
      {
        /* Prevent other access to the file which would lead to data 
           corruption */
        iter->hidden_file = h->data.file;
        h->data.file = NULL;
      }
      ASSERT_ALWAYS (iter->hidden_file != NULL);
      iter->active_buffer = 0;
      memset (&iter->cb, 0, sizeof(struct aiocb));
      iter->cb.aio_fildes = fileno (iter->hidden_file);
      iter->cb.aio_reqprio = 0;
      iter->cb.aio_sigevent.sigev_notify = SIGEV_NONE;
#else
      iter->buf = malloc (iter->bufsize * iter->handle->words * 
                          sizeof(file_word_t));
      if (iter->buf == NULL) 
        {
          free (iter);
          return NULL;
        }
#endif
    }

  return iter;
}


listz_iterator_t *  
listz_iterator_init (listz_handle_t h, const uint64_t firstres)
{
  const size_t listz_iterator_nr_buffered = 1<<16;
  return listz_iterator_init2 (h, firstres, listz_iterator_nr_buffered);
}


void  
listz_iterator_clear (listz_iterator_t *iter)
{
  if (iter->handle->storage == 1)
    {
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
      if (iter->valid > 0)
        {
          listz_iterator_suspend (&iter->cb);
        }
#endif
      listz_iterator_flush (iter);
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
      iter->handle->data.file = iter->hidden_file;
      free (iter->buf[0]);
      free (iter->buf[1]);
#else
      free (iter->buf);
#endif
    }
  free (iter);
}


/* Outside of listz_iterator_newbuf() we have:
   If iter->valid == 0, then there is no outstanding read request
   If iter->valid != 0, then there is an outstanding read request for the 
     non-active buffer
   There is no outstanding write request
   If there is an outstanding read request, then iter->offset is the start 
     position in the file (in units of residues) of that read request
*/

static inline void 
listz_iterator_switchbuf (listz_iterator_t *iter)
{
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
  if (iter->valid == 0)
    {
      /* We never read data before. We have to fetch data for the 
         current buffer and wait for it to finish, then issue a fetch 
         for the next buffer */
      listz_iterator_fetch (iter, iter->offset);
      iter->active_buffer ^= 1;
    }
  iter->valid = listz_iterator_suspend (&iter->cb) 
      / sizeof(file_word_t) / iter->handle->words;
#endif
  listz_iterator_flush (iter);
  listz_iterator_fetch (iter, iter->offset + iter->valid);
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
  iter->active_buffer ^= 1;
#endif
  iter->readptr = 0;
  ASSERT_ALWAYS (iter->valid > 0);
}


void
listz_iterator_read (listz_iterator_t *iter, mpz_t r)
{
  if (iter->handle->storage == 0)
    {
      mpz_set (r, iter->handle->data.mem[iter->readptr]);
#if defined(TRACE_ITER)
      gmp_printf ("%s(): readptr = %" PRIu64 ", r = %Zd (in memory)\n", 
                  __func__, (uint64_t) iter->readptr, r);
#endif
    }
  else
    {
      /* Try to detect incorrect use of iterator. We allow either read-only, 
         in which case we must have writeptr == 0 at all times, or sequential
         update (read-then-write) of each residue, in which case we must have
         writeptr == readptr here */
      ASSERT (iter->writeptr == 0 || iter->readptr == iter->writeptr);

      if (iter->readptr == iter->valid)
        listz_iterator_switchbuf (iter);

#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
      mpz_import (r, iter->handle->words, -1, sizeof(file_word_t), 0, 0, 
                  &iter->buf[iter->active_buffer][iter->readptr * iter->handle->words]);
#else
      mpz_import (r, iter->handle->words, -1, sizeof(file_word_t), 0, 0, 
                  &iter->buf[iter->readptr * iter->handle->words]);
#endif
#if defined(TRACE_ITER)
      gmp_printf ("%s(): offset = %" PRIu64 ", readptr = %" PRIu64 
                   ", r = %Zd (on disk)\n", 
                  __func__, iter->offset, (uint64_t) iter->readptr, r);
#endif
    }
  iter->readptr++;
}


void
listz_iterator_write (listz_iterator_t *iter, const mpz_t r)
{
  if (iter->handle->storage == 0)
    {
#if defined(TRACE_ITER)
      gmp_printf ("%s(): writeptr = %"PRIu64", r = %Zd\n", 
                  __func__, (uint64_t) iter->writeptr, r);
#endif
      mpz_set (iter->handle->data.mem[iter->writeptr], r);
    }
  else
    {
      /* Try to detect incorrect use of iterator. We allow either write-only, 
         in which case we must have readptr == 0 at all times, or sequential
         update (read-then-write) of each residue, in which case we must have
         writeptr + 1 == readptr */
      ASSERT (iter->readptr == 0 || iter->writeptr + 1 == iter->readptr);
      ASSERT (iter->writeptr <= iter->bufsize);
#if defined(TRACE_ITER)
      gmp_printf ("%s(): offset = %"PRIu64", writeptr = %"PRIu64
                  ", r = %Zd (on disk)\n", 
                  __func__, iter->offset, (uint64_t) iter->writeptr, r);
#endif
      if (iter->writeptr == iter->bufsize)
        {
#if defined(TRACE_ITER)
          printf ("%s(): flushing %"PRIu64" entries from buffer %d\n", 
                  iter->writeptr, iter->active_buffer);
#endif
          listz_iterator_flush (iter);
          iter->offset += iter->bufsize;
        }
      ASSERT_ALWAYS (mpz_sgn (r) >= 0);
      /* TODO: we may want to allow residues that are not fully reduced 
         (mod modulus), but only as far as the ECRT reduces them. */
#if defined(HAVE_AIO_READ) && defined(WANT_AIO)
      export_residue (&iter->buf[iter->active_buffer][iter->writeptr * iter->handle->words], 
                      iter->handle->words, r);
#else
      export_residue (&iter->buf[iter->writeptr * iter->handle->words], 
                      iter->handle->words, r);
#endif
    }
  iter->writeptr++;
}

/* Functions that can be used as callbacks to listz_iterator_read() and 
   listz_iterator_write() in mpzspv_fromto_mpzv(). Note that calling, e.g., 
   listz_iterator_read() by de-referencing a pointer of type mpz_producerfunc_t
   leads to undefined program behavior according to the C standard, even though
   it happens to work fine on x86[_64] architectures at least. On other 
   architectures, a function pointer may carry contextual information which 
   could be incorrect when de-referencing a function pointer of the wrong type. 
*/
void
listz_iterator_read_callback (void *iter, mpz_t r)
{
  listz_iterator_read ((listz_iterator_t *) iter, r);
}

void
listz_iterator_write_callback (void *iter, const mpz_t r)
{
  listz_iterator_write ((listz_iterator_t *) iter, r);
}
