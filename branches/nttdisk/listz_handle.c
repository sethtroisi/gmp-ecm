#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#include "listz_handle.h"

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
  
  /* Find out how many file_word_t's  m has */
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
      F->data.file = fopen (F->filename, "w+");
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
    }

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


/* Fetches one entry from F (either in memory or file) and stores it in r. */

void
listz_handle_get (listz_handle_t F, mpz_t r, file_word_t *buf, 
    const size_t index)
{
  if (F->storage == 0)
    mpz_set (r, F->data.mem[index]);
  else
    {
      size_t nr;

#ifdef _OPENMP
#pragma omp critical 
#endif
      {
        fseek (F->data.file, sizeof(file_word_t) * F->words * index, SEEK_SET);
        nr = fread (buf, sizeof(file_word_t), F->words, F->data.file);
      }

      ASSERT_ALWAYS (nr == F->words);
      mpz_import (r, F->words, -1, sizeof(file_word_t), 0, 0, buf);
    }
}

void
listz_handle_get2 (listz_handle_t F, mpz_t r, const size_t index)
{
  file_word_t *buf = NULL;
  if (F->storage == 1)
    buf = malloc (F->words * sizeof (file_word_t));
  listz_handle_get (F, r, buf, index);
  free(buf);
}


/* Stores the value of r in an entry of F (either in memory or file) */

void
listz_handle_set (listz_handle_t F, const mpz_t r, file_word_t *buf,
    const size_t index)
{
  if (F->storage == 0)
    mpz_set (F->data.mem[index], r);
  else
    {
      size_t nr;

      ASSERT_ALWAYS (mpz_sgn (r) >= 0);
      export_residue (buf, F->words, r);
#ifdef _OPENMP
#pragma omp critical 
#endif
      {
        fseek (F->data.file, sizeof(file_word_t) * F->words * index, SEEK_SET);
        nr = fwrite (buf, sizeof(file_word_t), F->words, F->data.file);
      }
      ASSERT_ALWAYS (nr == F->words);
    }
}


void
listz_handle_output_poly (const listz_handle_t l, const uint64_t len, 
                          const int monic, const int symmetric, 
                          const char *prefix, const char *suffix, 
                          const int verbosity)
{
  uint64_t i;
  mpz_t m;
  file_word_t *buf = NULL;

  if (!test_verbose(verbosity))
    return;
  
  if (l->storage != 0)
    buf = (file_word_t *) malloc (l->words * sizeof(file_word_t));

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
  if (monic)
    {
      if (symmetric)
	outputf (verbosity, "(x^%" PRIu64 " + x^-%" PRIu64 ") + ", len, len);
      else
	outputf (verbosity, "x^%" PRIu64 " + ", len);
    }
  for (i = len - 1; i > 0; i--)
    {
      listz_handle_get (l, m, buf, i);
      if (symmetric)
        outputf (verbosity, "Mod(%Zd,N) * (x^%" PRIu64 " + x^-%" PRIu64 ") + ", 
                 m, i, i);
      else
        outputf (verbosity, "Mod(%Zd,N) * x^%" PRIu64 " + ", m, i);
    }
  listz_handle_get (l, m, buf, 0);
  outputf (verbosity, "Mod(%Zd,N)", m);
  if (suffix != NULL)
    outputf (verbosity, suffix);
  free (buf);
  buf = NULL;
  mpz_clear (m);
}

static void 
listz_iterator_fetch (listz_iterator_t *iter, size_t offset)
{
  iter->offset = offset;
  ASSERT (iter->handle->storage == 1);
#ifdef _OPENMP
#pragma omp critical 
#endif
  {
    fseek (iter->handle->data.file, offset * sizeof(file_word_t), SEEK_SET);
    iter->valid = fread (iter->buf, sizeof(file_word_t), iter->bufsize, 
                         iter->handle->data.file);
  }
  if (iter->valid != iter->bufsize)
    {
      /* Check we read an integral number of residues */
      ASSERT_ALWAYS (iter->valid % iter->handle->words == 0);
    }
}


static void 
listz_iterator_flush (listz_iterator_t *iter)
{
  size_t written;

  ASSERT (iter->handle->storage == 1);
  if (!iter->dirty)
    return;

#ifdef _OPENMP
#pragma omp critical 
#endif
  {
    fseek (iter->handle->data.file, iter->offset * sizeof(file_word_t), 
           SEEK_SET);
    written = fwrite (iter->buf, sizeof(file_word_t), iter->writeptr, 
                      iter->handle->data.file);
  }
  ASSERT_ALWAYS (written == iter->writeptr);
  iter->dirty = 0;
}


listz_iterator_t *  
listz_iterator_init2 (listz_handle_t h, const size_t firstres, 
                      const size_t nr_buffered)
{
  listz_iterator_t *iter;
  
  iter = (listz_iterator_t *) malloc (sizeof(listz_iterator_t));
  if (iter == NULL)
    return NULL;
  iter->handle = h;

  if (iter->handle->storage == 0)
    iter->readptr = iter->writeptr = firstres;
  else
    {
      iter->offset = firstres * iter->handle->words;
      iter->readptr = iter->writeptr = iter->valid = 0;
      iter->dirty = 0;
      iter->bufsize = nr_buffered * iter->handle->words;
      iter->buf = malloc (iter->bufsize * sizeof(file_word_t));
      if (iter->buf == NULL) 
        {
          free (iter);
          return NULL;
        }
    }

  return iter;
}


listz_iterator_t *  
listz_iterator_init (listz_handle_t h, const size_t firstres)
{
  const size_t listz_iterator_nr_buffered = 65536;
  return listz_iterator_init2 (h, firstres, listz_iterator_nr_buffered);
}


void  
listz_iterator_clear (listz_iterator_t *iter)
{
  if (iter->handle->storage == 1)
    {
      listz_iterator_flush (iter);
      free (iter->buf);
    }
  free (iter);
}


void
listz_iterator_read (listz_iterator_t *iter, mpz_t r)
{
  if (iter->handle->storage == 0)
    {
      mpz_set (r, iter->handle->data.mem[iter->readptr++]);
    }
  else
    {
      if (iter->readptr == iter->valid)
        {
          listz_iterator_flush (iter);
          listz_iterator_fetch (iter, iter->offset + iter->valid);
          ASSERT_ALWAYS (iter->valid > 0);
          iter->readptr = iter->writeptr = 0;
        }
      mpz_import (r, iter->handle->words, -1, sizeof(file_word_t), 0, 0, 
                  &iter->buf[iter->readptr]);
      iter->readptr += iter->handle->words;
    }
}


void
listz_iterator_write (listz_iterator_t *iter, const mpz_t r)
{
  if (iter->handle->storage == 0)
    {
      mpz_set (iter->handle->data.mem[iter->writeptr++], r);
    }
  else
    {
      if (iter->writeptr == iter->bufsize)
        {
          listz_iterator_flush (iter);
          iter->writeptr = 0;
        }
      ASSERT_ALWAYS (mpz_sgn (r) >= 0);
      export_residue (&iter->buf[iter->writeptr], iter->handle->words, r);
      iter->writeptr += iter->handle->words;
      iter->dirty = 1;
    }
}


