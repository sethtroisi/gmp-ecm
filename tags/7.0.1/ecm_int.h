#ifndef ECM_INT_H_
#define ECM_INT_H_

#include "config.h"

#if defined UINT64_MAX || defined uint64_t
typedef int64_t ecm_int;
typedef uint64_t ecm_uint;
#define ECM_INT_MAX INT64_MAX
#define ECM_UINT_MAX UINT64_MAX
#elif defined HAVE_LONG_LONG_INT
typedef long long ecm_int;
typedef unsigned long long ecm_uint;
#define ECM_INT_MAX LLONG_MAX
#define ECM_UINT_MAX ULLONG_MAX
#else
typedef long ecm_int;
typedef unsigned long ecm_uint;
#define ECM_INT_MAX LONG_MAX
#define ECM_UINT_MAX ULONG_MAX
#endif

#endif /* ECM_INT_H_ */
