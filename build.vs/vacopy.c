
#include <stdarg.h>
#include <string.h>

void _vacopy(va_list *pap, va_list ap)
{
    *pap = ap;
}
