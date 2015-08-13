
#include "libntt.h"

int main(int argc, char **argv)
{
#if GMP_LIMB_BITS == 32
  test_main_sp30w32(argc, argv);
  test_main_sp30w32sse2(argc, argv);
  test_main_sp31w32(argc, argv);
  test_main_sp31w32sse2(argc, argv);
  test_main_sp62w32(argc, argv);
  test_main_sp62w32sse2(argc, argv);
#else
  test_main_sp30w64(argc, argv);
  test_main_sp30w64sse2(argc, argv);
  test_main_sp62w64(argc, argv);
  test_main_sp62w64sse2(argc, argv);
#endif
  return 0;
}
