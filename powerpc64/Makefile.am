MULREDC = mulredc1.asm mulredc2.asm mulredc3.asm mulredc4.asm mulredc5.asm \
          mulredc6.asm mulredc7.asm mulredc8.asm mulredc9.asm mulredc10.asm \
          mulredc11.asm mulredc12.asm mulredc13.asm mulredc14.asm \
          mulredc15.asm mulredc16.asm mulredc17.asm mulredc18.asm \
          mulredc19.asm mulredc20.asm

EXTRA_DIST = Makefile.dev README mulredc_1_2.m4 mulredc.m4 generate_all \
             powerpc-defs.m4

noinst_LTLIBRARIES = libmulredc.la
noinst_HEADERS = mulredc.h

# This library definition also causes the mulredc[n].asm and redc.asm files 
# to go in the distribution - no need for having them in EXTRA_DIST
libmulredc_la_SOURCES = $(MULREDC) redc.asm

.asm.s:
	$(M4) -I../ -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.s
.asm.S:
	$(M4) -I../ -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.S
