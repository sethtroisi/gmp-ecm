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

# It's actually the .s files that depend on config.m4, but automake
# knows them only as intermediate files, not as targets. Adding the
# dependency to libmulredc.la should work so long as no stale .s
# files exist.
libmulredc_la_DEPENDENCIES = $(top_builddir)/config.m4
  
# Rules for generating the .asm files from the .m4 scripts
mulredc1.asm: mulredc_1_2.m4
	$(M4) -DLENGTH=1 $< > $@
mulredc2.asm: mulredc_1_2.m4
	$(M4) -DLENGTH=2 $< > $@
mulredc3.asm: mulredc.m4
	$(M4) -DLENGTH=3 $< > $@
mulredc4.asm: mulredc.m4
	$(M4) -DLENGTH=4 $< > $@
mulredc5.asm: mulredc.m4
	$(M4) -DLENGTH=5 $< > $@
mulredc6.asm: mulredc.m4
	$(M4) -DLENGTH=6 $< > $@
mulredc7.asm: mulredc.m4
	$(M4) -DLENGTH=7 $< > $@
mulredc8.asm: mulredc.m4
	$(M4) -DLENGTH=8 $< > $@
mulredc9.asm: mulredc.m4
	$(M4) -DLENGTH=9 $< > $@
mulredc10.asm: mulredc.m4
	$(M4) -DLENGTH=10 $< > $@
mulredc11.asm: mulredc.m4
	$(M4) -DLENGTH=11 $< > $@
mulredc12.asm: mulredc.m4
	$(M4) -DLENGTH=12 $< > $@
mulredc13.asm: mulredc.m4
	$(M4) -DLENGTH=13 $< > $@
mulredc14.asm: mulredc.m4
	$(M4) -DLENGTH=14 $< > $@
mulredc15.asm: mulredc.m4
	$(M4) -DLENGTH=15 $< > $@
mulredc16.asm: mulredc.m4
	$(M4) -DLENGTH=16 $< > $@
mulredc17.asm: mulredc.m4
	$(M4) -DLENGTH=17 $< > $@
mulredc18.asm: mulredc.m4
	$(M4) -DLENGTH=18 $< > $@
mulredc19.asm: mulredc.m4
	$(M4) -DLENGTH=19 $< > $@
mulredc20.asm: mulredc.m4
	$(M4) -DLENGTH=20 $< > $@

.asm.s:
	$(M4) -I../ -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.s
#.asm.S:
#	$(M4) -I../ -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.S
