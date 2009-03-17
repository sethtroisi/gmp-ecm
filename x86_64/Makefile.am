MULREDC = mulredc1.asm mulredc2.asm mulredc3.asm mulredc4.asm mulredc5.asm mulredc6.asm \
          mulredc7.asm mulredc8.asm mulredc9.asm mulredc10.asm mulredc11.asm \
          mulredc12.asm mulredc13.asm mulredc14.asm mulredc15.asm mulredc16.asm \
          mulredc17.asm mulredc18.asm mulredc19.asm mulredc20.asm
MULREDC1 = mulredc1_2.asm mulredc1_3.asm mulredc1_4.asm mulredc1_5.asm \
           mulredc1_6.asm mulredc1_7.asm mulredc1_8.asm mulredc1_9.asm mulredc1_10.asm \
           mulredc1_11.asm mulredc1_12.asm mulredc1_13.asm mulredc1_14.asm \
           mulredc1_15.asm mulredc1_16.asm mulredc1_17.asm mulredc1_18.asm \
           mulredc1_19.asm mulredc1_20.asm
EXTRA_DIST = autogen.py generate_all mulredc.m4 mulredc1.m4 mulredc.h

noinst_LTLIBRARIES = libmulredc.la

# This library definition also causes the mulredc[n].asm, mulredc1_[n].asm, 
# and redc.asm files to go in the distribution - no need for having 
# them in EXTRA_DIST
libmulredc_la_SOURCES = $(MULREDC) $(MULREDC1) redc.asm

# There has to be a way of making this automatically
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

mulredc1_2.asm: mulredc1.m4
	$(M4) -DLENGTH=2 $< > $@
mulredc1_3.asm: mulredc1.m4
	$(M4) -DLENGTH=3 $< > $@
mulredc1_4.asm: mulredc1.m4
	$(M4) -DLENGTH=4 $< > $@
mulredc1_5.asm: mulredc1.m4
	$(M4) -DLENGTH=5 $< > $@
mulredc1_6.asm: mulredc1.m4
	$(M4) -DLENGTH=6 $< > $@
mulredc1_7.asm: mulredc1.m4
	$(M4) -DLENGTH=7 $< > $@
mulredc1_8.asm: mulredc1.m4
	$(M4) -DLENGTH=8 $< > $@
mulredc1_9.asm: mulredc1.m4
	$(M4) -DLENGTH=9 $< > $@
mulredc1_10.asm: mulredc1.m4
	$(M4) -DLENGTH=10 $< > $@
mulredc1_11.asm: mulredc1.m4
	$(M4) -DLENGTH=11 $< > $@
mulredc1_12.asm: mulredc1.m4
	$(M4) -DLENGTH=12 $< > $@
mulredc1_13.asm: mulredc1.m4
	$(M4) -DLENGTH=13 $< > $@
mulredc1_14.asm: mulredc1.m4
	$(M4) -DLENGTH=14 $< > $@
mulredc1_15.asm: mulredc1.m4
	$(M4) -DLENGTH=15 $< > $@
mulredc1_16.asm: mulredc1.m4
	$(M4) -DLENGTH=16 $< > $@
mulredc1_17.asm: mulredc1.m4
	$(M4) -DLENGTH=17 $< > $@
mulredc1_18.asm: mulredc1.m4
	$(M4) -DLENGTH=18 $< > $@
mulredc1_19.asm: mulredc1.m4
	$(M4) -DLENGTH=19 $< > $@
mulredc1_20.asm: mulredc1.m4
	$(M4) -DLENGTH=20 $< > $@

.asm.s:
	$(M4) -I../ -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.s
.asm.S:
	$(M4) -I../ -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.S
