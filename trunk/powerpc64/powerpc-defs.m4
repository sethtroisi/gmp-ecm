divert(-1)

dnl  m4 macros for PowerPC assembler (32 and 64).
dnl  Inspired from GMP 4.1.4

dnl  Copyright 2000 Free Software Foundation, Inc.
dnl
dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or
dnl  modify it under the terms of the GNU Lesser General Public License as
dnl  published by the Free Software Foundation; either version 2.1 of the
dnl  License, or (at your option) any later version.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful,
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl  Lesser General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public
dnl  License along with the GNU MP Library; see the file COPYING.LIB.  If
dnl  not, write to the Free Software Foundation, Inc., 59 Temple Place -
dnl  Suite 330, Boston, MA 02111-1307, USA.


dnl  Usage: r0 ... r31, cr0 ... cr7
dnl
dnl  Registers names, either left as "r0" etc or mapped to plain 0 etc,
dnl  according to the result of GMP_ASM_POWERPC_REGISTERS.

define(r0,0)
define(r1,1)
define(r3,3)
define(r4,4)
define(r5,5)
define(r6,6)
define(r7,7)
define(r8,8)
define(r9,9)
define(r10,10)
define(r11,11)
define(r12,12)
define(r13,13)
define(r14,14)
define(r15,15)
define(r16,16)

divert
