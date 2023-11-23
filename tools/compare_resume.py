#!/usr/bin/env python
#
# Script to compare two resume files
#
# Copyright 2023
# Seth Troisi
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, see
# http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
#
# Verify line from <FN>
#  $ python find_small_sigma.sage <FN>
#
# To pass extra args not present in resume file
#  $ python find_small_sigma.sage <FN> -- -k 3 -v


import argparse
import os

parser = argparse.ArgumentParser(description='Compare two resume files')

parser.add_argument('resumefile', type=str,
    help='1st resume file to verify')

parser.add_argument('verifyfile', type=str,
    help='2nd resume file to verify')


def get_lines(fn):
    if not os.path.exists(fn):
        raise ValueEroor(f"resume file {fn!r} does not exist")

    with open(fn) as f:
        return f.readlines()


def parse_line(line):
    values = {}
    for pair in line.strip().split(";"):
        pair = pair.strip()
        if not pair:
            continue
        k,v = pair.strip().split("=")
        values[k] = v
    return values


def value_to_keys(values):
    """Group values into (inputs, outputs). Drop some values (WHO, TIME)."""

    input_names = ("METHOD", "B1", "B2", "X0", "Y0", "A", "SIGMA", "N")
    output_names = ("X", "CHECKSUM")
    drop_names = ("PROGRAM", "WHO", "TIME")

    # Drop keys that have value = "0x0", these seem to just gum up the works
    for key in ("X0", "Y0", "Y"):
        if values.get(key, None) == "0x0":
            values.pop(key)

    # Possible don't include X0/Y0 if zero
    inputs = tuple((name, values.pop(name)) for name in input_names if name in values)
    outputs = tuple((name, values.pop(name)) for name in output_names if name in values)

    for name in drop_names:
        if name in values:
            values.pop(name)

    for key, value in values.items():
        print(f"{key!r} unhandled value: {value!r}")
    
    return inputs, outputs


def fn_to_grouped(fn):
    grouped = {}
    for line in get_lines(fn):
        values = parse_line(line)
        inputs, outputs = value_to_keys(values)
        assert inputs not in grouped, "Two lines had same inputs: " + str(inputs)
        grouped[inputs] = outputs

    return grouped


def shorten(value):
    if value.isnumeric():
        if len(value) > 20:
            return "{}...{}<{}>".format(value[:3], value[-3:], len(value))
        return value
    return value


def compare(args):
    a = fn_to_grouped(args.resumefile)
    b = fn_to_grouped(args.verifyfile)
    verified = 0

    for in_b, out_b in b.items():
        out_a = a.get(in_b, None)
        if out_a is None:
            print("Couldn't find input in resumefile:")
            print(in_b)
            print()
            for in_a in a:
                if in_b[-1] == in_a[-1]:
                    print("Here")
                    print(in_a)
            continue

        if out_a == out_b:
            verified += 1
            str_in = " ".join(f"{name}={shorten(value)}" for name, value in in_b if value != '0x0')
            str_out = " ".join(f"{name}={shorten(value)}" for name, value in out_b)
            print(f"Verified {str_in}  matches  {str_out}")
        else:
            print("Output mismatch for")
            print(in_b)
            print(out_a)
            print("vs")
            print(out_b)
            print()

    # TODO verify same number of lines
    print("Verified {} of {} lines from verifyfile".format(
        verified, len(b), len(a)))

if __name__ == '__main__':
    args = parser.parse_args()

    compare(args)

