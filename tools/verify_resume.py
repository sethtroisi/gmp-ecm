#!/usr/bin/env python
#
# Script to verify first, middle, and last line of a resume file
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
#  $ python verify_resume.py <FN>
#
# To pass extra args not present in resume file
#  $ python verify_resume.py <FN> -- -k 3 -v


import argparse
import os

parser = argparse.ArgumentParser(description='Spot check a resume file')

parser.add_argument('resumefile', type=str,
    help='resume file to verify')


def get_lines(fn):
    if not os.path.exists(fn):
        raise ValueEroor(f"resume file {fn!r} does not exist")

    with open(fn) as f:
        lines = f.readlines()
        assert len(lines) > 0
        if len(lines) <= 3:
            return lines
        return (lines[0], lines[len(lines)//2], lines[-1])

    assert False

def parse_line(line):
    values = {}
    for pair in line.strip().split(";"):
        pair = pair.strip()
        if not pair:
            continue
        k,v = pair.strip().split("=")
        values[k] = v
    return values

def build_command(values, saveafn, extra_args):
    # TODO support different ecm command
    command = ["ecm"]

    N = values.pop("N")
    B1 = values.pop("B1")
    B2 = values.pop("B2", "0")

    IGNORED = ["X", "CHECKSUM", "PROGRAM", "WHO", "TIME"]
    for ignore in IGNORED:
        values.pop(ignore, None)

    for key, value in values.items():
        if key == "METHOD":
            METHOD_LOOKUP = {"P-1": "-pm1", "P+1": "-pp1"}
            if value != "ECM":
                 command.append(METHOD_LOOKUP[value])
        elif key == "SIGMA":
            command.append("-sigma")
            command.append(value)
        elif key == "A":
            command.append("-A")
            command.append(value)
        elif key == "X0":
            if value not in ("0", "0x0"):
                command.append("-x0")
                command.append(value)
        elif key == "Y":
            assert value in ("0", "0x0"), "Y not supported"
        elif key == "Y0":
            assert value in ("0", "0x0"), "Y0 not supported"
        else:
            print(f"{key!r} unhandled value: {value!r}")
    
    command.append("-savea")
    command.append(saveafn)
    command.extend(extra_args)
    command.append(B1)
    command.append(B2)

    return (N, command)

if __name__ == '__main__':
    args, extra = parser.parse_known_args()

    #print(args)
    #print(extra)

    VERIFY_TMP_FN = "resume.verify_resume.txt"

    lines = get_lines(args.resumefile)
    parsed = [parse_line(line) for line in lines]
    commands = [build_command(parse, VERIFY_TMP_FN, extra) for parse in parsed]

    # Add savea to commands
    print("Now run these commands:")
    for N, command in commands:
        print(f"echo '{N}' | " + " ".join(command))
    print()
    print("Then compare with")
    print(f"python compare_resume.py {args.resumefile!r} {VERIFY_TMP_FN!r}")
