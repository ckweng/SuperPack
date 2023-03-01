#!/usr/bin/env python3

import os

header = """\
/**
 * @file {filename}
 *
 * SCL --- Secure Computation Library
 * Copyright (C) 2022 Anders Dalskov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 * USA
 */
"""

def check_file(filename, path):
    expected_header = header.format(filename=filename).rstrip().split("\n")
    with open(path, 'r') as f:
        lines = f.readlines()
        for a, b in zip(expected_header, lines):
            if not a.rstrip() == b.rstrip():
                return False
    return True


all_good = True

for path, __, names in os.walk("include"):
    for n in names:
        if n.endswith(".h"):
            full_name = os.path.join(path, n)
            if not check_file(n, full_name):
                print(f"{full_name} missing copyright header")
                all_good = False


exit(0 if all_good else 1)
