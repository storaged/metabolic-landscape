#!/usr/bin/python

import sys

# sys.argv[1] - file to be extended

splitted = sys.argv[1].split(".")
res_filename = ''.join([splitted[0], "_ext.tab"])

print res_filename

with open(res_filename, "w") as res_file:
    with open(sys.argv[1]) as f:
        content = f.readlines()
        content[0] = "first_col\t" + content[0]
        for x in content:
            res_file.write(x)


