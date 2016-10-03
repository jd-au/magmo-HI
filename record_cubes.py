

from __future__ import print_function

import os
import sys

# Hack to get parent folder in path
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import magmo

def main():
    #for day in range(11, 30+1):
    for day in range(21, 21+1):
        sources = magmo.get_day_obs_data(day)
        print (sources)
        for src in sources:
            cube_path = "day{}/1420/magmo-{}_1420_sl_restor.fits" \
                .format(day, src['source'])
            cube_found = os.path.exists(cube_path)
            print("{},{},{}".format(day, src['source'],
                                    "Y" if cube_found else "N"))
    return 0


# Run the script if it is called from the command line
if __name__ == "__main__":
    exit(main())

