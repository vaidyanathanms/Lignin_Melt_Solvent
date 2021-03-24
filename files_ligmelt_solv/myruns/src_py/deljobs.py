# To delete jobs

import os
import sys
import math
import subprocess

deljob1 = 846558
deljob2 = 846572

for jobnum in range(deljob1,deljob2+1):
    subprocess.call(["bkill", str(jobnum)])
