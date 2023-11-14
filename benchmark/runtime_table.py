import matplotlib.pyplot as plt
import numpy as np

plt.style.use("default")
seqlens = [1007,5055, 10122,20039,30007,400030,50120,60156,70107,80042,90003,100048,150043]


havac = [6.06,6.31,6.766,6.88,7.41,8.02,8.339,8.88,9.38,9.91,10.86,11.61,14.16]
nhmmer_32_threads = [2.36,8.32,20.53,49.75,70.72,101.33,130.17,151.93,177.86,209.37, 242.86,281.54, 434.84]

plt.xlabel("Model Database Length")
plt.ylabel("Time Taken (seconds)")
plt.yscale("log")
plt.plot(seqlens,nhmmer_32_threads, label = "nhmmer SSV (32 threads)")
plt.plot(seqlens, havac, label = "HAVAC")
plt.legend(loc='upper left')
plt.show()
