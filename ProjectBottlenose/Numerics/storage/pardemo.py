import numpy as np
from joblib import Parallel, delayed
import multiprocessing
# what are your inputs, and what operation do you want to
# perform on each input. For example...
def processInput(i):
	N = 10000000
	for j in range(1,N):
		y = np.sin(i/N)
	return y

if __name__ == '__main__':
    # do stuff with imports and functions defined about
	num_cores = multiprocessing.cpu_count()
	results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in range(num_cores*4))
