import multiprocessing as mp
import random
import string

# Define an output queue
output = mp.Queue()

# define a example function
def rand_string(length, pos, output):
    """ Generates a random string of numbers, lower- and uppercase chars. """
    rand_str = ''.join( random.choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for i in range(length) )
    output.put( (pos, rand_str) )

# Setup a list of processes that we want to run
processes = [mp.Process(target=rand_string, args=(5, x, output)) for x in range(4)]

# Run processes
for p in processes:
    p.start()

# Exit the completed processes
for p in processes:
    p.join()

# Get process results from the output queue
results = [output.get() for p in processes]

print(results)

xs = [results[i][1] for i in range(len(results))]

xs


''.join( random.choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for i in range(8) )

random.seed(1)


import random
print ("Random number with seed 30")
random.seed( 10 )
print ("first - ", random.randint(25,50))
#will generate a same random number as previous
random.seed( 10 )
print ("Second - ", random.randint(25,50))
#will generate a same random number as previous
random.seed( 10 )
print ("Third - ", random.randint(25,50))

random()


from matplotlib.pylab import *  # for random()
random()

rr
##
