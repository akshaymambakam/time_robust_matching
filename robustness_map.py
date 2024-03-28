import robust_tre

import matplotlib.pyplot as plt
import numpy as np

# Initialize scaling factor
scaling = 100
# Define the duration of signal
duration = 20

# Function that returns the time robustness (translation) oracle from the match set
def gen_rob_oracle(zset):
	return lambda x,y: robust_tre.trobustness(zset,x,y)

# Create a match set of choice
zset = robust_tre.zone_set()
zset.add_from_period(0,10)

# Filter out points where robustness is less than rob_bound = 0.1
rob_bound = 0.1
zset_rob = robust_tre.trmtrans(zset, rob_bound)

# Print the robust zone set
robust_tre.zsetprint(zset_rob)

# Get the time robustness (translation) oracle
rob_oracle = gen_rob_oracle(zset)

# generate 2 2d grids for the x & y bounds
y, x = np.meshgrid(np.linspace(0, duration, scaling), np.linspace(0, duration, scaling))

# Map the time robustness with respect to translation
z = np.zeros((scaling,scaling))
for i in range(scaling):
	for j in range(scaling):
		z[i,j] = rob_oracle(x[i,j],y[i,j])

# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
z_min, z_max = -np.abs(z).max(), np.abs(z).max()

fig, ax = plt.subplots()

c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
ax.set_title('Time Robustness (translation) map')
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax)

plt.show()