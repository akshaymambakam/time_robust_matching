import robust_tre
from robust_tre import zone_set

import matplotlib.pyplot as plt
import numpy as np

# Function that returns the time robustness (translation) oracle from the match set
def gen_rob_oracle(zset: zone_set, start_scope: float, end_scope: float) -> callable:
	return lambda x,y: robust_tre.trobustness(zset,x,start_scope,end_scope)

# Function that returns the optimal robustness (translation) oracle from the match set
def gen_rob_oracle_opt(zset: zone_set, start_scope: float, end_scope: float) -> callable:
	return lambda x,y: robust_tre.trobustness_opt(zset,x,y,start_scope,end_scope)

# Initialize scaling factor
scaling = 1000

# Initialize heat map display scope
start_scope = 0.0
end_scope = 20.0

# Define the duration of signal
duration = end_scope - start_scope

# Create a match set of choice
zset = robust_tre.zone_set()
zset.add_from_period(start_scope, 10.0)
zset.add_from_period(10.0, 15.0)
zset.add_from_period(15.1,end_scope)
robust_tre.zsetprint(zset)

# Filter out points where robustness is less than rob_bound = 0.1
rob_bound = 0.1
zset_rob = robust_tre.trmtrans(zset, rob_bound)

# Print the robust zone set
robust_tre.zsetprint(zset_rob)

# Get the time robustness (translation) oracle
rob_oracle = gen_rob_oracle_opt(zset, start_scope, end_scope)

# generate 2 2d grids for the x & y bounds
y, x = np.meshgrid(np.linspace(start_scope, end_scope, scaling), np.linspace(start_scope, end_scope, scaling))

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
corners = (x.min(), x.max(), y.min(), y.max())
ax.axis(corners)
fig.colorbar(c, ax=ax)

plt.show()