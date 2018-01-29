# diskload

This is MATLAB and Python code to compute the elastic response to a disk load.

The `diskload` MATLAB function computes the response to a uniform surface pressure load imposed in a disc of a given angular radius and height. The elastic response is found at one or more points located on the surface of the Earth at specified angular distance(s) from the center of the disc load.

The elastic response is computed using user-supplied elastic loading Love numbers (h,k,l) generated using a specific elastic structure model for the Earth. 

The original MATLAB code was described in the following publication:

Bevis, M., Melini, D., Spada, G., 2016. On computing the geoelastic response to a disk load, Geophys. J. Int., 205 (3), 1,804-1,812, doi:10.1093/gji/ggw115.

## Python Example

```
import diskload

# Load the love numbers and extrapolate out farther
love = diskload.love_numbers.read()
love = diskload.love_numbers.extrapolate( love, 100000 )

# Compute the displacement (at 0.2 degrees) of a
# disk with angular radius 0.1 degrees

uncompensated = diskload.Compensation.UNCOMPENSATED
alpha = 0.1
theta = 0.2
w = 17.0
u, v, g = diskload.truncated( alpha, uncompensated, theta, w, 100000, love )
```
