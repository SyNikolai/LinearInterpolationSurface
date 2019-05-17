# LinearInterpolationSurface
This is a simple linear interpolation algorithm implemented in C++.
From a dataset of cartesian coordinates (x,y,z) that correspond to surface points of an object ( case: asteroid Didymos),
we determine the maximum distance between two neighboring points. 
Then, all adjacent distances (distances between a point and its neighbours) are calculated using index data provided in another file.
In this implimentation, everytime an adj_distance is found greater than the half of the maximum distance, a new point is added.
This program can be easily customized in terms of how and when a new point will be added in our data.

Important: The dataset of vertices is not interpolated, after running the program the new data will replace the original one.
To avoid that look at how the output function is defined and make any necessary changes.

