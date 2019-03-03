# Dubins' Curves implemented in Python

In 1957, Lester Dubins mathematically proved that the shorts distance between two waypoints (x, y, psi) can be calculated as a combination of straight lines and circular arcs, requiring only three segments.

# Running the Code

Pretty straight forward...pull or copy the python file, and download the necessary libraries (numpy, matplotlib, etc)
Then, at the bottom, you can configure your own waypoints. The waypoint structure is (x, y, psi), where Psi is a bearing in NED coordinates (Due north = 0 deg)


Here's an example output

![Dubins Curve](https://github.com/fgabbert/dubins_py/blob/master/dubins_example.png)
