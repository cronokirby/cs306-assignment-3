# 3.1.1

We just followed the example math given, and calculated the normal,
by taking the cross product of the two edge vectors, taken in the right
order. To do this, we did (v2 - v1) x (v3 - v1). For the angle weights,
we weren't sure how to do this at first, but then we realized that
since the calculation using them normalizes the result anyways, it doesn't matter
as long as its proportional to the angle. So, we just used the angle
(there was already a function to calculate the angle between two vectors).

# 3.1.2

We followed the recommendation, where for each triangle, we go over each
of the vertices of that triangle, and add in a contribution for that vertex
from the triangle, according to the angle weight of that vertex in this triangle,
and the normal calculated in the previous step.

Then we go over all of the vertices, and normalize what we've accumulated so far.

# 3.2.1

We pulled out the pen and paper in order to calculate how to solve
the system. In order to use Cramer's rule, we need to add our own determinant
function. To do this, we made use of the fact that taking the transpose
doesn't change the determinant, so we did dot(col1, cross(col2, col3)).
Then, we wrote a function to solve a 3x3 system using cramer's rule.
At first we simply ignored gamma, solving a system with columns
-d A B, but then we realized that we needed to include the fact
that gamma was (1 - alpha - beta), giving us a system
-d (A - C) (B - C), with an offset of O - C.

Solving this gave us a value of alpha and beta, as well as t. We then
calculate gamma as 1 - alpha - beta, and check that all three of these
are contained within 0 and 1.

# 3.2.2

For phong shading, this was a pretty straightforward extension, interpolating
the normals using alpha, beta, and gamma.

One issue that we had, and never managed to solve, was that in one of the scenes,
there was jaggedness along what was a shadow in the "unmeshed" version. We
looked at just the coloring of the normals, and got the exact same results
between both the meshed and unmeshed sphere, so we knew that our phong calculation
was correct. Even going over this with the TAs, we weren't able to figure out why this
was happening, so oh well.

# 3.3

We were a bit stumped at how to calculate the intersection with the box at first,
and then we realized that a box consisted of three intervals, one for each coordinate,
and it suffices to figure out if the ray ends up in those intervals coordinate
wise.

This gives a simple equation:

t_x in [(a_x - o_x) / d_x, (b_x - o_x) / d_x]

And we can take the intersection of each of the coordinate intervals, the get the final
one.

With a final interval, we realized that we just had to check that the lower bound wasn't
greater than the upper bound, and that the upper bound was greater than 0, to prevent
intersections behind the camera.

One issue we had was that we forgot to make sure that the newly scaled interval
was oriented properly, by making sure that the calculated lower
bound wasn't actually flipped. This happens when the ray direction is negative in
some coordinate. In this case, our interval was "backwards", leading us to
see only a section of the sphere mesh.

Once we realized this, everything worked as expected.

## Workload

Lúcás 33%
Mohamed 33%
Noah 33%
