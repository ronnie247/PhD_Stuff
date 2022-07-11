from math import degrees
import numpy as np
from scipy import spatial
from stl import mesh  #need to use python instead of python3 to run this code
from matplotlib.pyplot import plot_mesh

# Import list of vertices from PDB file
vertices = np.array([
[-3, -3, 0],
[+3, -3, 0],
[+3, +3, 0],
[-3, +3, 0],
[+0, +0, +3]])

hull = spatial.ConvexHull(vertices)
faces = hull.simplices

myramid_mesh = mesh.Mesh(
  np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype)
)
for i, f in enumerate(faces):
  for j in range(3):
    myramid_mesh.vectors[i][j] = vertices[f[j],:]
    plot_mesh(myramid_mesh)

myramid_mesh.save('numpy_stl_example_02.stl')

# See ifpoint is inside mesh object:
# Tolerance is in degrees
# from math import pi, acos
# def is_inside(target_pt_global, mesh_obj, tolerance=0.02):
#     # Convert the point from global space to mesh local space
#     target_pt_local = mesh_obj.matrix_world.inverted() * target_pt_global
#     # Find the nearest point on the mesh and the nearest face normal
#     _, pt_closest, face_normal, _ = mesh_obj.closest_point_on_mesh(target_pt_local)
#     # Get the target-closest pt vector
#     target_closest_pt_vec = (pt_closest - target_pt_local).normalized()
#     # Compute the dot product = |a||b|*cos(angle)
#     dot_prod = target_closest_pt_vec.dot(face_normal)
#     # Get the angle between the normal and the target-closest-pt vector (from the dot prod)
#     angle = acos(min(max(dot_prod, -1), 1)) * 180 / pi
#     # Allow for some rounding error
#     inside = angle < 90-tolerance
#     return inside