# Import NURBS-Python and the visualization component
from geomdl import BSpline
from geomdl import utilities
# from geomdl.visualization import VisPlotly as vis
from geomdl.visualization import VisMPL as vis
# Create a BSpline surface instance
surf = BSpline.Surface()
# Set evaluation delta
surf.delta = 0.05
# Set degrees
surf.degree_u = 3
surf.degree_v = 3
# Set control points
	
	
control_points = [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 0, 1.5], [1, 0, -0.5], [0, 1, 0], [0, 1, 0], [-1, 0, -0.5], [0, 0, 1.5], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]]
surf.set_ctrlpts(control_points, 4, 4)
# Auto generate knot vectors
surf.knotvector_u = utilities.generate_knot_vector(surf.degree_u, 4)
surf.knotvector_v = utilities.generate_knot_vector(surf.degree_v, 4)
# Evaluate surface
surf.evaluate()
# Visualization config
vis_config = vis.VisConfig(figure_size=[8,8])
# Visualization component
vis_comp = vis.VisSurface(vis_config)
# Set visualization component of the surface
surf.vis = vis_comp
# Render the surface
surf.render()