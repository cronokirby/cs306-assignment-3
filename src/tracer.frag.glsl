precision highp float;

#define MAX_RANGE 1e6
//#define NUM_REFLECTIONS
//#define F_VISUALIZE_AABB

//#define NUM_SPHERES
#if NUM_SPHERES != 0
uniform vec4 spheres_center_radius[NUM_SPHERES]; // ...[i] = [center_x, center_y, center_z, radius]
#endif

//#define NUM_PLANES
#if NUM_PLANES != 0
uniform vec4 planes_normal_offset[NUM_PLANES]; // ...[i] = [nx, ny, nz, d] such that dot(vec3(nx, ny, nz), point_on_plane) = d
#endif

//#define NUM_CYLINDERS
struct Cylinder {
	vec3 center;
	vec3 axis;
	float radius;
	float height;
};
#if NUM_CYLINDERS != 0
uniform Cylinder cylinders[NUM_CYLINDERS];
#endif


//#define NUM_TRIANGLES
//#define NUM_MESHES
//#define FLAT_SHADING_STRATEGY
//#define PHONG_SHADING_STRATEGY

#if defined PHONG_SHADING_STRATEGY
struct Triangle {
	mat3 vertices;
 	mat3 vertex_normals;
};
#else
struct Triangle {
	mat3 vertices;
};
#endif

struct AABB {
	vec3 corner_min;
	vec3 corner_max;
};
#if NUM_TRIANGLES != 0
//uniform Triangle triangles[NUM_TRIANGLES];
uniform AABB mesh_extent;
uniform sampler2D mesh_vert_pos_normals;
uniform int mesh_material_id;
#endif

// materials
//#define NUM_MATERIALS
struct Material {
	vec3 color;
	float ambient;
	float diffuse;
	float specular;
	float shininess;
	float mirror;
};
uniform Material materials[NUM_MATERIALS];
#if (NUM_SPHERES != 0) || (NUM_PLANES != 0) || (NUM_CYLINDERS != 0) || (NUM_TRIANGLES != 0)
uniform int object_material_id[NUM_SPHERES+NUM_PLANES+NUM_CYLINDERS+NUM_MESHES];
#endif

// lights
//#define NUM_LIGHTS
struct Light {
	vec3 color;
	vec3 position;
};
#if NUM_LIGHTS != 0
uniform Light lights[NUM_LIGHTS];
#endif
uniform vec3 light_color_ambient;


varying vec3 v2f_ray_origin;
varying vec3 v2f_ray_direction;

/*
	Solve the quadratic a*x^2 + b*x + c = 0. The method returns the number of solutions and store them
	in the argument solutions.
*/
int solve_quadratic(float a, float b, float c, out vec2 solutions) {

	// Linear case: bx+c = 0
	if (abs(a) < 1e-12) {
		if (abs(b) < 1e-12) {
			// no solutions
			return 0; 
		} else {
			// 1 solution: -c/b
			solutions[0] = - c / b;
			return 1;
		}
	} else {
		float delta = b * b - 4. * a * c;

		if (delta < 0.) {
			// no solutions in real numbers, sqrt(delta) produces an imaginary value
			return 0;
		} 

		// Avoid cancellation:
		// One solution doesn't suffer cancellation:
		//      a * x1 = 1 / 2 [-b - bSign * sqrt(b^2 - 4ac)]
		// "x2" can be found from the fact:
		//      a * x1 * x2 = c

		// We do not use the sign function, because it returns 0
		// float a_x1 = -0.5 * (b + sqrt(delta) * sign(b));
		float sqd = sqrt(delta);
		if (b < 0.) {
			sqd = -sqd;
		}
		float a_x1 = -0.5 * (b + sqd);


		solutions[0] = a_x1 / a;
		solutions[1] = c / a_x1;

		// 2 solutions
		return 2;
	} 
}

/*
	Check for intersection of the ray with a given sphere in the scene.
*/
bool ray_sphere_intersection(
		vec3 ray_origin, vec3 ray_direction, 
		vec3 sphere_center, float sphere_radius, 
		out float t, out vec3 normal) 
{
	vec3 oc = ray_origin - sphere_center;

	vec2 solutions; // solutions will be stored here

	int num_solutions = solve_quadratic(
		// A: t^2 * ||d||^2 = dot(ray_direction, ray_direction) but ray_direction is normalized
		1., 
		// B: t * (2d dot (o - c))
		2. * dot(ray_direction, oc),	
		// C: ||o-c||^2 - r^2				
		dot(oc, oc) - sphere_radius*sphere_radius,
		// where to store solutions
		solutions
	);

	// result = distance to collision
	// MAX_RANGE means there is no collision found
	t = MAX_RANGE+10.;
	bool collision_happened = false;

	if (num_solutions >= 1 && solutions[0] > 0.) {
		t = solutions[0];
	}
	
	if (num_solutions >= 2 && solutions[1] > 0. && solutions[1] < t) {
		t = solutions[1];
	}

	if (t < MAX_RANGE) {
		vec3 intersection_point = ray_origin + ray_direction * t;
		normal = (intersection_point - sphere_center) / sphere_radius;

		return true;
	} else {
		return false;
	}	
}

/*
	Check for intersection of the ray with a given plane in the scene.
*/
bool ray_plane_intersection(
		vec3 ray_origin, vec3 ray_direction, 
		vec3 plane_normal, float plane_offset, 
		out float t, out vec3 normal) 
{
		// can use the plane center if you need it
	vec3 plane_center = plane_normal * plane_offset;
	// If direction is orthogonal to normal, we're parallel to plane
	t = MAX_RANGE + 10.0;
	if (abs(dot(ray_direction, plane_normal)) < 1e-12) {
		return false;
	}
	t = dot(plane_normal, plane_center - ray_origin) / dot(plane_normal, ray_direction);
	normal = plane_normal;
	if (dot(ray_direction, normal) >= 0.0) {
		normal = -normal;
	}
	return t > 0.0;
}

/*
	Check for intersection of the ray with a given cylinder in the scene.
*/
bool ray_cylinder_intersection(
		vec3 ray_origin, vec3 ray_direction, 
		Cylinder cyl,
		out float t, out vec3 normal) 
{
	vec3 co = ray_origin - cyl.center;
	vec3 a = cyl.axis / dot(cyl.axis, cyl.axis);
	vec3 u = ray_direction - dot(a, ray_direction) * a;
	vec3 v = co - dot(a, co) * a;
	vec2 solutions;
	int num_solutions = solve_quadratic(
		dot(u, u),
		2.0 * dot(u, v),
		dot(v, v) - cyl.radius * cyl.radius,
		solutions
	);

	t = MAX_RANGE+10.;

	if (num_solutions >= 1 && solutions[0] > 0.) {
		t = solutions[0];
	}
	
	if (num_solutions >= 2 && solutions[1] > 0. && solutions[1] < t) {
		vec3 intersection_point = ray_origin + ray_direction * solutions[1];
		vec3 cyl_point = intersection_point - cyl.center;
		float z = dot(cyl_point, a);
		if (abs(z) <= cyl.height / 2.0) {
			t = solutions[1];
			normal = (cyl_point - z * a) / cyl.radius;
			if (dot(ray_direction, normal) >= 0.0) {
				normal = -normal;
			}
			return true;
		}
	}

	if (t < MAX_RANGE) {
		vec3 intersection_point = ray_origin + ray_direction * t;
		vec3 cyl_point = intersection_point - cyl.center;
		float z = dot(cyl_point, a);
		if (abs(z) > cyl.height / 2.0) {
			t = MAX_RANGE + 10.0;
			return false;
		}
		normal = (cyl_point - z * a) / cyl.radius;
		if (dot(ray_direction, normal) >= 0.0) {
			normal = -normal;
		}
		return true;
	} else {
		return false;
	}
}


bool ray_AABB_filter(
	vec3 ray_origin, vec3 ray_direction, AABB aabb)
{
	/** TODO 3.3: 
	- construct the range of the ray parameter (`t`) values which are included in the bounding box
	- check that this range is non-empty
	- return whether the bounding box is intersected by the ray or not
	*/

	return true;
}


#if NUM_TRIANGLES != 0
Triangle get_triangle(int idx) {
	Triangle tri;

	// Texture is sampled in 0..1 coordinates
	float idx_norm = float(idx) / (float(NUM_TRIANGLES)-1.);

	tri.vertices[0] = texture2D(mesh_vert_pos_normals, vec2(0./5., idx_norm)).xyz;
	tri.vertices[1] = texture2D(mesh_vert_pos_normals, vec2(1./5., idx_norm)).xyz;
	tri.vertices[2] = texture2D(mesh_vert_pos_normals, vec2(2./5., idx_norm)).xyz;

	#ifdef PHONG_SHADING_STRATEGY
	tri.vertex_normals[0] = texture2D(mesh_vert_pos_normals, vec2(3./5., idx_norm)).xyz;
	tri.vertex_normals[1] = texture2D(mesh_vert_pos_normals, vec2(4./5., idx_norm)).xyz;
	tri.vertex_normals[2] = texture2D(mesh_vert_pos_normals, vec2(5./5., idx_norm)).xyz;
	#endif
	return tri;
}
#endif

float det3(mat3 A) {
	return dot(A[0], (cross(A[1], A[2])));
}

vec3 solve(mat3 A, vec3 b) {
	float under = det3(A);

	vec3 ret;

	mat3 scratch = A;
	scratch[0] = b;	
	ret.x = det3(scratch) / under;

	scratch = A;
	scratch[1] = b;	
	ret.y = det3(scratch) / under;

	scratch = A;
	scratch[2] = b;	
	ret.z = det3(scratch) / under;

	return ret;
}

bool ray_triangle_intersection(
		vec3 ray_origin, vec3 ray_direction, 
		Triangle tri,
		out float t, out vec3 normal) 
{	
	/** TODO 3.2.1: 
	- compute the ray's first valid intersection with the triangle
	- check that the intersection occurs in front of the viewer: `t > 0`
	- store ray parameter in `t`
	- store normal at intersection_point in `normal`.
	- return whether there is an intersection with `t > 0`
	*/

	/** TODO 3.2.2: 
	- test which of `FLAT_SHADING_STRATEGY` or `PHONG_SHADING_STRATEGY` is defined
	- implement flat shading using the triangle normal you computed in task 3.2.1
	- implement Phong shading using the vertex normals stored in `tri.vertex_normals`

	Hint: vertex normals are stored in the same order as the vertex positions
	*/
	vec3 intersection_point;
	t = MAX_RANGE + 10.;

	vec3 p0 = tri.vertices[0];
	vec3 p1 = tri.vertices[1];
	vec3 p2 = tri.vertices[2];
	mat3 A;
	A[0] = -ray_direction;
	A[1] = p0 - p2;
	A[2] = p1 - p2;
	vec3 solution = solve(A, ray_origin - p2);
	float alpha = solution.y;
	float beta = solution.z;
	float gamma = 1.0 - alpha - beta;;
	if (alpha < 0.0 || beta < 0.0 || gamma < 0.0 || alpha > 1.0 || beta > 1.0 || gamma > 1.0) {
		return false;
	}
	if (solution.x < 0.0) {
		return false;
	}
	t = solution.x;

	intersection_point = ray_origin + t * ray_direction;
	#if defined PHONG_SHADING_STRATEGY
		normal = alpha * tri.vertex_normals[0] + beta * tri.vertex_normals[1] + gamma * tri.vertex_normals[2];
	#else
		normal = normalize(cross(p1 - p0, p2 - p0));
	#endif
	if (dot(ray_direction, normal) >= 0.0) {
		normal = -normal;
	}
	return true;
}



/*
	Check for intersection of the ray with any object in the scene.
*/
bool ray_intersection(
		vec3 ray_origin, vec3 ray_direction, 
		out float col_distance, out vec3 col_normal, out int material_id) 
{
	col_distance = MAX_RANGE + 10.;
	col_normal = vec3(0., 0., 0.);

	float object_distance;
	vec3 object_normal;

	// Check for intersection with each sphere
	#if NUM_SPHERES != 0 // only run if there are spheres in the scene
	for(int i = 0; i < NUM_SPHERES; i++) {
		bool b_col = ray_sphere_intersection(
			ray_origin, 
			ray_direction, 
			spheres_center_radius[i].xyz, 
			spheres_center_radius[i][3], 
			object_distance, 
			object_normal
		);

		// choose this collision if its closer than the previous one
		if (b_col && object_distance < col_distance) {
			col_distance = object_distance;
			col_normal = object_normal;
			material_id =  object_material_id[i];
		}
	}
	#endif

	// Check for intersection with each plane
	#if NUM_PLANES != 0 // only run if there are planes in the scene
	for(int i = 0; i < NUM_PLANES; i++) {
		bool b_col = ray_plane_intersection(
			ray_origin, 
			ray_direction, 
			planes_normal_offset[i].xyz, 
			planes_normal_offset[i][3], 
			object_distance, 
			object_normal
		);

		// choose this collision if its closer than the previous one
		if (b_col && object_distance < col_distance) {
			col_distance = object_distance;
			col_normal = object_normal;
			material_id =  object_material_id[NUM_SPHERES+i];
		}
	}
	#endif

	// Check for intersection with each cylinder
	#if NUM_CYLINDERS != 0 // only run if there are cylinders in the scene
	for(int i = 0; i < NUM_CYLINDERS; i++) {
		bool b_col = ray_cylinder_intersection(
			ray_origin, 
			ray_direction,
			cylinders[i], 
			object_distance, 
			object_normal
		);

		// choose this collision if its closer than the previous one
		if (b_col && object_distance < col_distance) {
			col_distance = object_distance;
			col_normal = object_normal;
			material_id =  object_material_id[NUM_SPHERES+NUM_PLANES+i];
		}
	}
	#endif

	// Check for intersection with each triangle
	#if NUM_TRIANGLES != 0 // only run if there are triangles in the scene
	if( ray_AABB_filter(ray_origin, ray_direction, mesh_extent) ) {
	// if (true){
		for(int i = 0; i < NUM_TRIANGLES; i++) {
			bool b_col = ray_triangle_intersection(
				ray_origin, 
				ray_direction,
				get_triangle(i),
				object_distance, 
				object_normal
			);

			// choose this collision if its closer than the previous one
			if (b_col && object_distance < col_distance) {
				col_distance = object_distance;
				col_normal = object_normal;
				// there is only one mesh, and it has the last entry in the array
				material_id = object_material_id[NUM_SPHERES+NUM_PLANES+NUM_CYLINDERS];
				//material_id = mesh_material_id;
			}
		}
	}
	#endif

	return col_distance < MAX_RANGE;
}

/*
	Return the color at an intersection point given a light and a material, exluding the contribution
	of potential reflected rays.
*/
vec3 lighting(
		vec3 object_point, vec3 object_normal, vec3 direction_to_camera, 
		Light light, Material mat) {
  vec3 contribution = vec3(0.0);

	/** TODO 2.1: 
	- compute the diffuse component
	- make sure that the light is located in the correct side of the object
	- compute the specular component 
	- make sure that the reflected light shines towards the camera
	- return the ouput color
	*/
	vec3 l = normalize(light.position - object_point);
	float diffuse_amount = dot(object_normal, l);
	vec3 r = reflect(-l, object_normal);
	vec3 v = normalize(direction_to_camera);
	if (diffuse_amount < 0.0) {
		return vec3(0.0);
	}
	contribution += mat.color * mat.diffuse * diffuse_amount;
	if (dot(r, v) >= 0.0) {
		contribution += mat.color * mat.specular * pow(dot(r, v), mat.shininess);
	}


	/** TODO 2.2: 
	- shoot a shadow ray from the intersection point to the light
	- check whether it intersects an object from the scene
	- update the lighting accordingly
	*/
	float col_distance;
	vec3 col_normal = vec3(0.);
	int mat_id      = 0;
	ray_intersection(object_point + 1e-4 * object_normal, l, col_distance, col_normal, mat_id);
	if (col_distance < length(light.position - object_point)) {
		return vec3(0.0);
	}
	

	return contribution;
}

/*
	Get the material corresponding to mat_id from the list of materials.
*/
Material get_mat2(int mat_id) {
	Material m = materials[0];
	for(int mi = 1; mi < NUM_MATERIALS; mi++) {
		if(mi == mat_id) {
			m = materials[mi];
		}
	}
	return m;
}

void main() {
	vec3 ray_origin = v2f_ray_origin;
	vec3 ray_direction = normalize(v2f_ray_direction);

	vec3 pix_color = vec3(0.);
	float reflection_weight = 1.0;

	for (int i_reflection = 0; i_reflection < NUM_REFLECTIONS + 1; i_reflection++) {
		float col_distance;
		vec3 col_normal = vec3(0.);
		int mat_id      = 0;
		ray_intersection(ray_origin, ray_direction, col_distance, col_normal, mat_id);
		if (col_distance >= MAX_RANGE) {
			break;
		}

		Material mat = get_mat2(mat_id);
		vec3 color = mat.color * mat.ambient * light_color_ambient;

		vec3 object_point = col_distance * ray_direction + ray_origin;
		#if NUM_LIGHTS != 0
		for (int i = 0; i < NUM_LIGHTS; ++i) {
			color += lights[i].color * lighting(object_point, col_normal, ray_origin - object_point, lights[i], mat);
		}
		#endif

		pix_color += (1.0 - mat.mirror) * reflection_weight * color;
		reflection_weight *= mat.mirror;
		ray_direction = reflect(normalize(ray_direction), col_normal);
		ray_origin = object_point + 1e-4 * col_normal;
	}



	gl_FragColor = vec4(pix_color, 1.0);
}
