import * as vec3 from "../lib/gl-matrix_3.3.0/esm/vec3.js"

function get_vert(mesh, vert_id) {

    const offset = vert_id*3
    const scaled_vertex = mesh.tris.vertices.slice(offset, offset+3).map(function(x) {return mesh.scale*x;})
    return vec3.add([0., 0., 0.], scaled_vertex, mesh.offset)
}

export function compute_triangle_normals_and_angle_weights(mesh) {

    /** TODO 3.1.1: 
	- compute the normal vector to each triangle in the mesh
    - push it into the array `tri_normals`
    - compute the angle weights for vert1, vert2, then vert3 and store it into an array [w1, w2, w3]
    - push this array into `angle_weights`

    Hint: you can use `vec3` specific methods such as `normalize()`, `add()`, `cross()`, `angle()`, or `subtract()`.
          The absolute value of a float is given by `Math.abs()`.
	*/

	const num_faces     = (mesh.tris.indices.length / 3) | 0
    const tri_normals   = []
    const angle_weights = []
    for(let i_face = 0; i_face < num_faces; i_face++) {
        const vert1 = get_vert(mesh, mesh.tris.indices[3*i_face + 0])
        const vert2 = get_vert(mesh, mesh.tris.indices[3*i_face + 1])
        const vert3 = get_vert(mesh, mesh.tris.indices[3*i_face + 2])

        let normal = vec3.cross([0, 0, 0], vec3.sub([0, 0, 0], vert2, vert1), vec3.sub([0, 0, 0], vert3, vert1))
        normal = vec3.normalize([0, 0, 0], normal)
        
        const angle1 = Math.abs(vec3.angle(vec3.sub([0, 0, 0], vert3, vert1), vec3.sub([0, 0, 0], vert2, vert1)))
        const angle2 = Math.abs(vec3.angle(vec3.sub([0, 0, 0], vert3, vert2), vec3.sub([0, 0, 0], vert1, vert2)))
        const angle3 = Math.abs(vec3.angle(vec3.sub([0, 0, 0], vert2, vert3), vec3.sub([0, 0, 0], vert1, vert3)))
        // Modify the way triangle normals and angle_weights are computed
        tri_normals.push(normal)
        angle_weights.push([angle1, angle2, angle3])
    }
    return [tri_normals, angle_weights]
}

export function compute_vertex_normals(mesh, tri_normals, angle_weights) {

    /** TODO 3.1.2: 
	- go through the triangles in the mesh
    - add the contribution of the current triangle to its vertices' normal
    - normalize the obtained vertex normals
	*/

	const num_faces    = (mesh.tris.indices.length / 3) | 0
    const num_vertices = (mesh.tris.vertices.length / 3) | 0
    let vertex_normals = Array(num_vertices).fill([0., 0., 0.])

    for(let i_face = 0; i_face < num_faces; i_face++) {
        const iv1 = mesh.tris.indices[3*i_face + 0]
        const iv2 = mesh.tris.indices[3*i_face + 1]
        const iv3 = mesh.tris.indices[3*i_face + 2]

        const scratch = [0, 0, 0];
        const [w1, w2, w3] = angle_weights[i_face]
        vec3.scale(scratch, tri_normals[i_face], w1)
        vec3.add(vertex_normals[iv1], scratch, vertex_normals[iv1])
        vec3.scale(scratch, tri_normals[i_face], w2)
        vec3.add(vertex_normals[iv2], scratch, vertex_normals[iv2])
        vec3.scale(scratch, tri_normals[i_face], w3)
        vec3.add(vertex_normals[iv3], scratch, vertex_normals[iv3])
    }

    for(let i_vertex = 0; i_vertex < num_vertices; i_vertex++) {
        // Normalize the vertex normals
        vec3.normalize(vertex_normals[i_vertex], vertex_normals[i_vertex]);
    }

    return vertex_normals
}