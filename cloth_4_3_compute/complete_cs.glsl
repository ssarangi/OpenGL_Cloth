#version 430

#define DAMPING 0.01f // how much to damp the cloth simulation each frame
#define TIME_STEPSIZE2 0.5f*0.5f // how large time step each particle takes each frame

struct Particle
{
    vec4 position;
    vec4 accumulated_normal;
    vec2 uv;
    vec2 id;
    vec4 old_pos; // the position of the particle in the previous time step, used as part of the verlet numerical integration scheme
    vec4 acceleration; // a vector representing the current acceleration of the particle
    uint t1Idx;
    uint t2Idx;
    uint t3Idx;
    uint t4Idx;
    int movable; // can the particle move or not ? used to pin parts of the cloth
    float mass; // the mass of the particle (is always 1 in this example)
    int padding1;
    int padding2;
};

layout (std430, binding = 0) buffer ParticleBuffer
{
    Particle particles[];
} particleBuffer;

uint get_invocation()
{
   uint work_group = gl_WorkGroupID.x * gl_NumWorkGroups.y * gl_NumWorkGroups.z + gl_WorkGroupID.y * gl_NumWorkGroups.z + gl_WorkGroupID.z;
   return work_group * gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z + gl_LocalInvocationIndex;
}

void addForce(vec4 force_vec, unsigned int particleID)
{
    if (particleBuffer.particles[particleID].acceleration.y != 0.0)
        particleBuffer.particles[particleID].acceleration.w = particleID;
    else
    {
        vec4 f = vec4(0, -0.05, 0, 0);
        particleBuffer.particles[particleID].acceleration += f;
    }
}

// -------------------------------------------------- Gravity -------------------------------------------------- //
void addGravity(vec4 gravity, unsigned int particleID)
{
    addForce(gravity, particleID);
}

// -------------------------------------------------- Wind Force -------------------------------------------------- //
vec3 calcTriangleNormal(uint idx1, uint idx2, uint idx3)
{
    Particle p1 = particleBuffer.particles[idx1];
    Particle p2 = particleBuffer.particles[idx2];
    Particle p3 = particleBuffer.particles[idx3];

    vec4 pos1 = p1.position;
    vec4 pos2 = p2.position;
    vec4 pos3 = p3.position;

    vec4 v1 = pos2 - pos1;
    vec4 v2 = pos3 - pos1;

    return cross(vec3(v1), vec3(v2));
}

void addWindForcesForTriangle(vec4 direction, uint idx1, uint idx2, uint idx3)
{
    vec3 normal = calcTriangleNormal(idx1, idx2, idx3);
    vec3 d = normalize(normal);
    vec3 force = normal * (dot(d, vec3(direction)));

    barrier();
    addForce(vec4(force, 0), idx1);
    addForce(vec4(force, 0), idx2);
    addForce(vec4(force, 0), idx3);
}

void addWindForce(vec4 wind, unsigned int particleID)
{
    Particle p = particleBuffer.particles[particleID];

    addWindForcesForTriangle(wind, p.t1Idx, particleID, p.t2Idx);
    addWindForcesForTriangle(wind, p.t3Idx, particleID, p.t4Idx);
}

void satisfyConstraint(unsigned int particleID)
{

}

void vertlet(unsigned int particleID)
{
}

void ballCollision(unsigned int particleID)
{

}

layout (local_size_x = 4, local_size_y = 4, local_size_z = 1) in;

void main()
{
    uint flattened_id = get_invocation();

    particleBuffer.particles[flattened_id].id.y = flattened_id;

    // vec4 gravity = vec4(0, -0.2, 0, 0) * TIME_STEPSIZE2;
    vec4 gravity = vec4(0, -0.05, 0, 0);
    particleBuffer.particles[flattened_id].old_pos = gravity;
    addGravity(gravity, flattened_id);
    
    vec4 windForce = vec4(0.5, 0, 0.2, 0.0) * TIME_STEPSIZE2;
    //addWindForce(windForce, flattened_id);
    //vertlet(flattened_id);
    //ballCollision(flattened_id);
}