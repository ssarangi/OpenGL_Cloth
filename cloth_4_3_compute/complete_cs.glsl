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
    int movable; // can the particle move or not ? used to pin parts of the cloth
    float mass; // the mass of the particle (is always 1 in this example)
    int padding1;
    int padding2;
};

layout (std430, binding = 0) buffer ParticleBuffer
{
    Particle particles[];
} particleBuffer;

layout (std430, binding = 1) buffer IDBuffer
{
    unsigned int ids[];
} idBuffer;

void addForce(vec4 force_vec, unsigned int particleID)
{
    particleBuffer.particles[particleID].acceleration += force_vec / particleBuffer.particles[particleID].mass;
}

void addGravity(vec4 gravity, unsigned int particleID)
{
    addForce(gravity, particleID);
}

void addWindForce(vec4 wind, unsigned int particleID)
{

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

uint get_invocation()
{
   uint work_group = gl_WorkGroupID.x * gl_NumWorkGroups.y * gl_NumWorkGroups.z + gl_WorkGroupID.y * gl_NumWorkGroups.z + gl_WorkGroupID.z;
   return work_group * gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z + gl_LocalInvocationIndex;
}

void main()
{
    uint flattened_id = get_invocation();

    particleBuffer.particles[flattened_id].id.x = flattened_id;
    idBuffer.ids[flattened_id] = flattened_id;
    vec3 gravity = vec3(0, -0.2, 0) * TIME_STEPSIZE2;

    addGravity(vec4(gravity, 0), flattened_id);
    //addWindForce(vec3(0.5, 0, 0.2) * TIME_STEPSIZE2, flattened_id);
    //vertlet(flattened_id);
    //ballCollision(flattened_id);
}