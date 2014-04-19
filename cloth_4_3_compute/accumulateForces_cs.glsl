#version 430

#define TIME_STEPSIZE2 0.5f*0.5f // how large time step each particle takes each frame

struct Particle
{
    vec3 position;
    vec2 uv;
    vec3 accumulated_normal;

    bool movable; // can the particle move or not ? used to pin parts of the cloth

    float mass; // the mass of the particle (is always 1 in this example)
    vec3 old_pos; // the position of the particle in the previous time step, used as part of the verlet numerical integration scheme
    vec3 acceleration; // a vector representing the current acceleration of the particle
};

layout (std430, binding = 0) buffer ParticleBuffer
{
    Particle particles[];
} particleBuffer;

layout (std430, binding = 1) buffer LockBuffer
{
    uint locks[];
} lockBuffer;


uniform vec3 gravity = vec3(0, -0.2, 0);
uniform vec3 windForce = (vec3(0.5, 0, 0.2);

void addForce(int id)
{
    particleBuffer.particles[id].acceleration += (gravity * TIME_STEPSIZE2) / particleBuffer.particles[id].mass;
}

vec3 calculateTriangleNormal(Particle p1, Particle p2, Particle p3)
{
    vec3 pos1 = p1.position;
    vec3 pos2 = p2.position;
    vec3 pos3 = p3.position;

    vec3 v1 = pos2 - pos1;
    vec3 v2 = pos3 - pos1;

    return cross(v1, v2);
}

vec3 getWindForceForTriangle(int idx1, int idx2, int idx3)
{
    Particle p1 = particleBuffer.particles[idx1];
    Particle p2 = particleBuffer.particles[idx2];
    Particle p3 = particleBuffer.particles[idx3];
    vec3 normal = calculateTriangleNormal(p1, p2, p3);
    vec3 normalizedNormal = normalize(normal);
    vec3 force = normalizedNormal * (dot(normal, windForce));
    return force;
}

void windForce()
{
    int idx1 = gl_GlobalInvocationID.y * gl_WorkGroupSize.x * num_groups_x + gl_GlobalInvocationID.x + 1;
    int idx2 = gl_GlobalInvocationID.y * gl_WorkGroupSize.x * num_groups_x + gl_GlobalInvocationID.x;
    int idx3 = (gl_GlobalInvocationID.y + 1) * gl_WorkGroupSize.x * num_groups_x + gl_GlobalInvocationID.x;

    vec3 force = getWindForceForTriangle(idx1, idx2, idx3);

    while(lockBuffer[idx1] == 1)
    {
        if(lockBuffer[idx1] == 0)
        {
            atomicExchange(lockBuffer[idx1], 1);
            particleBuffer.particles[idx1].acceleration += force / mass;
            atomicExchange(lockBuffer[idx1], 0);
            break;
        }
    }

    atomicAdd(particleBuffer.particles[idx1].acceleration, force / mass);
    atomicAdd(particleBuffer.particles[idx2].acceleration, force / mass);
    atomicAdd(particleBuffer.particles[idx3].acceleration, force / mass);

    barrier();

    idx1 = (gl_GlobalInvocationID.y + 1) * gl_WorkGroupSize.x * num_groups_x + gl_GlobalInvocationID.x + 1;
    idx2 = gl_GlobalInvocationID.y * gl_WorkGroupSize.x * num_groups_x + gl_GlobalInvocationID.x + 1;
    idx3 = (gl_GlobalInvocationID.y + 1) * gl_WorkGroupSize.x * num_groups_x + gl_GlobalInvocationID.x;

    vec3 force = getWindForceForTriangle(idx1, idx2, idx3);

    atomicAdd(particleBuffer.particles[idx1].acceleration, force / mass);
    atomicAdd(particleBuffer.particles[idx2].acceleration, force / mass);
    atomicAdd(particleBuffer.particles[idx3].acceleration, force / mass);
}

layout (local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

void main()
{
    // particleBuffer.particles[gl_GlobalInvocationID.x].position += vec3(0.0f, 0.005f, 0.0f);

    int flattened_id = gl_LocalInvocationID.z * gl_WorkGroupSize.x * num_groups_x * gl_WorkGroupSize.y * num_groups_y +
                        gl_LocalInvocationID.y * gl_WorkGroupSize.x * num_groups_x +
                         gl_LocalInvocationID.x;

    addForce(flattened_id);

    barrier();

    windForce();
}