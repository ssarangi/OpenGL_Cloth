#version 430

#define DAMPING 0.01f // how much to damp the cloth simulation each frame
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

layout (local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

void vertlet(int flattened_id)
{
    if (particleBuffer.particles[flattened_id].movable)
    {
        vec3 temp = particleBuffer.particles[flattened_id].position;
        vec3 old_pos = particleBuffer.particles[flattened_id].old_pos;
        vec3 acceleration = particleBuffer.particles[flattened_id].acceleration;
        particleBuffer.particles[flattened_id].position = temp + (temp - old_pos) * (1.0f - DAMPING) + acceleration * TIME_STEPSIZE2;

        particleBuffer.particles[flattened_id].old_pos = temp;
        particleBuffer.particles[flattened_id].acceleration = vec3(0, 0, 0);
    }
}

void main()
{
    // particleBuffer.particles[gl_GlobalInvocationID.x].position += vec3(0.0f, 0.005f, 0.0f);

    int flattened_id = gl_LocalInvocationID.z * gl_WorkGroupSize.x * gl_WorkGroupSize.y +
                        gl_LocalInvocationID.y * gl_WorkGroupSize.x +
                         gl_LocalInvocationID.x;

    
}