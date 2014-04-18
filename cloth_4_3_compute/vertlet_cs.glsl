#version 430
uniform float roll;

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

struct Constraint
{
    int p1Idx;
    int p2Idx;
    float rest_distance;
};

layout (std430, binding = 0) buffer ParticleBuffer
{
    Particle particles[];
} particleBuffer;

layout (std430, binding = 1) buffer ConstraintBuffer
{
    Constraint constraints[];
} constraintBuffer;

layout (local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

void main()
{
    particleBuffer.particles[gl_GlobalInvocationID.x].position += vec3(0.0f, 0.005f, 0.0f);
}