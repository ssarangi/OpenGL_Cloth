#version 430

#define DAMPING 0.01f // how much to damp the cloth simulation each frame
#define TIME_STEPSIZE2 0.5f*0.5f // how large time step each particle takes each frame

layout (local_size_x = 16, local_size_y = 1) in;

uniform vec3 ball_pos;

struct Particle
{
    vec4 position;
    vec4 accumulated_normal;
    vec2 uv;
    vec2 id;
    vec4 old_pos; // the position of the particle in the previous time step, used as part of the verlet numerical integration scheme
    vec4 acceleration; // a vector representing the current acceleration of the particle
    int t1Idx;
    int t2Idx;
    int t3Idx;
    int t4Idx;

    int constraint1;
    int constraint2;
    int constraint3;
    int constraint4;
    int constraint5;
    int constraint6;
    int constraint7;
    int constraint8;
    float constraint1_rest_distance;
    float constraint2_rest_distance;
    float constraint3_rest_distance;
    float constraint4_rest_distance;
    float constraint5_rest_distance;
    float constraint6_rest_distance;
    float constraint7_rest_distance;
    float constraint8_rest_distance;

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
   //uint work_group = gl_WorkGroupID.x * gl_NumWorkGroups.y * gl_NumWorkGroups.z + gl_WorkGroupID.y * gl_NumWorkGroups.z + gl_WorkGroupID.z;
   //return work_group * gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z + gl_LocalInvocationIndex;

   // uint work_group = gl_WorkGroupID.y * gl_NumWorkGroups.x * gl_WorkGroupSize.x + gl_WorkGroupID.x * gl_WorkGroupSize.x + gl_LocalInvocationIndex;
   uint work_group = gl_GlobalInvocationID.y * (gl_NumWorkGroups.x * gl_WorkGroupSize.x) + gl_GlobalInvocationID.x;
   return work_group;
}

void offsetPos(vec4 v, unsigned int particleID)
{
    if (particleBuffer.particles[particleID].movable == 1)
    {
        particleBuffer.particles[particleID].position += v;
    }
}

void vertlet(unsigned int particleID)
{
    if (particleBuffer.particles[particleID].movable == 1)
    {
        vec4 temp = particleBuffer.particles[particleID].position;
        vec4 old_pos = particleBuffer.particles[particleID].old_pos;
        vec4 acceleration = particleBuffer.particles[particleID].acceleration;

        particleBuffer.particles[particleID].position = temp + (temp - old_pos) * (1.0f - DAMPING) + acceleration * TIME_STEPSIZE2;
        particleBuffer.particles[particleID].old_pos = temp;
        particleBuffer.particles[particleID].acceleration = vec4(0, 0, 0, 0);
    }
}

void ball_collision(unsigned int particleID)
{
    float radius = 2.0f;

                //vec4 v = (*particle).getPos() - glm::vec4(center, 0.0);
                //float l = length(v);
                //if (length(v) < radius) // if the particle is inside the ball
                //{
                //    (*particle).offsetPos(normalize(v)*(radius - l)); // project the particle to the surface of the ball
                //}

    vec4 v = particleBuffer.particles[particleID].position - vec4(ball_pos, 0.0);
    float l = length(v);

    if (length(v) < radius)
    {
        offsetPos(normalize(v) * (radius - l), particleID);
    }
}

void main()
{
    uint flattened_id = get_invocation();

    vertlet(flattened_id);
    ball_collision(flattened_id);
}