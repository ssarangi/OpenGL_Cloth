#version 430

layout (local_size_x = 16, local_size_y = 1) in;

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

struct Constraint
{
    int p1Idx;
    int p2Idx;
    float rest_distance;
};

layout (std430, binding = 1) buffer ConstraintBuffer
{
    Constraint constraints[];
} constraintBuffer;

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

void satisfyConstraintForConnectedConstraint(unsigned int particleID, unsigned int constraintID, float rest_distance)
{
    if (constraintID == -1)
        return;

    vec4 p1Pos = particleBuffer.particles[particleID].position;
    vec4 p2Pos = particleBuffer.particles[constraintID].position;

    vec4 p1_to_p2 = p2Pos - p1Pos;
    float current_distance = length(p1_to_p2);
    vec4 correctionVector = p1_to_p2 * (1 - rest_distance / current_distance);
    vec4 correctionVectorHalf = correctionVector * 0.5f;

    barrier();
    offsetPos(correctionVectorHalf, particleID);
    offsetPos(-correctionVectorHalf, constraintID);
    memoryBarrier();
    barrier();
}

void satisfyConstraint(unsigned int constraintID)
{
    Constraint c = constraintBuffer.constraints[constraintID];
    satisfyConstraintForConnectedConstraint(c.p1Idx, c.p2Idx, c.rest_distance);
}

void main()
{
    uint flattened_id = get_invocation();
    satisfyConstraint(flattened_id);
}
