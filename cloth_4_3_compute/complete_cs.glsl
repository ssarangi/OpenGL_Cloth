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

    uint constraint1;
    uint constraint2;
    uint constraint3;
    uint constraint4;
    uint constraint5;
    uint constraint6;
    uint constraint7;
    uint constraint8;
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

layout (std430, binding = 1) buffer IDBuffer
{
    vec2 ids[];
} idBuffer;

uint get_invocation()
{
   //uint work_group = gl_WorkGroupID.x * gl_NumWorkGroups.y * gl_NumWorkGroups.z + gl_WorkGroupID.y * gl_NumWorkGroups.z + gl_WorkGroupID.z;
   //return work_group * gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z + gl_LocalInvocationIndex;

   // uint work_group = gl_WorkGroupID.y * gl_NumWorkGroups.x * gl_WorkGroupSize.x + gl_WorkGroupID.x * gl_WorkGroupSize.x + gl_LocalInvocationIndex;
   uint work_group = gl_GlobalInvocationID.y * (gl_NumWorkGroups.x * gl_WorkGroupSize.x) + gl_GlobalInvocationID.x;
   return work_group;
}

void addForce(vec4 force_vec, unsigned int particleID)
{
    particleBuffer.particles[particleID].acceleration += force_vec / particleBuffer.particles[particleID].mass;
}

void offsetPos(vec4 v, unsigned int particleID)
{
    if (particleBuffer.particles[particleID].movable == 1)
    {
        particleBuffer.particles[particleID].position += v;
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

// ------------------------------------------------------ Satisfy Constraint -----------------------------------------//
void satisfyConstraintForConnectedConstraint(unsigned int particleID, unsigned int constraintID, float rest_distance)
{
    if (constraintID == -1)
        return;

    vec4 p1Pos = particleBuffer.particles[particleID].position;
    vec4 p2Pos = particleBuffer.particles[constraintID].position;

    vec4 p1_to_p2 = p2Pos - p1Pos;
    float current_distance = length(p1_to_p2);
    vec4 correctionVector = p1_to_p2 * (1 - rest_distance / current_distance);
    vec4 correctionVectorHalf = (correctionVector * 0.5f) * 0.00001;

    barrier();
    // offsetPos(correctionVectorHalf, particleID);
    // offsetPos(correctionVectorHalf, constraintID);
}

void satisfyConstraint(unsigned int particleID)
{
    uint constraint1 = particleBuffer.particles[particleID].constraint1;
    uint constraint2 = particleBuffer.particles[particleID].constraint2;
    uint constraint3 = particleBuffer.particles[particleID].constraint3;
    uint constraint4 = particleBuffer.particles[particleID].constraint4;
    uint constraint5 = particleBuffer.particles[particleID].constraint5;
    uint constraint6 = particleBuffer.particles[particleID].constraint6;
    uint constraint7 = particleBuffer.particles[particleID].constraint7;
    uint constraint8 = particleBuffer.particles[particleID].constraint8;

    float rest_distance_1 = particleBuffer.particles[particleID].constraint1_rest_distance;
    float rest_distance_2 = particleBuffer.particles[particleID].constraint2_rest_distance;
    float rest_distance_3 = particleBuffer.particles[particleID].constraint3_rest_distance;
    float rest_distance_4 = particleBuffer.particles[particleID].constraint4_rest_distance;
    float rest_distance_5 = particleBuffer.particles[particleID].constraint5_rest_distance;
    float rest_distance_6 = particleBuffer.particles[particleID].constraint6_rest_distance;
    float rest_distance_7 = particleBuffer.particles[particleID].constraint7_rest_distance;
    float rest_distance_8 = particleBuffer.particles[particleID].constraint8_rest_distance;

    satisfyConstraintForConnectedConstraint(particleID, constraint1, rest_distance_1);
    satisfyConstraintForConnectedConstraint(particleID, constraint2, rest_distance_2);
    satisfyConstraintForConnectedConstraint(particleID, constraint3, rest_distance_3);
    satisfyConstraintForConnectedConstraint(particleID, constraint4, rest_distance_4);
    satisfyConstraintForConnectedConstraint(particleID, constraint5, rest_distance_5);
    satisfyConstraintForConnectedConstraint(particleID, constraint6, rest_distance_6);
    satisfyConstraintForConnectedConstraint(particleID, constraint7, rest_distance_7);
    satisfyConstraintForConnectedConstraint(particleID, constraint8, rest_distance_8);
}

//--------------------------------------------------------------------------------------------------------------------//

void vertlet(unsigned int particleID)
{
    if (particleBuffer.particles[particleID].movable)
    {
        vec4 temp = particleBuffer.particles[particleID].position;
        vec4 old_pos = particleBuffer.particles[particleID].old_pos;
        vec4 acceleration = particleBuffer.particles[particleID].acceleration;

        particleBuffer.particles[particleID].position = temp + (temp - old_pos) * (1.0f - DAMPING) + acceleration * TIME_STEPSIZE2;
        particleBuffer.particles[particleID].old_pos = temp;
        particleBuffer.particles[particleID].acceleration = vec4(0, 0, 0, 0);
    }
}

void ballCollision(unsigned int particleID)
{

}

layout (local_size_x = 16, local_size_y = 1) in;

void main()
{
    uint flattened_id = get_invocation();

    particleBuffer.particles[flattened_id].id.y = flattened_id;
    idBuffer.ids[flattened_id] = vec2(flattened_id, flattened_id);

    vec4 gravity = vec4(0, -0.2, 0, 0) * TIME_STEPSIZE2;
    addGravity(gravity, flattened_id);
    
    vec4 windForce = vec4(0.5, 0, 0.2, 0.0) * TIME_STEPSIZE2;
    addWindForce(windForce, flattened_id);

    // satisfyConstraint(flattened_id);

    //vertlet(flattened_id);
    //ballCollision(flattened_id);
}