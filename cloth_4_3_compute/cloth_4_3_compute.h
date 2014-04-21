/**
* Based on "Mosegaards Cloth Simulation Coding Tutorial" ( http://cg.alexandra.dk/2009/06/02/mosegaards-cloth-simulation-coding-tutorial/ )
*/

#define GLM_FORCE_RADIANS
#define _USE_MATH_DEFINES
#include "GL/gl3w.h"
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <GL/gl.h>
#include <unordered_map>

#include <cmath>
#include <vector>
#include <iostream>
#include <iterator>
#include "TextResource.h"

namespace cloth_4_3_compute
{
    /* Some physics constants */
#define NUM_PARTICLES_WIDTH 48
#define NUM_PARTICLES_HEIGHT 48

#define DAMPING 0.01f // how much to damp the cloth simulation each frame
#define TIME_STEPSIZE2 0.5f*0.5f // how large time step each particle takes each frame
#define CONSTRAINT_ITERATIONS 15 // how many iterations of constraint satisfaction each frame (more is rigid, less is soft)

    void ExitOnGLError(const char* error_message)
    {
        const GLenum ErrorValue = glGetError();

        if (ErrorValue != GL_NO_ERROR)
        {
            const char* APPEND_DETAIL_STRING = ": %s\n";
            const size_t APPEND_LENGTH = strlen(APPEND_DETAIL_STRING) + 1;
            const size_t message_length = strlen(error_message);
            char* display_message = (char*)malloc(message_length + APPEND_LENGTH);

            memcpy(display_message, error_message, message_length);
            memcpy(&display_message[message_length], APPEND_DETAIL_STRING, APPEND_LENGTH);

            fprintf(stderr, display_message, gluErrorString(ErrorValue));

            free(display_message);
            exit(EXIT_FAILURE);
        }
    }

    using namespace glm;

    class Vertex
    {
    public:
        vec4 position;
        vec4 normal;
        vec2 uv;
    };

    GLuint litShader;
    GLuint unlitShader;

    GLuint computeShader;

    mat4 projection;
    mat4 view;
    vec4 lightPos0; // light position in eye space
    vec4 lightPos1;

    GLuint buildCandyColorTexture(vec4 color1, vec4 color2, int width)
    {
        std::vector<vec4> textureData;
        for (int i = 0; i < width; i++)
        {
            if (i % 2 == 0)
            {
                textureData.push_back(color1);
            }
            else
            {
                textureData.push_back(color2);
            }
        }
        
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        GLuint textureName;
        glGenTextures(1, &textureName);
        glBindTexture(GL_TEXTURE_2D, textureName);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width,
            1, 0, GL_RGBA, GL_FLOAT,
            value_ptr(textureData[0]));
        return textureName;
    }

    /* The particle class represents a particle of mass that can move around in 3D space*/
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

    public:

        Particle(vec4 pos, int x, int y, int width, int height) : old_pos(pos), acceleration(vec4(0, 0, 0, 0)), mass(1), movable(true), accumulated_normal(vec4(0, 0, 0, 0))
        {
            id = vec2(-1, -1);
            position = pos;

            t1Idx = y * width + (x - 1);
            t2Idx = (y - 1) * width + x;
            t3Idx = y * width + (x + 1);
            t4Idx = (y + 1) * width + x;
            constraint1 = -1;
            constraint2 = -1;
            constraint3 = -1;
            constraint4 = -1;
            constraint5 = -1;
            constraint6 = -1;
            constraint7 = -1;
            constraint8 = -1;
            constraint1_rest_distance = 0.0f;
            constraint2_rest_distance = 0.0f;
            constraint3_rest_distance = 0.0f;
            constraint4_rest_distance = 0.0f;
            constraint5_rest_distance = 0.0f;
            constraint6_rest_distance = 0.0f;
            constraint7_rest_distance = 0.0f;
            constraint8_rest_distance = 0.0f;
        }

        Particle(){}

        void addForce(vec3 f)
        {
            acceleration += glm::vec4(f / mass, 0);
        }

        /* This is one of the important methods, where the time is progressed a single step size (TIME_STEPSIZE)
        The method is called by Cloth.time_step()
        Given the equation "force = mass * acceleration" the next position is found through verlet integration*/
        void timeStep()
        {
            if (movable)
            {
                vec4 temp = position;
                position = position + (position - old_pos)*(1.0f - DAMPING) + acceleration*TIME_STEPSIZE2;
                old_pos = temp;
                acceleration = vec4(0, 0, 0, 0); // acceleration is reset since it HAS been translated into a change in position (and implicitely into velocity)	
            }
        }

        vec4& getPos() { return position; }

        void resetAcceleration() { acceleration = vec4(0, 0, 0, 0); }

        void offsetPos(const vec4 v) { if (movable) position += v; }

        void makeUnmovable() { movable = false; }

        void addToNormal(vec4 normal)
        {
            accumulated_normal += normalize(normal);
        }

        void setUV(vec2 uv) { this->uv = uv; }

        vec4& getNormal() { return accumulated_normal; } // notice, the normal is not unit length

        void resetNormal() { accumulated_normal = vec4(0, 0, 0, 0); }

    };

    class Constraint
    {
    public:
        int p1Idx;
        int p2Idx;
        float rest_distance; // the length between particle p1 and p2 in rest configuration

        Constraint(Particle *p1, Particle *p2, uint p1idx, uint p2idx) : p1Idx(p1idx), p2Idx(p2idx)
        {
            vec4 vec = p1->getPos() - p2->getPos();
            rest_distance = length(vec);
        }
    };

    struct Cloth
    {
    public:
        GLuint vertexArrayObject = 0;
        GLuint vertexBuffer = 0;
        GLuint texture;
        int elementSize;

        int num_particles_width; // number of particles in "width" direction
        int num_particles_height; // number of particles in "height" direction
        // total number of particles is num_particles_width*num_particles_height

        std::vector<Particle> particles; // all particles that are part of this cloth
        std::vector<Constraint> constraints; // all constraints between particles as part of this cloth
        std::unordered_map<uint, Particle*> m_index2particle;
        std::unordered_map<Particle*, uint> m_particle2index;

        int getParticleIndex(int x, int y) { return y*num_particles_width + x; }

        Particle* getParticle(int x, int y) { return &particles[getParticleIndex(x, y)]; }

        void addConstraintToParticle(Particle* p, int index, float rest_distance)
        {
            if (p->constraint1 == -1)
            {
                p->constraint1 = index;
                p->constraint1_rest_distance = rest_distance;
            }
            else if (p->constraint2 == -1)
            {
                p->constraint2 = index;
                p->constraint2_rest_distance = rest_distance;
            }
            else if (p->constraint3 == -1)
            {
                p->constraint3 = index;
                p->constraint3_rest_distance = rest_distance;
            }
            else if (p->constraint4 == -1)
            {
                p->constraint4 = index;
                p->constraint4_rest_distance = rest_distance;
            }
            else if (p->constraint5 == -1)
            {
                p->constraint5 = index;
                p->constraint5_rest_distance = rest_distance;
            }
            else if (p->constraint6 == -1)
            {
                p->constraint6 = index;
                p->constraint6_rest_distance = rest_distance;
            }
            else if (p->constraint7 == -1)
            {
                p->constraint7 = index;
                p->constraint7_rest_distance = rest_distance;
            }
            else if (p->constraint8 == -1)
            {
                p->constraint8 = index;
                p->constraint8_rest_distance = rest_distance;
            }
        }

        void makeConstraint(Particle *p1, Particle *p2)
        {
            uint p1Idx = m_particle2index[p1];
            uint p2Idx = m_particle2index[p2];

            Constraint c = Constraint(p1, p2, p1Idx, p2Idx);

            addConstraintToParticle(p1, p2Idx, c.rest_distance);
            // addConstraintToParticle(p2, p1Idx, c.rest_distance);
            constraints.push_back(c);
        }

        /* This is one of the important methods, where a single constraint between two particles p1 and p2 is solved
        the method is called by Cloth.time_step() many times per frame*/
        void satisfyConstraint(Constraint* pConstraint)
        {
            Particle* p1 = m_index2particle[pConstraint->p1Idx];
            Particle* p2 = m_index2particle[pConstraint->p2Idx];
            vec4 p1_to_p2 = p2->getPos() - p1->getPos(); // vector from p1 to p2
            float current_distance = length(p1_to_p2); // current distance between p1 and p2
            vec4 correctionVector = p1_to_p2 * (1 - pConstraint->rest_distance / current_distance); // The offset vector that could moves p1 into a distance of rest_distance to p2
            vec4 correctionVectorHalf = correctionVector*0.5f; // Lets make it half that length, so that we can move BOTH p1 and p2.
            p1->offsetPos(correctionVectorHalf); // correctionVectorHalf is pointing from p1 to p2, so the length should move p1 half the length needed to satisfy the constraint.
            p2->offsetPos(-correctionVectorHalf); // we must move p2 the negative direction of correctionVectorHalf since it points from p2 to p1, and not p1 to p2.	
        }


        /* A private method used by drawShaded() and addWindForcesForTriangle() to retrieve the
        normal vector of the triangle defined by the position of the particles p1, p2, and p3.
        The magnitude of the normal vector is equal to the area of the parallelogram defined by p1, p2 and p3
        */
        vec3 calcTriangleNormal(Particle *p1, Particle *p2, Particle *p3)
        {
            vec4 pos1 = p1->getPos();
            vec4 pos2 = p2->getPos();
            vec4 pos3 = p3->getPos();

            vec4 v1 = pos2 - pos1;
            vec4 v2 = pos3 - pos1;

            return cross(glm::vec3(v1), glm::vec3(v2));
        }

        /* A private method used by windForce() to calcualte the wind force for a single triangle
        defined by p1,p2,p3*/
        void addWindForcesForTriangle(Particle *p1, Particle *p2, Particle *p3, const vec3 direction)
        {
            vec3 normal = calcTriangleNormal(p1, p2, p3);
            vec3 d = normalize(normal);
            vec3 force = normal*(dot(d, direction));
            p1->addForce(force);
            p2->addForce(force);
            p3->addForce(force);
        }

        /* A private method used by drawShaded(), that draws a single triangle p1,p2,p3 with a color*/
        void insertTriangle(Particle *p1, const vec2 uv, std::vector<Vertex> &vertexData)
        {
            Vertex v1 = { p1->getPos(), p1->getNormal(), uv };
            vertexData.push_back(v1);
        }

    public:
        GLuint vertex_vbo_storage;
        GLuint id_storage;

        /* This is a important constructor for the entire system of particles and constraints*/
        Cloth(float width, float height, int num_particles_width, int num_particles_height) : num_particles_width(num_particles_width), num_particles_height(num_particles_height)
        {
            particles.resize(num_particles_width*num_particles_height); //I am essentially using this vector as an array with room for num_particles_width*num_particles_height particles

            // creating particles in a grid of particles from (0,0,0) to (width,-height,0)
            for (int x = 0; x < num_particles_width; x++)
            {
                for (int y = 0; y < num_particles_height; y++)
                {
                    vec3 pos = vec3(width * (x / (float)num_particles_width),
                        -height * (y / (float)num_particles_height),
                        0);
                    
                    particles[y*num_particles_width + x] = Particle(glm::vec4(pos.x, pos.y, pos.z, 0.0), x, y, num_particles_width, num_particles_height); // insert particle in column x at y'th row
                    particles[y*num_particles_width + x].id.x = y * num_particles_width + x;
                    m_index2particle[y*num_particles_width + x] = &particles[y*num_particles_width + x];
                    m_particle2index[&particles[y*num_particles_width + x]] = y*num_particles_width + x;
                }
            }

            for (int y = 0; y < num_particles_height; y++)
            {
                for (int x = 0; x < num_particles_width; x++)
                {
                    vec2 uv(x / (num_particles_width - 1.0f), y / (num_particles_height - 1.0f));

                    Particle* p = getParticle(x, y);
                    p->setUV(uv);
                }
            }

            // Connecting immediate neighbor particles with constraints (distance 1 and sqrt(2) in the grid)
            for (int x = 0; x < num_particles_width; x++)
            {
                for (int y = 0; y < num_particles_height; y++)
                {
                    if (x < num_particles_width - 1) makeConstraint(getParticle(x, y), getParticle(x + 1, y));
                    if (y < num_particles_height - 1) makeConstraint(getParticle(x, y), getParticle(x, y + 1));
                    if (x < num_particles_width - 1 && y < num_particles_height - 1) makeConstraint(getParticle(x, y), getParticle(x + 1, y + 1));
                    if (x < num_particles_width - 1 && y < num_particles_height - 1) makeConstraint(getParticle(x + 1, y), getParticle(x, y + 1));
                }
            }


            // Connecting secondary neighbors with constraints (distance 2 and sqrt(4) in the grid)
            for (int x = 0; x < num_particles_width; x++)
            {
                for (int y = 0; y < num_particles_height; y++)
                {
                    if (x < num_particles_width - 2) makeConstraint(getParticle(x, y), getParticle(x + 2, y));
                    if (y < num_particles_height - 2) makeConstraint(getParticle(x, y), getParticle(x, y + 2));
                    if (x < num_particles_width - 2 && y < num_particles_height - 2) makeConstraint(getParticle(x, y), getParticle(x + 2, y + 2));
                    if (x < num_particles_width - 2 && y < num_particles_height - 2) makeConstraint(getParticle(x + 2, y), getParticle(x, y + 2));
                }
            }


            // making the upper left most three and right most three particles unmovable
            for (int i = 0; i < 3; i++)
            {
                getParticle(0 + i, 0)->offsetPos(vec4(0.5, 0.0, 0.0, 0.0)); // moving the particle a bit towards the center, to make it hang more natural - because I like it ;)
                getParticle(0 + i, 0)->makeUnmovable();

                getParticle(0 + i, 0)->offsetPos(vec4(-0.5, 0.0, 0.0, 0.0)); // moving the particle a bit towards the center, to make it hang more natural - because I like it ;)
                getParticle(num_particles_width - 1 - i, 0)->makeUnmovable();
            }
        }

        ~Cloth()
        {
            //for (int i = 0; i < particles.size(); ++i)
            //{
            //    delete particles[i];
            //}

            //particles.clear();
        }

        /* drawing the cloth as a smooth shaded (and colored according to column) OpenGL triangular mesh
        Called from the display() method
        The cloth is seen as consisting of triangles for four particles in the grid as follows:

        (x,y)   *--* (x+1,y)
        | /|
        |/ |
        (x,y+1) *--* (x+1,y+1)

        */
        void drawShaded()
        {
            // reset normals (which where written to last frame)
            std::vector<Particle>::iterator particle;
            for (particle = particles.begin(); particle != particles.end(); particle++)
            {
                (*particle).resetNormal();
            }

            //create smooth per particle normals by adding up all the (hard) triangle normals that each particle is part of
            for (int x = 0; x < num_particles_width - 1; x++)
            {
                for (int y = 0; y < num_particles_height - 1; y++)
                {
                    vec3 normal = calcTriangleNormal(getParticle(x + 1, y), getParticle(x, y), getParticle(x, y + 1));
                    getParticle(x + 1, y)->addToNormal(glm::vec4(normal, 0.0));
                    getParticle(x, y)->addToNormal(glm::vec4(normal, 0.0));
                    getParticle(x, y + 1)->addToNormal(glm::vec4(normal, 0.0));

                    normal = calcTriangleNormal(getParticle(x + 1, y + 1), getParticle(x + 1, y), getParticle(x, y + 1));
                    getParticle(x + 1, y + 1)->addToNormal(glm::vec4(normal, 0.0));
                    getParticle(x + 1, y)->addToNormal(glm::vec4(normal, 0.0));
                    getParticle(x, y + 1)->addToNormal(glm::vec4(normal, 0.0));
                }
            }

            if (vertexArrayObject == 0)
            {
                glGenVertexArrays(1, &vertexArrayObject);
                glBindVertexArray(vertexArrayObject);

                glGenBuffers(1, &vertexBuffer);
                glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);

                GLuint positionAttributeLocation = glGetAttribLocation(litShader, "position");
                GLuint uvAttributeLocation = glGetAttribLocation(litShader, "uv");
                GLuint normalAttributeLocation = glGetAttribLocation(litShader, "normal");
                glEnableVertexAttribArray(positionAttributeLocation);
                glEnableVertexAttribArray(uvAttributeLocation);
                glEnableVertexAttribArray(normalAttributeLocation);
                glVertexAttribPointer(positionAttributeLocation, 4, GL_FLOAT, GL_FALSE, sizeof(Particle), (const GLvoid *)0);
                glVertexAttribPointer(normalAttributeLocation, 4, GL_FLOAT, GL_FALSE, sizeof(Particle), (const GLvoid *)sizeof(vec4));
                glVertexAttribPointer(uvAttributeLocation, 2, GL_FLOAT, GL_FALSE, sizeof(Particle), (const GLvoid *)(sizeof(vec4)+sizeof(vec4)));

                std::vector<int> indices;


                for (int j = 0; j < num_particles_height - 1; j++)
                {
                    int index;
                    if (j > 0)
                    {
                        indices.push_back(j * num_particles_width); // make degenerate
                    }
                    
                    for (int i = 0; i <= num_particles_width - 1; i++)
                    {
                        index = j * num_particles_width + i;
                        indices.push_back(index);
                        indices.push_back(index + num_particles_width);
                    }
                    
                    if (j + 1 < num_particles_height - 1)
                    {
                        indices.push_back(index + num_particles_width); // make degenerate
                    }
                }

                elementSize = indices.size();

                GLuint elementArrayBuffer;
                glGenBuffers(1, &elementArrayBuffer);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementArrayBuffer);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, elementSize * sizeof(int), &(indices[0]), GL_STATIC_DRAW);

                vec4 color1 = vec4(1.0f, 1.0f, 1.0f, 1.0f);
                vec4 color2 = vec4(0.6f, 0.2f, 0.2f, 1.0f);
                texture = buildCandyColorTexture(color1, color2, num_particles_width - 1);
            }

            glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
            glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(Particle), &particles[0], GL_DYNAMIC_COPY);

            mat4 modelView = view;
            mat4 mvp = projection * modelView;
            glUniformMatrix4fv(glGetUniformLocation(litShader, "mvp"), 1, false, value_ptr(mvp));
            mat3 normalMatrix = inverse(transpose(mat3(modelView)));
            glUniformMatrix3fv(glGetUniformLocation(litShader, "normalMatrix"), 1, false, value_ptr(normalMatrix));

            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, texture);
            glUniform1i(glGetUniformLocation(litShader, "mainTexture"), 0);

            glBindVertexArray(vertexArrayObject);
            glDrawElements(GL_TRIANGLE_STRIP, elementSize, GL_UNSIGNED_INT, 0);
        }

        /* this is an important methods where the time is progressed one time step for the entire cloth.
        This includes calling satisfyConstraint() for every constraint, and calling timeStep() for all particles
        */
        void timeStep()
        {
            std::vector<Constraint>::iterator constraint;
            for (int i = 0; i < CONSTRAINT_ITERATIONS; i++) // iterate over all constraints several times
            {
                for (constraint = constraints.begin(); constraint != constraints.end(); constraint++)
                {
                    satisfyConstraint(&(*constraint)); // satisfy constraint.
                }
            }

            std::vector<Particle>::iterator particle;
            for (particle = particles.begin(); particle != particles.end(); particle++)
            {
                (*particle).timeStep(); // calculate the position of each particle at the next time step.
            }
        }

        /* used to add gravity (or any other arbitrary vector) to all particles*/
        void addForce(const vec3 direction)
        {
            std::vector<Particle>::iterator particle;
            for (particle = particles.begin(); particle != particles.end(); particle++)
            {
                (*particle).addForce(direction); // add the forces to each particle
            }

        }

        /* used to add wind forces to all particles, is added for each triangle since the final force is proportional to the triangle area as seen from the wind direction*/
        void windForce(const vec3 direction)
        {
            for (int x = 0; x < num_particles_width - 1; x++)
            {
                for (int y = 0; y < num_particles_height - 1; y++)
                {
                    addWindForcesForTriangle(getParticle(x + 1, y), getParticle(x, y), getParticle(x, y + 1), direction);
                    addWindForcesForTriangle(getParticle(x + 1, y + 1), getParticle(x + 1, y), getParticle(x, y + 1), direction);
                }
            }
        }

        /* used to detect and resolve the collision of the cloth with the ball.
        This is based on a very simples scheme where the position of each particle is simply compared to the sphere and corrected.
        This also means that the sphere can "slip through" if the ball is small enough compared to the distance in the grid bewteen particles
        */
        void ballCollision(const vec3 center, const float radius)
        {
            std::vector<Particle>::iterator particle;
            for (particle = particles.begin(); particle != particles.end(); particle++)
            {
                vec4 v = (*particle).getPos() - glm::vec4(center, 0.0);
                float l = length(v);
                if (length(v) < radius) // if the particle is inside the ball
                {
                    (*particle).offsetPos(normalize(v)*(radius - l)); // project the particle to the surface of the ball
                }
            }
        }

        void doFrame()
        {

        }
    };

    /***** Above are definition of classes; vec3, Particle, Constraint, and Cloth *****/




    // Just below are three global variables holding the actual animated stuff; Cloth and Ball 
    Cloth cloth1(14, 10, NUM_PARTICLES_WIDTH, NUM_PARTICLES_HEIGHT); // one Cloth object of the Cloth class
    vec3 ball_pos(7, -5, 0); // the center of our one ball
    float ball_radius = 2; // the radius of our one ball



    /***** Below are functions Init(), display(), reshape(), keyboard(), arrow_keys(), main() *****/

    /* This is where all the standard Glut/OpenGL stuff is, and where the methods of Cloth are called;
    addForce(), windForce(), timeStep(), ballCollision(), and drawShaded()*/

    void init(GLvoid)
    {
        // glShadeModel(GL_SMOOTH);
        glClearColor(0.2f, 0.2f, 0.4f, 0.5f);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        lightPos0 = vec4(-1.0, 1.0, 0.5, 0.0);
        vec4 lightAmbient0 = vec4(0.2, 0.2, 0.2, 1.0);
        vec4 lightDiffuse0 = vec4(0.8, 0.8, 0.8, 1.0);

        lightPos1 = vec4(1.0, 0.0, -0.2, 0.0);
        vec4 lightAmbient1 = vec4(0.0, 0.0, 0.0, 0.0);
        vec4 lightDiffuse1 = vec4(0.5, 0.5, 0.3, 0.0);

        vec4 ambient[2] = { lightAmbient0, lightAmbient1 };
        vec4 diffuse[2] = { lightDiffuse0, lightDiffuse1 };
        glUseProgram(litShader);
        glUniform4fv(glGetUniformLocation(litShader, "lightAmbient"), 2, value_ptr(ambient[0]));
        glUniform4fv(glGetUniformLocation(litShader, "lightDiffuse"), 2, value_ptr(diffuse[0]));

        vec4 lightModelAmbient = vec4(0.2, 0.2, 0.2, 1.0);
        glUniform4fv(glGetUniformLocation(litShader, "lightModelAmbient"), 1, value_ptr(lightModelAmbient));

        glGenBuffers(1, &cloth1.vertex_vbo_storage);
        glGenBuffers(1, &cloth1.id_storage);
    }


    float ball_time = 0; // counter for used to calculate the z position of the ball below

    void drawSolidSphere(vec3& position)
    {
        glUseProgram(litShader);
        static GLuint vertexArrayObject = 0;
        static int elementCount;
        static GLuint sphereTex;
        if (vertexArrayObject == 0)
        {
            std::vector<Vertex> vertexData;
            int slices = 64;
            int stacks = 32;
            float radius = 1.9;
            int vertexCount = (stacks + 1) * (slices + 1);
            float piDivStacks = M_PI / stacks;
            float PIDiv2 = M_PI / 2;
            float PI2 = M_PI * 2;

            for (int j = 0; j <= stacks; j++)
            {
                float latitude1 = piDivStacks * j - PIDiv2;
                float sinLat1 = sin(latitude1);
                float cosLat1 = cos(latitude1);
                for (int i = 0; i <= slices; i++)
                {
                    float longitude = (PI2 / slices) * i;
                    float sinLong = sin(longitude);
                    float cosLong = cos(longitude);
                    vec3 normal = vec3(cosLong * cosLat1, sinLat1, sinLong * cosLat1);
                    vec3 position = normal * radius;
                    Vertex v = 
                    { 
                        glm::vec4(position, 0.0), 
                        glm::vec4(normal, 0.0),
                        vec2(j / (float)stacks, i / (float)slices)
                    };

                    vertexData.push_back(v);
                }
            }
            std::vector<GLuint> indices;
            // create indices
            for (int j = 0; j < stacks; j++)
            {
                int index;
                if (j > 0)
                {
                    indices.push_back(j * (slices + 1)); // make degenerate
                }
                
                for (int i = 0; i <= slices; i++)
                {
                    index = j * (slices + 1) + i;
                    indices.push_back(index);
                    indices.push_back(index + slices + 1);
                }
                
                if (j + 1 < stacks)
                {
                    indices.push_back(index + slices + 1); // make degenerate
                }
            }

            glGenVertexArrays(1, &vertexArrayObject);
            glBindVertexArray(vertexArrayObject);

            GLuint vertexBuffer;
            glGenBuffers(1, &vertexBuffer);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
            glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(Vertex), value_ptr(vertexData[0].position), GL_STATIC_DRAW);

            GLuint positionAttributeLocation = glGetAttribLocation(litShader, "position");
            GLuint uvAttributeLocation = glGetAttribLocation(litShader, "uv");
            GLuint normalAttributeLocation = glGetAttribLocation(litShader, "normal");
            glEnableVertexAttribArray(positionAttributeLocation);
            glEnableVertexAttribArray(uvAttributeLocation);
            glEnableVertexAttribArray(normalAttributeLocation);
            glVertexAttribPointer(positionAttributeLocation, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid *)0);
            glVertexAttribPointer(normalAttributeLocation, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid *)sizeof(vec4));
            glVertexAttribPointer(uvAttributeLocation, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid *)(sizeof(vec4)+sizeof(vec4)));

            GLuint elementArrayBuffer;
            glGenBuffers(1, &elementArrayBuffer);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementArrayBuffer);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), &(indices[0]), GL_STATIC_DRAW);
            elementCount = indices.size();
            vec4 color = vec4(0.4, 0.8, 0.5, 1.0);
            sphereTex = buildCandyColorTexture(color, color, 1);
        }

        mat4 modelView = view;
        modelView = translate(modelView, position);
        mat4 mvp = projection * modelView;
        glUniformMatrix4fv(glGetUniformLocation(litShader, "mvp"), 1, false, value_ptr(mvp));
        mat3 normalMatrix = inverse(transpose(mat3(modelView)));
        glUniformMatrix3fv(glGetUniformLocation(litShader, "normalMatrix"), 1, false, value_ptr(normalMatrix));

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, sphereTex);
        glUniform1i(glGetUniformLocation(litShader, "mainTexture"), 0);

        glBindVertexArray(vertexArrayObject);
        glDrawElements(GL_TRIANGLE_STRIP, elementCount, GL_UNSIGNED_INT, 0);
    }

    /* display method called each frame*/
    void display()
    {
        bool verify = true;

        // calculating positions
        ball_time++;
        ball_pos[2] = (float)cos(ball_time / 50.0f) * 7;

        vec2 gg[2304];
        memset(gg, 0, 2304 * sizeof(vec2));

        std::vector<Particle> particle_copy;
        particle_copy.resize(cloth1.particles.size());
        std::copy(cloth1.particles.begin(), cloth1.particles.end(), particle_copy.begin());
        
        if (verify)
            cloth1.addForce(vec3(0, -0.2, 0) * TIME_STEPSIZE2); // add gravity each frame, pointing down

        glUseProgram(computeShader);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, cloth1.vertex_vbo_storage);
        
        if (verify)
            glBufferData(GL_SHADER_STORAGE_BUFFER, cloth1.particles.size() * sizeof(Particle), &(particle_copy[0]), GL_DYNAMIC_COPY);
        else
            glBufferData(GL_SHADER_STORAGE_BUFFER, cloth1.particles.size() * sizeof(Particle), &(cloth1.particles[0]), GL_DYNAMIC_COPY);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, cloth1.id_storage);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 2304 * sizeof(vec2), &(gg[0]), GL_DYNAMIC_COPY);

        glDispatchCompute(144, 1, 1);

        glMemoryBarrier(GL_ALL_BARRIER_BITS);

        {
            GLenum err = gl3wGetError();
        }

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, 0);

        glBindBuffer(GL_ARRAY_BUFFER, cloth1.vertex_vbo_storage);
        Particle * ptr = reinterpret_cast<Particle *>(glMapBufferRange(GL_ARRAY_BUFFER, 0, cloth1.particles.size() * sizeof(Particle), GL_MAP_READ_BIT));

        if (verify)
            memcpy(&particle_copy[0], ptr, particle_copy.size()*sizeof(Particle));
        else
            memcpy(&cloth1.particles[0], ptr, cloth1.particles.size()*sizeof(Particle));

        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glBindBuffer(GL_ARRAY_BUFFER, cloth1.id_storage);
        int * i = reinterpret_cast<int *>(glMapBufferRange(GL_ARRAY_BUFFER, 0, 2304 * sizeof(vec2), GL_MAP_READ_BIT));

        memcpy(&gg, i, 2304*sizeof(vec2));
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);


        if (verify)
        {
            for (int i = 0; i < particle_copy.size(); ++i)
            {
                //std::cout << particle_copy[i].acceleration.x << " <---> " << cloth1.particles[i].acceleration.x << std::endl;
                //std::cout << particle_copy[i].acceleration.y << " <---> " << cloth1.particles[i].acceleration.y << std::endl;
                //std::cout << particle_copy[i].acceleration.z << " <---> " << cloth1.particles[i].acceleration.z << std::endl;
                //std::cout << particle_copy[i].acceleration.w << " <---> " << cloth1.particles[i].acceleration.w << std::endl;
                //assert(particle_copy[i].position == cloth1.particles[i].position);
                assert(particle_copy[i].acceleration == cloth1.particles[i].acceleration);
            }
        }

        // cloth1.windForce(vec3(0.5, 0, 0.2) * TIME_STEPSIZE2); // generate some wind each frame
        cloth1.timeStep(); // calculate the particle positions of the next frame
        cloth1.ballCollision(ball_pos, ball_radius); // resolve collision with the ball

        // drawing
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        view = mat4(1.0f);
        view = translate(view, vec3(-6.5, 6, -15.0f));
        view = rotate(view, 25.0f, vec3(0, 1, 0));
        
        // setup light 
        glUseProgram(litShader);

        vec4 eyeSpaceLight[2] = { lightPos0, lightPos1 };
        glUniform4fv(glGetUniformLocation(litShader, "lightPosition"), 2, value_ptr(eyeSpaceLight[0]));
        cloth1.drawShaded();

        drawSolidSphere(ball_pos);

        glutSwapBuffers();
        glutPostRedisplay();
    }

    void reshape(int w, int h)
    {
        glViewport(0, 0, w, h);
        if (h == 0)
            projection = perspective(80.0f, (float)w, 1.0f, 5000.0f);
        else
            projection = perspective(80.0f, (float)w / (float)h, 1.0f, 5000.0f);
    }

    void keyboard(unsigned char key, int x, int y)
    {
        switch (key)
        {
        case 27:
            exit(0);
            break;
        case 'w':
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            break;
        case 'W':
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            break;
        default:
            break;
        }
    }

    void arrow_keys(int a_keys, int x, int y)
    {
        switch (a_keys)
        {
        case GLUT_KEY_UP:
            glutFullScreen();
            break;
        case GLUT_KEY_DOWN:
            glutReshapeWindow(1280, 720);
            break;
        default:
            break;
        }
    }

    void checkCompileStatus(GLuint shader, const char *shadername)
    {
        GLint  compiled;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
        if (!compiled)
        {
            std::cerr << shadername << " failed to compile:" << std::endl;
            GLint  logSize;
            glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logSize);
            char* logMsg = new char[logSize];
            glGetShaderInfoLog(shader, logSize, NULL, logMsg);
            std::cerr << logMsg << std::endl;
            delete[] logMsg;
        }
    }

    void checkLinkStatus(GLuint program, const char * programName)
    {
        GLint  linked;
        glGetProgramiv(program, GL_LINK_STATUS, &linked);
        if (!linked)
        {
            std::cerr << "Shader program " << programName << " failed to link" << std::endl;
            GLint  logSize;
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logSize);
            char* logMsg = new char[logSize];
            glGetProgramInfoLog(program, logSize, NULL, logMsg);
            std::cerr << logMsg << std::endl;
            delete[] logMsg;
            system("Pause");
            exit(0);
        }
    }

    GLuint loadShader(const char* vertexShaderName, const char* fragmentShaderName)
    {
        const char *vertexShaderSource = TextResource::load(vertexShaderName);
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);
        checkCompileStatus(vertexShader, vertexShaderName);

        const char *fragmentShaderSource = TextResource::load(fragmentShaderName);
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);
        checkCompileStatus(fragmentShader, fragmentShaderName);

        GLuint program = glCreateProgram();
        glAttachShader(program, vertexShader);
        glAttachShader(program, fragmentShader);
        glLinkProgram(program);
        checkLinkStatus(program, vertexShaderName);
        return program;
    }

    GLuint loadComputeShader(const char* computeShaderName)
    {
        const char *computeShaderSource = TextResource::load(computeShaderName);
        GLuint computeShader = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(computeShader, 1, &computeShaderSource, NULL);
        glCompileShader(computeShader);
        checkCompileStatus(computeShader, computeShaderName);

        GLuint program = glCreateProgram();
        glAttachShader(program, computeShader);
        glLinkProgram(program);
        checkLinkStatus(program, computeShaderName);
        return program;
    }

    int main(int &argc, char** argv)
    {
        glutInit(&argc, argv);
        glutInitContextVersion(4, 2);
        glutInitContextProfile(GLUT_CORE_PROFILE);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
        glutInitWindowSize(1280, 720);

        glutCreateWindow("Cloth Tutorial Refactoring OpenGL 4.3 Compute");

        if (gl3wInit())
        {
            fprintf(stderr, "failed to initialize OpenGL\n");
            return -1;
        }
        
        if (!gl3wIsSupported(4, 2)) {
            fprintf(stderr, "OpenGL 4.2 not supported\n");
            return -1;
        }

        litShader = loadShader("../cloth_4_3_compute/lambert.vert", "../cloth_4_3_compute/lambert.frag");
        unlitShader = loadShader("../cloth_4_3_compute/unlit.vert", "../cloth_4_3_compute/unlit.frag");
        computeShader = loadComputeShader("../cloth_4_3_compute/complete_cs.glsl");
        init();

        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        glutKeyboardFunc(keyboard);
        glutSpecialFunc(arrow_keys);

        glutMainLoop();

        return 0;
    }

}
