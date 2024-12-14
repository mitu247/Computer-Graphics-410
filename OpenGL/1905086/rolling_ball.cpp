#include <bits/stdc++.h>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <GL/glut.h>
#endif

#define DEG_INC 10 // perclick degree increase of arrow
#define DEG_ROT 15 // perclick degree increase of sphere (rotation of the ball)
int angleOfArrow = 0;

struct TSF
{
    GLdouble t, max_t;
    TSF()
    {
        t = 1.0;
        max_t = 1.0;
    }
};
struct PT
{
    double x, y, z;

    PT() {}
    PT(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    struct PT operator=(const struct PT &p)
    {
        this->x = p.x;
        this->y = p.y;
        this->z = p.z;

        return *this;
    }

    // vector + vector
    struct PT operator+(const struct PT &p)
    {
        double newx = this->x + p.x;
        double newy = this->y + p.y;
        double newz = this->z + p.z;

        return PT(newx, newy, newz);
    }

    // vector * number
    struct PT operator*(const double num)
    {
        double newx = this->x * num;
        double newy = this->y * num;
        double newz = this->z * num;

        return PT(newx, newy, newz);
    }

    // vector / number
    struct PT operator/(const double num)
    {
        double newx = this->x / num;
        double newy = this->y / num;
        double newz = this->z / num;

        return PT(newx, newy, newz);
    }

    // vector x vector
    struct PT operator*(const struct PT &p)
    {
        double newx = (this->y * p.z) - (p.y * this->z);
        double newy = (p.x * this->z) - (this->x * p.z);
        double newz = (this->x * p.y) - (p.x * this->y);

        return PT(newx, newy, newz);
    }

    double dot_product(struct PT vect)
    {
        return (this->x * vect.x) + (this->y * vect.y) + (this->z * vect.z);
    }
};

struct Sphere
{
    double x, y, z;
    double radius;
    bool goForward;
    bool ballInit, simulation, xcollision, ycollision;
    double wallheight, walllength;
    struct PT direction_vector;
    Sphere()
    {
        x = y = z = 0;
        radius = 1.0; // default radius
        goForward = ballInit = true;
        xcollision = ycollision = false;
        wallheight = 1.7;
        walllength = 8.0;
        simulation = false;
        // wrt x-axis
        struct PT t(cos(0), sin(0), 0);
        direction_vector = t;
    }
    Sphere(double radius)
    {
        x = y = z = 0;
        this->radius;
        goForward = ballInit = true;
        xcollision = ycollision = false;
        wallheight = 1.7;
        walllength = 8.0;
        simulation = false;
        // wrt x-axis
        struct PT t(cos(0), sin(0), 0);
        direction_vector = t;
    }
};
int callCount = 0;

struct Camera
{
    struct PT pos, right, up, look;

    Camera() {}

    Camera(struct PT pos, struct PT right, struct PT up, struct PT look)
    {
        this->pos = pos;
        this->right = right;
        this->up = up;
        this->look = look;
    }
};

struct Camera cam;
GLdouble sectorCount = 50;
GLdouble stackCount = 50;
struct PT points[500][500];

struct Sphere sphere;

double deg2rad(double deg)
{
    return (deg * M_PI) / 180;
}

void normalize_vector(struct PT &vect)
{
    vect = vect / sqrt((vect.x * vect.x) + (vect.y * vect.y) + (vect.z * vect.z));
}
void rotation_orth(struct PT &v1, struct PT &v2, struct PT a, double an)
{
    v1 = v1 * cos(an) + (a * v1) * sin(an);
    v2 = v2 * cos(an) + (a * v2) * sin(an);
}
void rotate_general(struct PT &v, PT a, double angle)
{
    // Rodrigues Formulae
    angle = deg2rad(angle);
    normalize_vector(a);
    v = v * cos(angle) + (a * (1 - cos(angle)) * a.dot_product(v)) + (a * v) * sin(angle);
}
struct PT getDirectionVector()
{
    struct PT veloceVect(cos(deg2rad(angleOfArrow)), sin(deg2rad(angleOfArrow)), 0);
    normalize_vector(veloceVect);
    return veloceVect;
}
void updateAngleOfArrow(struct PT direction_vector)
{
    double x = direction_vector.x;
    double y = direction_vector.y;
    if (x > 0 and y >= 0)
        angleOfArrow = atan(y / x) * 180 / M_PI;
    else if (x < 0 and y >= 0)
        angleOfArrow = 180 - (atan(y / -x) * 180 / M_PI);
    else if (x < 0 and y <= 0)
        angleOfArrow = 180 + (atan(-y / -x) * 180 / M_PI);
    else
        angleOfArrow = 360 - (atan(-y / x) * 180 / M_PI);
}
void axes()
{
    glBegin(GL_LINES);
    {
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(0, 0, 0);
        glVertex3f(100, 0, 0);

        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 100, 0);

        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, 100);
    }
    glEnd();
}

void drawArrow(GLdouble sector_size, GLdouble r_cy, GLdouble h_cy)
{

    GLdouble h_c = h_cy;
    GLdouble r_c = 1.5 * r_cy;
    GLdouble del_c = (2 * M_PI * r_c) / sector_size;
    GLdouble del_cy = (2 * M_PI * r_cy) / sector_size;
    GLdouble theta_c = del_c / r_c;
    GLdouble theta_cy = del_cy / r_cy;
    glPushMatrix();

    glTranslatef(sphere.x, sphere.y, sphere.z);
    glRotatef(angleOfArrow, 0, 0, 1);
    glTranslatef(sphere.radius, 0, 0);
    glRotatef(90, 0, 1, 0);
    glColor3f(0, 0, 1);

    for (int i = 0; i < sector_size; i++)
    {
        glPushMatrix();
        glRotated(i * 180 * theta_c / M_PI, 0, 0, 1);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0, 0, h_c + h_cy);
            glVertex3f(r_c, 0, h_cy);
            glVertex3f(r_c, del_c, h_cy);
        }
        glEnd();

        glPopMatrix();
        glPushMatrix();
        glRotated(i * 180 * theta_cy / M_PI, 0, 0, 1);
        glBegin(GL_QUADS);
        {
            glVertex3f(r_cy, 0, h_cy);
            glVertex3f(r_cy, 0, 0);
            glVertex3f(r_cy, del_cy, 0);
            glVertex3f(r_cy, del_cy, h_cy);
        }
        glEnd();
        glPopMatrix();
    }

    glPopMatrix();
}
void drawChecker(GLdouble radius, GLdouble wall_height, GLdouble wall_length)
{
    glPushMatrix();
    int color = 1;
    for (int i = -50; i < 50; i += 1)
    {
        for (int j = 50; j > -50; j -= 1)
        {
            glColor3f(color, color, color);
            glBegin(GL_QUADS);
            {
                glVertex3f(i, j, -radius);
                glVertex3f(i, j - 1, -radius);
                glVertex3f(i + 1, j - 1, -radius);
                glVertex3f(i + 1, j, -radius);
            }
            glEnd();
            color = !color;
        }
        color = !color;
    }
    for (int i = 0; i < 4; i++)
    {
        glRotatef(90 * i, 0, 0, 1);
        glColor3f(1, 0, 0);
        glBegin(GL_QUADS);
        {
            glVertex3f(wall_length, wall_length, -radius);
            glVertex3f(wall_length, wall_length, -radius + wall_height);
            glVertex3f(wall_length, -wall_length, -radius + wall_height);
            glVertex3f(wall_length, -wall_length, -radius);
        }
        glEnd();
    }
    glPopMatrix();
}
void drawSphere(float radius)
{
    if (sphere.ballInit)
    {
        float x, y, z, xy;                           // vertex position
        float nx, ny, nz, lengthInv = 1.0f / radius; // vertex normal
        float s, t;                                  // vertex texCoord

        float sectorStep = 2 * M_PI / sectorCount;
        float stackStep = M_PI / stackCount;
        float sectorAngle, stackAngle;

        for (int i = 0; i <= stackCount; ++i)
        {
            stackAngle = M_PI / 2 - i * stackStep; // starting from pi/2 to -pi/2
            xy = radius * cosf(stackAngle);        // r * cos(u)
            z = radius * sinf(stackAngle);         // r * sin(u)

            for (int j = 0; j <= sectorCount; ++j)
            {
                sectorAngle = j * sectorStep; // starting from 0 to 2pi

                // vertex position (x, y, z)
                x = xy * cosf(sectorAngle); // r * cos(u) * cos(v)
                y = xy * sinf(sectorAngle); // r * cos(u) * sin(v)

                points[i][j].x = x;
                points[i][j].y = y;
                points[i][j].z = z;
            }
        }
        sphere.ballInit = false;
    }

    glPushMatrix();

    glTranslatef(sphere.x, sphere.y, sphere.z);

    int colorCount = 0;
    for (int i = 0, color = 1; i < stackCount; i++)
    {
        for (int j = 0; j < sectorCount; j++)
        {
            if (j % 8 == 0 or j % 8 == 1 or j % 8 == 2 or j % 8 == 3)
            {
                if (i > stackCount / 2)
                    glColor3f(1.0, 0.0, 0.0);
                else
                    glColor3f(0.0, 1.0, 0.0);
            }
            else
            {
                if (i > stackCount / 2)
                    glColor3f(0.0, 1.0, 0.0);
                else
                    glColor3f(1.0, 0.0, 0.0);
            }
            colorCount++;
            if (colorCount > sectorCount)
                colorCount = 0;

            glBegin(GL_QUADS);
            {
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
            }
            glEnd();
        }
    }
    glPopMatrix();
}
void sphereMove()
{
    // sphere rotation + translation = sphere moving illusion
    // Rotation vector of the sphere (glRotatef at drawSphere() + glTranslatef at drawArrow())
    struct PT z_axis = PT(0, 0, 1);
    struct PT a = getDirectionVector();
    rotate_general(a, z_axis, 90);

    double degInc = DEG_ROT;
    if (sphere.goForward == false)
        degInc = -degInc;
    for (int i = 0; i <= stackCount; i++)
    {
        for (int j = 0; j <= sectorCount; j++)
        {
            struct PT temp = points[i][j];
            rotate_general(temp, a, degInc);
            points[i][j] = temp;
        }
    }
    
    struct PT vv = getDirectionVector();
    sphere.x += vv.x * sphere.radius * deg2rad(degInc);
    sphere.y += vv.y * sphere.radius * deg2rad(degInc);

    if (sphere.simulation == 0)
    {
        // Reflection after collision
        
        // Collision Detection x axis
        double x1_dist = sphere.walllength - sphere.x - sphere.radius;
        if (x1_dist <= 0.01)
        {
            vv.x *= -1;
            updateAngleOfArrow(vv);
        }
        double x2_dist = sphere.x + sphere.walllength - sphere.radius;
        if (x2_dist <= 0.01)
        {
            vv.x *= -1;
            updateAngleOfArrow(vv);
        }
        // Collision Dectection y axis
        double y1_dist = sphere.walllength - sphere.y - sphere.radius;
        if (y1_dist <= 0.01)
        {
            vv.y *= -1;
            updateAngleOfArrow(vv);
        }
        double y2_dist = sphere.y + sphere.walllength - sphere.radius;
        if (y2_dist <= 0.01)
        {
            vv.y *= -1;
            updateAngleOfArrow(vv);
        }

        angleOfArrow %= 360;
        sphere.direction_vector = getDirectionVector(); // Updating the direction vector of the sphere
    }
}
void display()
{
    glEnable(GL_DEPTH_TEST);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(
        cam.pos.x, cam.pos.y, cam.pos.z,
        cam.pos.x + cam.look.x, cam.pos.y + cam.look.y, cam.pos.z + cam.look.z,
        cam.up.x, cam.up.y, cam.up.z);

    axes();
    drawChecker(sphere.radius, sphere.wallheight, sphere.walllength);
    drawSphere(sphere.radius);
    drawArrow(100, 0.1, 0.4);

    // glFlush();

    glutSwapBuffers();
}
void calculateTimesForCollision();

void Timer(int value)
{
    // Schedule next collision
    
    if (sphere.simulation == true)
    {   
        
        struct PT vv = getDirectionVector();
        if (sphere.xcollision == true)
        {
            // if (vv.x > 0) {
            //     sphere.x = sphere.walllength - sphere.radius;
            // }
            // else if (vv.x < 0) sphere.x = -sphere.walllength + sphere.radius;
            
            vv.x *= -1;
            sphere.xcollision = false;
            updateAngleOfArrow(vv);
        }
        else if (sphere.ycollision == true)
        {
            // if (vv.y > 0) {
            //     sphere.y = sphere.walllength - sphere.radius;
            // }
            // else if (vv.y < 0.001) sphere.y = -sphere.walllength + sphere.radius;
            vv.y *= -1;
            sphere.ycollision = false;
            updateAngleOfArrow(vv);
        }
        calculateTimesForCollision();
    }
}
int inctime = 45;
void rollSphere(int value) {
    if (sphere.simulation == true) {
        sphere.goForward = true;
        sphereMove();
        glutTimerFunc(inctime, rollSphere, 0);
    } 
}

void calculateTimesForCollision() {
        // 0 - x+
        // 1 - y+
        // 2 - x-
        // 3 - y-
        struct PT vv = getDirectionVector();
        double timesOfCollisions[4];
        double xdistpos = (sphere.walllength - sphere.x - sphere.radius);
        double ydistpos = (sphere.walllength - sphere.y - sphere.radius);
        double xdistneg = (-sphere.walllength - sphere.x - sphere.radius);
        double ydistneg = (-sphere.walllength - sphere.y - sphere.radius);
        double perRotationDistance = sphere.radius * deg2rad(DEG_ROT);

        timesOfCollisions[0] = xdistpos*1.0/(vv.x * perRotationDistance);
        timesOfCollisions[1] = ydistpos*1.0/(vv.y * perRotationDistance);
        timesOfCollisions[2] = xdistneg*1.0/(vv.x * perRotationDistance);
        timesOfCollisions[3] = ydistneg*1.0/(vv.y * perRotationDistance);

        int savedIdx = 0;
        double minTime = INT_MAX;
        for (int i=0; i<4; i++) {
            if (timesOfCollisions[i] > 0) {
                if (minTime > timesOfCollisions[i]) {
                    minTime = timesOfCollisions[i];
                    savedIdx = i;
                }
            }
        }
        if (savedIdx % 2 == 0) sphere.xcollision = true;
        else sphere.ycollision = true;

        minTime = minTime * inctime;

        glutTimerFunc(minTime, Timer, 0);
}

void init()
{
    // glClearColor(0.1f, .0f, 0.0f, 1.0f); // Set background color to black and opaque

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50, 1, 1, 100);

    cam.pos = PT(10, 0, 0);
    cam.right = PT(0, 1, 0);
    cam.up = PT(0, 0, 1);
    cam.look = PT(-1, 0, 0);
}

void keyboardHandler(unsigned char key, int x, int y)
{
    double degInc = DEG_INC;
    switch (key)
    {
    case '1':
        rotation_orth(cam.look, cam.right, cam.up, deg2rad(degInc));
        break;

    case '2':
        rotation_orth(cam.look, cam.right, cam.up, -deg2rad(degInc));
        break;

    case '3':
        rotation_orth(cam.look, cam.up, cam.right, deg2rad(degInc));
        break;

    case '4':
        rotation_orth(cam.look, cam.up, cam.right, -deg2rad(degInc));
        break;

    case '5':
        rotation_orth(cam.right, cam.up, cam.look, deg2rad(degInc));
        break;

    case '6':
        rotation_orth(cam.right, cam.up, cam.look, -deg2rad(degInc));
        break;

    case 'j':
        angleOfArrow += DEG_INC;
        angleOfArrow = angleOfArrow % 360;
        break;

    case 'l':
        angleOfArrow -= DEG_INC;
        angleOfArrow = (angleOfArrow + 360) % 360;
        break;

    case 'i':
        sphere.goForward = true;
        sphereMove();
        break;

    case 'k':
        sphere.goForward = false;
        sphereMove();
        break;

    case ' ':
        sphere.simulation = !sphere.simulation;
        calculateTimesForCollision();
        glutTimerFunc(0, rollSphere, 0);
        break;

    case 27:
        exit(0);
        break;

    default:
        printf("ulta palta key press kora bondho koren!\n");
        break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP:
        cam.pos.x += cam.look.x;
        cam.pos.y += cam.look.y;
        cam.pos.z += cam.look.z;
        break;

    case GLUT_KEY_DOWN:
        cam.pos.x -= cam.look.x;
        cam.pos.y -= cam.look.y;
        cam.pos.z -= cam.look.z;
        break;

    case GLUT_KEY_RIGHT:
        cam.pos.x += cam.right.x;
        cam.pos.y += cam.right.y;
        cam.pos.z += cam.right.z;
        break;

    case GLUT_KEY_PAGE_UP:
        cam.pos.x += cam.up.x;
        cam.pos.y += cam.up.y;
        cam.pos.z += cam.up.z;
        break;

    case GLUT_KEY_LEFT:
        cam.pos.x -= cam.right.x;
        cam.pos.y -= cam.right.y;
        cam.pos.z -= cam.right.z;
        break;

    case GLUT_KEY_PAGE_DOWN:
        cam.pos.x -= cam.up.x;
        cam.pos.y -= cam.up.y;
        cam.pos.z -= cam.up.z;
        break;
    }
}

void idle()
{
    // printf("Currently we don't have any job!!!\n");

    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    printf("Hello World\n");
    glutInit(&argc, argv);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(780, 780);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("OpenGL Demo");

    init();

    glutDisplayFunc(display);

    glutKeyboardFunc(keyboardHandler);
    glutSpecialFunc(specialKeyListener);
    glutIdleFunc(idle);
    glutMainLoop();
    return 0;
}