#include <bits/stdc++.h>
#define _USE_MATH_DEFINES
#ifdef __linux__
#include <GL/glut.h>
#include <GL/gl.h>
#elif WIN32
#include <GL/glut.h>
#include <GL/gl.h>
#endif
using namespace std;
GLdouble degInc = 5,initDeg = 0; // perclick degree increase


GLdouble deg2rad(GLdouble deg)
{
    return (deg * M_PI) / 180;
}

void axes()
{
    glLineWidth(3);
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
struct TSF
{
    GLdouble t, max_t;
    TSF()
    {
        t = 1.0;
        max_t = 1.0;
    }
} tsf;
struct PT
{
    GLdouble x, y, z;

    PT() {}
    PT(GLdouble x, GLdouble y, GLdouble z)
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
        GLdouble newx = this->x + p.x;
        GLdouble newy = this->y + p.y;
        GLdouble newz = this->z + p.z;

        return PT(newx, newy, newz);
    }

    // vector * number
    struct PT operator*(GLdouble num)
    {
        GLdouble newx = this->x * num;
        GLdouble newy = this->y * num;
        GLdouble newz = this->z * num;

        return PT(newx, newy, newz);
    }

    // vector / number
    struct PT operator/(GLdouble num)
    {
        GLdouble newx = this->x / num;
        GLdouble newy = this->y / num;
        GLdouble newz = this->z / num;

        return PT(newx, newy, newz);
    }

    // vector x vector
    struct PT operator*(const struct PT &p)
    {
        GLdouble newx = (this->y * p.z) - (p.y * this->z);
        GLdouble newy = (p.x * this->z) - (this->x * p.z);
        GLdouble newz = (this->x * p.y) - (p.x * this->y);

        return PT(newx, newy, newz);
    }
};

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

struct SphereInfo
{
    struct PT points[500][500];
    GLdouble radius;
    GLdouble maxRadius;

    SphereInfo() {}

    SphereInfo(GLdouble radius, GLdouble maxRadius)
    {
        this->radius = radius;
        this->maxRadius = maxRadius;
    }
} sp;

struct CylinderInfo
{
    struct PT points[500][500];
    GLdouble radius;
    GLdouble height, maxHeight;
    GLdouble diAngle;

    CylinderInfo() {}
    CylinderInfo(GLdouble radius, GLdouble maxHeight, GLdouble diAngle)
    {
        this->radius = radius;
        this->maxHeight = maxHeight;
        this->height = maxHeight;
        this->diAngle = diAngle;
    }
} cy;

void triangle()
{
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(1, 0, 0);
        glVertex3f(0, 1, 0);
        glVertex3f(0, 0, 1);
    }
    glEnd();
}

void rotation(struct PT &v1, struct PT &v2, struct PT a, GLdouble an)
{
    v1 = v1 * cos(an) + (a * v1) * sin(an);
    v2 = v2 * cos(an) + (a * v2) * sin(an);
}
void generatePoints(int sectorCount)
{
    GLdouble st = -1, en = 1;
    for (int i = 0; i <= sectorCount; i++)
    {
        GLdouble x = st + ((en - st) / sectorCount) * i;
        for (int j = 0; j <= sectorCount; j++)
        {
            GLdouble y = st + ((en - st) / sectorCount) * j;
            sp.points[i][j].x = x;
            sp.points[i][j].y = y;
            sp.points[i][j].z = 1;
        }
    }
}
void getCurvedSurfacePoints(int sectorCount)
{
    GLdouble unit;
    for (int i = 0; i <= sectorCount; i++)
    {
        for (int j = 0; j <= sectorCount; j++)
        {
            unit = sqrt((sp.points[i][j].x * sp.points[i][j].x) + (sp.points[i][j].y * sp.points[i][j].y) + (sp.points[i][j].z * sp.points[i][j].z));
            sp.points[i][j] = sp.points[i][j] / unit;
            sp.points[i][j] = sp.points[i][j] * sp.radius;
        }
    }
    for (int i = 0; i < sectorCount; i++)
    {
        for (int j = 0; j < sectorCount; j++)
        {
            glBegin(GL_QUADS);
            {
                glVertex3f(sp.points[i][j].x, sp.points[i][j].y, sp.points[i][j].z);
                glVertex3f(sp.points[i + 1][j].x, sp.points[i + 1][j].y, sp.points[i + 1][j].z);
                glVertex3f(sp.points[i + 1][j + 1].x, sp.points[i + 1][j + 1].y, sp.points[i + 1][j + 1].z);
                glVertex3f(sp.points[i][j + 1].x, sp.points[i][j + 1].y, sp.points[i][j + 1].z);
            }
            glEnd();
        }
    }
}
void cylinder()
{
    glColor3f(1, 1, 0);
    GLdouble sectors = 200;
    GLdouble cylinderGridWidth = 2 * tan(deg2rad(cy.diAngle / 2));
    cy.height = sqrt(2) * tsf.t;
    cy.radius = sp.radius;


    for (int i = 0; i <= sectors; i++)
    {
        GLdouble z = -cy.height/2 + (cy.height / sectors) * i;
        for (int j = 0; j <= sectors; j++)
        {
            GLdouble x = -tan(deg2rad(cy.diAngle / 2)) + (cylinderGridWidth / sectors) * j;
            GLdouble y = cy.radius;

            cy.points[i][j].x = (cy.radius * x) / sqrt(x * x + y * y);
            cy.points[i][j].y = (cy.radius * y) / sqrt(x * x + y * y);
            cy.points[i][j].z = z;
        }
    }
    for (int i = 0; i < sectors; i++)
    {
        for (int j = 0; j < sectors; j++)
        {
            glBegin(GL_QUADS);
            {
                glVertex3f(cy.points[i][j].x, cy.points[i][j].y, cy.points[i][j].z);
                glVertex3f(cy.points[i + 1][j].x, cy.points[i + 1][j].y, cy.points[i + 1][j].z);
                glVertex3f(cy.points[i + 1][j + 1].x, cy.points[i + 1][j + 1].y, cy.points[i + 1][j + 1].z);
                glVertex3f(cy.points[i][j + 1].x, cy.points[i][j + 1].y, cy.points[i][j + 1].z);
            }
            glEnd();
        }
    }

}
void placeCylinders() {
        // Placing the cylinders in their positions
    glPushMatrix();

    for (int j=0; j<2; j++) {
        // Rotate the whole upper pyramid
        glRotatef(180*j, 0, 1, 0);
        for (int i=0; i<4; i++) {
            glPushMatrix();
            glRotatef(90*i, 0, 0, 1);
            glRotatef(45, 1, 0, 0);
            glTranslatef(0, tsf.t/sqrt(2), 0);
            cylinder();
            glPopMatrix();
        }
    }
    
    glPopMatrix();

    for (int i=0; i<4; i++) {
        glPushMatrix();
        glRotatef(45 + 90 * i, 0, 0, 1);
        glRotatef(90, 0, 1, 0);
        glTranslatef(0, tsf.t/sqrt(2), 0);
        cylinder();
        glPopMatrix();
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
    glPushMatrix();
    glRotatef(initDeg,0,0,1);

    glPushMatrix();
    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        glColor3f(i % 2, (i + 1) % 2, 1);
        glRotatef(90.0 * i, 1, 0, 0);
        glTranslatef((tsf.max_t - tsf.t) / 3, (tsf.max_t - tsf.t) / 3, (tsf.max_t - tsf.t) / 3);
        glScalef(tsf.t, tsf.t, tsf.t);
        triangle();
        glPopMatrix();
    }

    glRotatef(-90, 0, 1, 0);

    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        glColor3f((i + 1) % 2, (i) % 2, 0);
        glRotatef(90.0 * i, 0, 0, 1);
        glTranslatef((tsf.max_t - tsf.t) / 3, (tsf.max_t - tsf.t) / 3, (tsf.max_t - tsf.t) / 3);
        glScalef(tsf.t, tsf.t, tsf.t);
        triangle();
        glPopMatrix();
    }

   glPopMatrix();

    // sphere cubes
    for (int i = 0; i < 4; i++)
    {
        glPushMatrix();
        glColor3f((i % 2), (i + 1) % 2, 0);
        glRotatef(90 * i, 1, 0, 0);
        glTranslatef(0, 0, tsf.t);
        generatePoints(200);
        getCurvedSurfacePoints(200);
        glPopMatrix();
    }
    for (int i = 0; i < 2; i++)
    {
        glPushMatrix();
        glColor3f(0, 0, 1);
        glRotatef(90 + 180 * i, 0, 1, 0);
        glTranslatef(0, 0, tsf.t);
        generatePoints(200);
        getCurvedSurfacePoints(200);
        glPopMatrix();
    }

    // cylinders
    placeCylinders();
    
    glPopMatrix();

    glutSwapBuffers();
    glFlush();
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

    sp.radius = 0;
    sp.maxRadius = tsf.t / sqrt(3);

    cy.radius = sp.radius;
    cy.maxHeight = tsf.max_t*sqrt(2);
    cy.diAngle = 70.5287794;
    cy.height = cy.maxHeight;
}

void keyboardHandler(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '1':
        rotation(cam.look, cam.right, cam.up, deg2rad(degInc));
        break;

    case '2':
        rotation(cam.look, cam.right, cam.up, -deg2rad(degInc));
        break;

    case '3':
        rotation(cam.look, cam.up, cam.right, deg2rad(degInc));
        break;

    case '4':
        rotation(cam.look, cam.up, cam.right, -deg2rad(degInc));
        break;

    case '5':
        rotation(cam.right, cam.up, cam.look, deg2rad(degInc));
        break;

    case '6':
        rotation(cam.right, cam.up, cam.look, -deg2rad(degInc));
        break;

    case ',':
        tsf.t = tsf.t - 0.1;
        sp.radius += (sp.maxRadius / 10);
        if (tsf.t < 0)
            tsf.t = 0, sp.radius = sp.maxRadius;
        break;

    case '.':
        tsf.t = tsf.t + 0.1;
        sp.radius -= (sp.maxRadius / 10);
        if (tsf.t > tsf.max_t)
            tsf.t = tsf.max_t, sp.radius = 0;
        break;
    case 'a':
        initDeg-=degInc;
        break;
    case 'd':
        initDeg+=degInc;
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
    //printf("Currently we don't have any job!!!\n");

    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    printf("Hello World\n");
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
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