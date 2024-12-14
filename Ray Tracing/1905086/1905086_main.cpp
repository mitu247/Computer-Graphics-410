#include "1905086_classes.h"
#include <fstream>
#include "bitmap_image.hpp"
std::vector<Object *> objects;
std::vector<Light *> pointLights;
std::vector<Light *> spotLights;

int recursionLevel, pixelDimension;
bitmap_image image = bitmap_image(780, 780);
int imageCount = 1;
const double degInc = 10;
double deg2rad(double deg)
{
    return (deg * M_PI) / 180;
}

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

void colorClamp(vector<double> &color) {
    for (int i = 0; i < 3; ++i) {
        if (color[i] > 1) {
            color[i] = 1;
        } else if (color[i] < 0) {
            color[i] = 0;
        }
    }
}

void capture() {
    std::cout << "Capturing image..." << std::endl;
    int windowWidth = 780;
    int windowHeight = 780;
    double viewAngle = 50 * (M_PI / 180); // Field of view of the camera in radians

    // Initialize bitmap image and set background color
    image.set_all_channels(0, 0, 0); // Set background to black

    double planeDistance = (windowHeight / 2.0) / tan(viewAngle / 2.0);
    PT topleft = cam.pos + cam.look*planeDistance - cam.right*windowWidth/2 + cam.up*windowHeight/2;
    double du = windowWidth / (double)image.width();
    double dv = windowHeight / (double)image.height();

    // Choose middle of the grid cell
    topleft = topleft + cam.right*(0.5*du) - cam.up*(0.5*dv);

    double tMin;
    Object *nearest;
    for (int i = 0; i < image.width(); ++i) {
        for (int j = 0; j < image.height(); ++j) {
            // Calculate curPixel using topleft, r, u, i, j, du, dv
            PT curPixel = topleft + cam.right*i*du - cam.up*j*dv;
            Ray ray(cam.pos, curPixel - cam.pos);
            vector<double> color = {0, 0, 0};
            double t = std::numeric_limits<double>::infinity();
            nearest = nullptr;
            for (Object *obj : objects) {
                double intersect = obj->intersect(ray, color, 0);
                if (intersect < t && intersect > 0) {
                    t = intersect;
                    nearest = obj;
                }
            }

            if (nearest != nullptr) {
                color[0] = color[1] = color[2] = 0;
                tMin = nearest->intersect(ray, color, 1);
                colorClamp(color);
                image.set_pixel(i, j, color[0]*255, color[1]*255, color[2]*255);
            }
        }
    }
    std::string filename = "Output_1" + std::to_string(imageCount++) + ".bmp";
    image.save_image(filename);
    std::cout << "Image saved as " << filename << std::endl;
}




void rotation(struct PT &v1, struct PT &v2, struct PT a, double an)
{
    v1 = v1 * cos(an) + (a * v1) * sin(an);
    v2 = v2 * cos(an) + (a * v2) * sin(an);
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
    case '0':
        capture();
        break;

    default:
        printf("ulta palta key press kora bondho koren!\n");
        break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    double factor = 2;
    switch (key)
    {
    case GLUT_KEY_UP:
        cam.pos.x += cam.look.x*factor;
        cam.pos.y += cam.look.y*factor;
        cam.pos.z += cam.look.z*factor;
        break;

    case GLUT_KEY_DOWN:
        cam.pos.x -= cam.look.x*factor;
        cam.pos.y -= cam.look.y*factor;
        cam.pos.z -= cam.look.z*factor;
        break;

    case GLUT_KEY_RIGHT:
        cam.pos.x += cam.right.x*factor;
        cam.pos.y += cam.right.y*factor;
        cam.pos.z += cam.right.z*factor;
        break;

    case GLUT_KEY_PAGE_UP:
        cam.pos.x += cam.up.x*factor;
        cam.pos.y += cam.up.y*factor;
        cam.pos.z += cam.up.z*factor;
        break;

    case GLUT_KEY_LEFT:
        cam.pos.x -= cam.right.x*factor;
        cam.pos.y -= cam.right.y*factor;
        cam.pos.z -= cam.right.z*factor;
        break;

    case GLUT_KEY_PAGE_DOWN:
        cam.pos.x -= cam.up.x*factor;
        cam.pos.y -= cam.up.y*factor;
        cam.pos.z -= cam.up.z*factor;
        break;
    }
}


void readFromFile()
{
    std::ifstream file("scene.txt");
    file >> recursionLevel >> pixelDimension;

    int numObjects;
    file >> numObjects;
    for (int i = 0; i < numObjects; ++i)
    {
        std::string type;
        file >> type;
        if (type == "sphere")
        {
            std::vector<double> details(12);
            for (int j = 0; j < 12; ++j)
            {
                file >> details[j];
            }
            PT center = PT(details[0], details[1], details[2]);
            double radius = details[3];
            Object *sphere = new Sphere(center, radius);
            sphere->setColor(details[4], details[5], details[6]);
            sphere->setCoEfficients(details[7], details[8], details[9], details[10]);
            sphere->setShine(details[11]);
            objects.push_back(sphere);
        }
        else if (type == "triangle")
        {
            std::vector<double> details(17);
            for (int j = 0; j < 17; ++j)
            {
                file >> details[j];
            }
            PT a = PT(details[0], details[1], details[2]);
            PT b = PT(details[3], details[4], details[5]);
            PT c = PT(details[6], details[7], details[8]);
            Object *triangle = new Triangle(a, b, c);
            triangle->setColor(details[9], details[10], details[11]);
            triangle->setCoEfficients(details[12], details[13], details[14], details[15]);
            triangle->setShine(details[16]);
            objects.push_back(triangle);
        }
        else if (type == "general")
        {
            std::vector<double> details(24);
            for (int j = 0; j < 24; ++j)
            {
                file >> details[j];
            }
            Object *general = new General(details[0], details[1], details[2], details[3], details[4], details[5], details[6], details[7], details[8], details[9],
                                          details[10], details[11], details[12], details[13], details[14], details[15]);
            general->setColor(details[16], details[17], details[18]);
            general->setCoEfficients(details[19], details[20], details[21], details[22]);
            general->setShine(details[23]);
            objects.push_back(general);
        }
    }

    int numPointLights;
    file >> numPointLights;
    cout << "numPointLights: " << numPointLights << endl;
    for (int i = 0; i < numPointLights; ++i)
    {
        Light *light = new PointLight;
        file >> light->light_pos.x >> light->light_pos.y >> light->light_pos.z;
        for (int j = 0; j < 3; ++j) file >> light->color[j];
        pointLights.push_back(light);
    }

    int numSpotLights;
    file >> numSpotLights;
    for (int i = 0; i < numSpotLights; ++i)
    {
        Light *light = new SpotLight;
        file >> light->light_pos.x >> light->light_pos.y >> light->light_pos.z;
        for (int j = 0; j < 3; ++j) file >> light->color[j];
        file >> light->light_direction.x >> light->light_direction.y >> light->light_direction.z;
        file >> light->cutoff_angle;
        spotLights.push_back(light);
    }

    Object *floor = new Floor(1000, 20);
    floor->setColor(1, 1, 1);
    floor->setCoEfficients(0.4, 0.4, 0.4, 0.2);
    floor->setShine(5);
    objects.push_back(floor);

    file.close();
}

void init()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50, 1, 1, 1000);

    cam.pos = PT(100, 100, 50);
    cam.right = PT(-1 / sqrt(2), 1 / sqrt(2), 0);
    cam.up = PT(0, 0, 1);
    cam.look = PT(-1 / sqrt(2), -1 / sqrt(2), 0);
    readFromFile();
}

void display()
{

    glEnable(GL_DEPTH_TEST);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // set up the camera
    gluLookAt(cam.pos.x, cam.pos.y, cam.pos.z,
              cam.pos.x + cam.look.x, cam.pos.y + cam.look.y, cam.pos.z + cam.look.z,
              cam.up.x, cam.up.y, cam.up.z);

    // draw objects here
    for (auto object : objects)
    {
        object->draw();
    }

    // swap buffers
    glutSwapBuffers();
}

void idle()
{
    // printf("Currently we don't have any job!!!\n");

    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowPosition(650, 0);
    glutInitWindowSize(780, 780);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("OpenGL Demo");

    init();

    glutDisplayFunc(display);

    glutKeyboardFunc(keyboardHandler);
    glutSpecialFunc(specialKeyListener);
    glutIdleFunc(idle);
    glutMainLoop();
    objects.clear();
    pointLights.clear();
    spotLights.clear();
    image.clear();
    return 0;
}