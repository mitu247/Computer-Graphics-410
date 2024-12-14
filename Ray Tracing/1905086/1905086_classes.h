#ifndef _1905086_CLASSES_H_
#define _1905086_CLASSES_H_

#include <bits/stdc++.h>
#include <GL/glut.h>
using namespace std;

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

    PT operator=(const struct PT &p)
    {
        this->x = p.x;
        this->y = p.y;
        this->z = p.z;

        return *this;
    }

    // vector + vector
    PT operator+(const PT &p)
    {
        return PT(this->x + p.x, this->y + p.y, this->z + p.z);
    }

    // vector - vector
    PT operator-(const PT &p)
    {
        return PT(this->x - p.x, this->y - p.y, this->z - p.z);
    }

    // vector * number
    PT operator*(const double num)
    {
        return PT(this->x * num, this->y * num, this->z * num);
    }

    // vector / number
    PT operator/(const double num)
    {
        return PT(this->x / num, this->y / num, this->z / num);
    }

    // vector x vector
    PT operator*(const struct PT &p)
    {
        double newx = (this->y * p.z) - (p.y * this->z);
        double newy = (p.x * this->z) - (this->x * p.z);
        double newz = (this->x * p.y) - (p.x * this->y);

        return PT(newx, newy, newz);
    }
    // dot product
    double dot(const PT &p)
    {
        return (this->x * p.x) + (this->y * p.y) + (this->z * p.z);
    }
    // normalize
    PT normalize()
    {
        double mag = sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
        return PT(this->x / mag, this->y / mag, this->z / mag);
    }
};
class Object;
class Light;
extern std::vector<Object *> objects;
extern std::vector<Light *> pointLights;
extern std::vector<Light *> spotLights;
extern int recursionLevel, pixelDimension;

class Ray
{
public:
    PT origin;
    PT dir; // normalize for easier calculations

    Ray() {}

    Ray(PT origin, PT dir)
    {
        this->origin = origin;
        this->dir = dir.normalize();
    }
};

class Light
{
public:
    PT light_pos;
    double color[3];
    PT light_direction;
    double cutoff_angle;
    Light() {}
};
class PointLight : public Light
{
public:
    PointLight() {}
};

class SpotLight : public Light
{
public:
    SpotLight() {}
};

class Object
{
protected:
    PT reference_point;
    double height, width, length;
    vector<double> color = {0, 0, 0};
    double coEfficients[4] = {0, 0, 0, 0};
    int shine;

    void letThereBeLight(Ray r, PT intersection, vector<double> &color, vector<double> objectColor, vector<Light *> lights, int idx)
    {
        PT normal = getNormal(intersection, r);

        PT lightDir = intersection - lights[idx]->light_pos;
        double lightDistance = sqrt(lightDir.dot(lightDir));
        Ray lightToObjRay = Ray(lights[idx]->light_pos, lightDir);
        for (int j = 0; j < objects.size(); j++)
        {
            double t = objects[j]->getObjectParameter(lightToObjRay, color, 0);
            if (t > 0 && (lightDistance - t) > 0.00001)
            {
                return; // Object is obscured
            }
        }

        // Light to object direction
        lightDir = lightDir.normalize();
        lightDir = lightDir * -1;

        double lambertValue = normal.dot(lightDir);
        if (lambertValue > 0)
        {
            color[0] += objectColor[0] * coEfficients[1] * lights[idx]->color[0] * lambertValue;
            color[1] += objectColor[1] * coEfficients[1] * lights[idx]->color[1] * lambertValue;
            color[2] += objectColor[2] * coEfficients[1] * lights[idx]->color[2] * lambertValue;
        }
        // Phong value
        PT view = r.origin - intersection;
        view = view.normalize();
        // lightDir is considered now as Light->Object direction
        PT reflection = lightDir - normal * (2 * normal.dot(lightDir));
        reflection = reflection.normalize();
        double specularValue = std::pow(reflection.dot(view.normalize()), shine);
        if (specularValue > 0)
        {
            color[0] += objectColor[0] * coEfficients[2] * lights[idx]->color[0] * specularValue;
            color[1] += objectColor[1] * coEfficients[2] * lights[idx]->color[1] * specularValue;
            color[2] += objectColor[2] * coEfficients[2] * lights[idx]->color[2] * specularValue;
        }

        return; // Let there be light in the object
    }

public:
    Object() {}
    virtual void draw() {}
    virtual double getObjectParameter(Ray, vector<double> &, int) = 0;
    // getColoratIntersectionPoint will return color array
    virtual vector<double> getColorAtIntersectionPoint(PT point) = 0;
    virtual PT getNormal(PT, Ray) = 0;
    virtual double intersect(Ray r, vector<double> &color, int level)
    {
        double tMin = getObjectParameter(r, color, 0);
        if (tMin < 0)
            return -1.0;
        if (level == 0)
            return tMin;

        PT intersection = r.origin + r.dir * tMin;
        vector<double> objectColor = getColorAtIntersectionPoint(intersection);

        color[0] = objectColor[0] * coEfficients[0];
        color[1] = objectColor[1] * coEfficients[0];
        color[2] = objectColor[2] * coEfficients[0];

        // For Point Lights
        for (int i = 0; i < pointLights.size(); i++)
        {
            letThereBeLight(r, intersection, color, objectColor, pointLights, i);
        }

        // For Spot Lights
        for (int i = 0; i < spotLights.size(); i++)
        {
            Ray spotLightRay(spotLights[i]->light_pos, intersection - spotLights[i]->light_pos); // Light Ray

            double beta = acos(spotLightRay.dir.dot(spotLights[i]->light_direction.normalize())) * 180 / M_PI;

            // The object is in the hind side of the spotlight
            if (spotLights[i]->cutoff_angle < beta)
                continue;

            letThereBeLight(r, intersection, color, objectColor, spotLights, i);
        }

        if (recursionLevel <= level)
        {
            return tMin;
        }
        // Multilevel Reflection
        PT normal = getNormal(intersection, r);
        PT reflectionRayDir = r.dir - normal * (2 * normal.dot(r.dir));
        PT reflectionRayOrigin = intersection + reflectionRayDir * 0.000001; // to avoid self intersection
        Ray reflectionRay(reflectionRayOrigin, reflectionRayDir);

        Object *nearest = nullptr;
        double tMin2 = 1e10;
        for (int i = 0; i < objects.size(); i++)
        {
            double t = objects[i]->getObjectParameter(reflectionRay, color, 0);
            if (t > 0 && (tMin2 - t) > 0.00001)
            {
                tMin2 = t;
                nearest = objects[i];
            }
        }
        if (nearest != nullptr)
        {
            vector<double> color2 = {0, 0, 0};
            double t = nearest->intersect(reflectionRay, color2, level + 1);

            {
                color[0] += color2[0] * coEfficients[3];
                color[1] += color2[1] * coEfficients[3];
                color[2] += color2[2] * coEfficients[3];
            }
        }

        return tMin;
    }
    void setColor(double r, double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    void setShine(int s)
    {
        shine = s;
    }
    void setCoEfficients(double a, double d, double s, double r)
    {
        coEfficients[0] = a;
        coEfficients[1] = d;
        coEfficients[2] = s;
        coEfficients[3] = r;
    }
};

class Sphere : public Object
{
public:
    PT points[500][500];
    const int sectorCount = 100;
    const int stackCount = 100;

    Sphere(PT center, double radius)
    {
        reference_point = center;
        length = radius;
    }

    PT getNormal(PT point, Ray r) override
    {
        PT normal = point - reference_point;
        return normal.normalize();
    }

    void sphere(double rad)
    {
        double posX, posY, posZ, posXY; // vertex position
        double lengthInv = 1.0 / rad;   // vertex normal

        double sectorStep = 2 * M_PI / sectorCount;
        double stackStep = M_PI / sectorCount;
        double sectorAng, stackAng;

        for (int i = 0; i <= sectorCount; ++i)
        {
            stackAng = M_PI / 2 - i * stackStep; // starting from pi/2 to -pi/2
            posXY = rad * cos(stackAng);         // r * cos(u)
            posZ = rad * sin(stackAng);          // r * sin(u)

            for (int j = 0; j <= sectorCount; ++j)
            {
                sectorAng = j * sectorStep; // starting from 0 to 2pi

                // vertex position (x, y, z)
                posX = posXY * cos(sectorAng); // r * cos(u) * cos(v)
                posY = posXY * sin(sectorAng); // r * cos(u) * sin(v)

                points[i][j].x = posX;
                points[i][j].y = posY;
                points[i][j].z = posZ;
            }
        }
        glPushMatrix();

        for (int i = 0, color = 1; i < stackCount; i++)
        {
            for (int j = 0; j < sectorCount; j++)
            {
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

    void draw() override
    {
        glColor3f(color[0], color[1], color[2]);
        // Draw the sphere slices along the y-axis
        for (int i = 0; i < 2; i++)
        {
            glPushMatrix();
            glTranslatef(reference_point.x, reference_point.y, reference_point.z);
            glRotatef(i * 180 + 90, 0, 1, 0); // Rotate around y-axis
            sphere(length);
            glPopMatrix();
        }

        // Draw the sphere slices along the x-axis
        for (int i = 0; i < 4; i++)
        {
            glPushMatrix();
            glTranslatef(reference_point.x, reference_point.y, reference_point.z);
            glRotatef(i * 90, 1, 0, 0); // Rotate around x-axis

            sphere(length);
            glPopMatrix();
        }
    }
    double getObjectParameter(Ray r, vector<double> &color, int level) override
    {
        PT rayOriginToRefPoint = reference_point - r.origin;
        double projectionLength = rayOriginToRefPoint.dot(r.dir);

        if (projectionLength < 0)
            return -1.0;

        double perpendicularDistSquared = rayOriginToRefPoint.dot(rayOriginToRefPoint) - projectionLength * projectionLength;
        double radiusSquared = length * length;

        if (perpendicularDistSquared > radiusSquared)
            return -1.0;

        double halfChordLength = sqrt(radiusSquared - perpendicularDistSquared);
        double intersection1 = projectionLength - halfChordLength;
        double intersection2 = projectionLength + halfChordLength;

        if (intersection1 > intersection2)
            std::swap(intersection1, intersection2);

        if (intersection1 < 0)
        {
            intersection1 = intersection2;
            if (intersection1 < 0)
                return -1.0;
        }

        return intersection1;
    }

    vector<double> getColorAtIntersectionPoint(PT point) override
    {
        return color;
    }
};

class Floor : public Object
{
public:
    Floor(double floorWidth, double tileWidth)
    {
        reference_point = PT(-floorWidth / 2, -floorWidth / 2, 0);
        width = floorWidth;
        length = tileWidth;
    }

    // Normal of the floor
    PT getNormal(PT point, Ray r) override
    {
        return PT(0, 0, 1);
    }

    void draw() override
    {
        double width = -reference_point.x * 2;
        int tilesPerSide = width / length;

        int i = 0;
        int j = 0;

        while (i < tilesPerSide)
        {
            while (j < tilesPerSide)
            {
                bool isAlternateTile = (i + j) % 2 == 0;
                if (isAlternateTile)
                {
                    glColor3f(color[0], color[1], color[2]); // white
                }
                else
                {
                    glColor3f(0, 0, 0); // black
                }

                // Calculate tile corners
                double x1 = reference_point.x + i * length;
                double y1 = reference_point.y + j * length;
                double x2 = reference_point.x + (i + 1) * length;
                double y2 = reference_point.y + (j + 1) * length;

                // Draw tile
                glBegin(GL_QUADS);
                glVertex3d(x1, y1, reference_point.z);
                glVertex3d(x1, y2, reference_point.z);
                glVertex3d(x2, y2, reference_point.z);
                glVertex3d(x2, y1, reference_point.z);
                glEnd();
                j++;
            }

            j = 0;
            i++;
        }
    }

    double getObjectParameter(Ray r, vector<double> &color, int level) override
    {
        if (fabs(r.dir.z) < 1e-6)
        {
            return -1.0; // Ray is parallel to the floor
        }

        PT n = this->getNormal(r.origin, r);
        double t = -n.dot(r.origin - reference_point) / n.dot(r.dir);

        PT pointIntThePlane = r.origin + r.dir * t;

        // Check if the point is within the boundaries
        if (pointIntThePlane.x < reference_point.x || pointIntThePlane.x > -reference_point.x || pointIntThePlane.y < reference_point.y || pointIntThePlane.y > -reference_point.y)
        {
            return -1.0;
        }

        return t;
    }
    vector<double> getColorAtIntersectionPoint(PT point) override
    {
        int xCount = (point.x - reference_point.x) / length;
        int yCount = (point.y - reference_point.y) / length;
        int tilesPerSide = width / length;

        if (xCount < 0 || xCount >= tilesPerSide || yCount < 0 || yCount >= tilesPerSide)
        {
            return {0, 0, 0};
        }

        if ((xCount + yCount) % 2 == 0)
        {
            return {color[0], color[1], color[2]};
        }
        else
        {
            return {0, 0, 0}; // black
        }
    }
};

class Triangle : public Object
{
public:
    PT pointA, pointB, pointC;

    Triangle(PT a, PT b, PT c)
    {
        pointA = a;
        pointB = b;
        pointC = c;
    }

    PT getNormal(PT point, Ray r) override
    {
        PT v1 = pointB - pointA;
        PT v2 = pointC - pointA;
        PT normal = v1 * v2;
        normal = normal.normalize();

        if (normal.dot(r.dir) < 0)
            normal = normal * -1;

        return normal;
    }

    void draw()
    {
        glBegin(GL_TRIANGLES);
        glColor3f(color[0], color[1], color[2]);
        glVertex3d(pointA.x, pointA.y, pointA.z);
        glVertex3d(pointB.x, pointB.y, pointB.z);
        glVertex3d(pointC.x, pointC.y, pointC.z);
        glEnd();
    }

    bool isAnticlockwise(PT pointA, PT pointB, PT pointC)
    {
        // Calculate the cross product of the vectors AB and AC
        PT crossProduct = (pointB - pointA) * (pointC - pointA);

        // If the z-component of the cross product is positive, the points are in anticlockwise order
        return crossProduct.z > 0;
    }

    double getObjectParameter(Ray r, vector<double> &color, int level) override
    {
        // Ensure points are in anticlockwise order
        if (!isAnticlockwise(pointA, pointB, pointC))
        {
            std::swap(pointB, pointC);
        }

        PT firstEdge = pointB - pointA;
        PT secondEdge = pointC - pointA;
        PT crossProduct = r.dir * secondEdge;
        double dotProduct = firstEdge.dot(crossProduct);

        if (fabs(dotProduct) < 0.00001)
            return -1.0;

        double reciprocal = 1.0 / dotProduct;
        PT rayOriginVector = r.origin - pointA;
        double parameterU = reciprocal * rayOriginVector.dot(crossProduct);

        if (parameterU < 0.0 || parameterU > 1.0)
            return -1.0;

        PT crossProduct2 = rayOriginVector * firstEdge;
        double parameterV = reciprocal * (r.dir.dot(crossProduct2));

        if (parameterV < 0.0 || parameterU + parameterV > 1.0)
            return -1.0;

        double parameterT = reciprocal * (secondEdge.dot(crossProduct2));

        return (parameterT > 0.00001) ? parameterT : -1.0;
    }
    vector<double> getColorAtIntersectionPoint(PT point) override
    {
        return color;
    }
};

class General : public Object
{
private:
    bool cutDim(string dim)
    {
        if (dim == "length" and length > 0)
            return true;
        if (dim == "width" and width > 0)
            return true;
        if (dim == "height" and height > 0)
            return true;
        return false;
    }

    double clip(double min_t, double max_t, Ray r)
    {
        if (min_t < 0 and max_t < 0)
            return -1.0;
        if (min_t > 0)
        {
            PT newPoint = r.origin + r.dir * min_t;
            if (cutDim("length"))
            {
                if (newPoint.x > reference_point.x and newPoint.x < reference_point.x + length)
                    return min_t; // Check if the point is within the boundaries
            }
            if (cutDim("width"))
            {
                if (newPoint.y > reference_point.y and newPoint.y < reference_point.y + width)
                    return min_t; // Check if the point is within the boundaries
            }
            if (cutDim("height"))
            {
                if (newPoint.z > reference_point.z and newPoint.z < reference_point.z + height)
                    return min_t; // Check if the point is within the boundaries
            }
        }
        if (max_t > 0)
        {
            PT newPoint = r.origin + r.dir * max_t;
            if (cutDim("length"))
            {
                if (newPoint.x > reference_point.x and newPoint.x < reference_point.x + length)
                    return max_t; // Check if the point is within the boundaries
            }
            if (cutDim("width"))
            {
                if (newPoint.y > reference_point.y and newPoint.y < reference_point.y + width)
                    return max_t; // Check if the point is within the boundaries
            }
            if (cutDim("height"))
            {
                if (newPoint.z > reference_point.z and newPoint.z < reference_point.z + height)
                    return max_t; // Check if the point is within the boundaries
            }
        }
        return -1.0;
    }

public:
    double A, B, C, D, E, F, G, H, I, J;

    General(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j,
            double x, double y, double z, double l, double w, double p)
    {
        A = a;
        B = b;
        C = c;
        D = d;
        E = e;
        F = f;
        G = g;
        H = h;
        I = i;
        J = j;
        reference_point = PT(x, y, z);
        length = l;
        width = w;
        height = p;
    }

    double myPow(double base, int exp)
    {
        double result = 1.0;
        for (int i = 0; i < exp; i++)
        {
            result *= base;
        }
        return result;
    }

    double mySqrt(double value)
    {
        double x = value;
        double y = 1;
        double e = 0.000001; // precision
        while (x - y > e)
        {
            x = (x + y) / 2;
            y = value / x;
        }
        return x;
    }
    double myMin(double a, double b)
    {
        return (a < b) ? a : b;
    }
    double myMax(double a, double b)
    {
        return (a > b) ? a : b;
    }

    void draw() override
    {
        // Drawing a general quadric surface is complex and depends on the specific coefficients.
        // You might need to use a library or write a custom algorithm to draw it.
    }
    double getObjectParameter(Ray r, vector<double> &color, int level) override
    {
        // Coefficients for the quadratic equation
        double a = A * myPow(r.dir.x, 2) + B * myPow(r.dir.y, 2) + C * myPow(r.dir.z, 2) + D * r.dir.x * r.dir.y + E * r.dir.x * r.dir.z + F * r.dir.y * r.dir.z;
        double b = 2 * (A * r.origin.x * r.dir.x + B * r.origin.y * r.dir.y + C * r.origin.z * r.dir.z) + D * (r.origin.x * r.dir.y + r.dir.x * r.origin.y) + E * (r.origin.x * r.dir.z + r.dir.x * r.origin.z) + F * (r.origin.y * r.dir.z + r.dir.y * r.origin.z) + G * r.dir.x + H * r.dir.y + I * r.dir.z;
        double c = A * myPow(r.origin.x, 2) + B * myPow(r.origin.y, 2) + C * myPow(r.origin.z, 2) + D * r.origin.x * r.origin.y + E * r.origin.x * r.origin.z + F * r.origin.y * r.origin.z + G * r.origin.x + H * r.origin.y + I * r.origin.z + J;

        if (fabs(a) < 0.00001)
        {
            return -c / b; // The ray is parallel to the surface
        }
        // Calculate the discriminant
        double discriminant = myPow(b, 2) - 4 * a * c;

        // If the discriminant is negative, there are no real roots
        if (discriminant < 0)
        {
            return -1.0;
        }

        // Calculate the two roots of the quadratic equation
        double t0 = (-b - mySqrt(discriminant)) / (2 * a);
        double t1 = (-b + mySqrt(discriminant)) / (2 * a);

        // t0 would be the smallest one
        double min_t = myMin(t0, t1);
        double max_t = myMax(t0, t1);

        // Return the smallest positive root
        return clip(min_t, max_t, r);
    }
    vector<double> getColorAtIntersectionPoint(PT point) override
    {
        return color;
    }
    PT getNormal(PT point, Ray r) override
    {
        double x = 2 * A * point.x + D * point.y + E * point.z + G;
        double y = 2 * B * point.y + D * point.x + F * point.z + H;
        double z = 2 * C * point.z + E * point.x + F * point.y + I;
        return PT(x, y, z).normalize();
    }
};

#endif