#include <bits/stdc++.h>
#include <stack>
#include <cmath>
#include "1905086_z_buffer.cpp"
using namespace std;

fstream fin("scene.txt", ios::in);
fstream f_out("stage1.txt", ios::out);
fstream f_out2("stage2.txt", ios::out);
fstream f_out3("stage3.txt", ios::out);

struct Matrix4f
{
    float m[4][4];
};

struct Vector4f
{
    float v[4];
};

Vector4f operator*(float c, const Vector4f &x)
{
    Vector4f result;
    for (int i = 0; i < 4; ++i)
    {
        result.v[i] = c * x.v[i];
    }
    return result;
}
Vector4f operator+(const Vector4f &x, const Vector4f &y)
{
    Vector4f result;
    for (int i = 0; i < 4; ++i)
    {
        result.v[i] = x.v[i] + y.v[i];
    }
    return result;
}

Matrix4f IdentityMatrix()
{
    Matrix4f I = {1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1};
    return I;
}

Matrix4f multiplyMatrix(const Matrix4f &A, const Matrix4f &B)
{
    Matrix4f result;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            result.m[i][j] = 0;
            for (int k = 0; k < 4; ++k)
            {
                result.m[i][j] += A.m[i][k] * B.m[k][j];
            }
        }
    }
    return result;
}
Matrix4f generateTranslationMatrix(float tx, float ty, float tz)
{
    Matrix4f T = {1, 0, 0, tx,
                  0, 1, 0, ty,
                  0, 0, 1, tz,
                  0, 0, 0, 1};
    return T;
}

Matrix4f generateScalingMatrix(float sx, float sy, float sz)
{
    Matrix4f S = {sx, 0, 0, 0,
                  0, sy, 0, 0,
                  0, 0, sz, 0,
                  0, 0, 0, 1};
    return S;
}

Vector4f crossProduct(const Vector4f &a, const Vector4f &b)
{
    Vector4f result;
    result.v[0] = a.v[1] * b.v[2] - a.v[2] * b.v[1];
    result.v[1] = a.v[2] * b.v[0] - a.v[0] * b.v[2];
    result.v[2] = a.v[0] * b.v[1] - a.v[1] * b.v[0];
    result.v[3] = 0;
    return result;
}

float dotProduct(const Vector4f &a, const Vector4f &b)
{
    float result = 0;
    for (int i = 0; i < 3; ++i)
    {
        result += a.v[i] * b.v[i];
    }
    return result;
}

Vector4f R(const Vector4f &x, const Vector4f &a, float angle)
{
    float c = cos(angle * M_PI / 180.0);
    float s = sin(angle * M_PI / 180.0);
    return c * x + (1 - c) * dotProduct(a, x) * a + s * crossProduct(a, x);
}

Matrix4f generateRotationMatrix(float angle, float ax, float ay, float az)
{
    Vector4f a = {ax, ay, az, 0};
    float norm = sqrt(dotProduct(a, a));
    for (int i = 0; i < 3; ++i)
    {
        a.v[i] /= norm;
    }

    Vector4f i = {1, 0, 0, 0};
    Vector4f j = {0, 1, 0, 0};
    Vector4f k = {0, 0, 1, 0};

    Vector4f c1 = R(i, a, angle);
    Vector4f c2 = R(j, a, angle);
    Vector4f c3 = R(k, a, angle);

    Matrix4f R = {c1.v[0], c2.v[0], c3.v[0], 0,
                  c1.v[1], c2.v[1], c3.v[1], 0,
                  c1.v[2], c2.v[2], c3.v[2], 0,
                  0, 0, 0, 1};
    return R;
}
Vector4f transformPoint(const Matrix4f &M, const Vector4f &P)
{
    Vector4f P_prime;
    for (int i = 0; i < 4; ++i)
    {
        P_prime.v[i] = 0;
        for (int j = 0; j < 4; ++j)
        {
            P_prime.v[i] += M.m[i][j] * P.v[j];
        }
    }
    for (int i = 0; i < 4; ++i)
    {
        P_prime.v[i] /= P_prime.v[3];
    }
    return P_prime;
}

Matrix4f generateProjectionMatrix(float fovY, float aspectRatio, float near, float far)
{
    float fovX = fovY * aspectRatio;
    float t = near * tan((fovY * M_PI / 180.0) / 2);
    float r = near * tan((fovX * M_PI / 180.0) / 2);

    Matrix4f P = {near / r, 0, 0, 0,
                  0, near / t, 0, 0,
                  0, 0, -(far + near) / (far - near), -(2 * far * near) / (far - near),
                  0, 0, -1, 0};
    return P;
}

int main()
{
    stack<Matrix4f> S;
    Matrix4f M = IdentityMatrix();
    Vector4f eye, look, up, perspective;
    fin >> eye.v[0] >> eye.v[1] >> eye.v[2];
    eye.v[3] = 0;
    fin >> look.v[0] >> look.v[1] >> look.v[2];
    look.v[3] = 0;
    fin >> up.v[0] >> up.v[1] >> up.v[2];
    up.v[3] = 0;
    fin >> perspective.v[0] >> perspective.v[1] >> perspective.v[2] >> perspective.v[3];

    Vector4f l = {look.v[0] - eye.v[0], look.v[1] - eye.v[1], look.v[2] - eye.v[2], 0};
    float l_magnitude = sqrt(l.v[0] * l.v[0] + l.v[1] * l.v[1] + l.v[2] * l.v[2]);
    l.v[0] /= l_magnitude;
    l.v[1] /= l_magnitude;
    l.v[2] /= l_magnitude;

    Vector4f r = crossProduct(l, up);
    float r_magnitude = sqrt(r.v[0] * r.v[0] + r.v[1] * r.v[1] + r.v[2] * r.v[2]);
    r.v[0] /= r_magnitude;
    r.v[1] /= r_magnitude;
    r.v[2] /= r_magnitude;

    Vector4f u = crossProduct(r, l);

    Matrix4f T = generateTranslationMatrix(-eye.v[0], -eye.v[1], -eye.v[2]);
    Matrix4f R = {r.v[0], r.v[1], r.v[2], 0,
                  u.v[0], u.v[1], u.v[2], 0,
                  -l.v[0], -l.v[1], -l.v[2], 0,
                  0, 0, 0, 1};
    Matrix4f V = multiplyMatrix(R, T);

    Matrix4f P = generateProjectionMatrix(perspective.v[0], perspective.v[1], perspective.v[2], perspective.v[3]);

    string command;
    while (fin >> command)
    {
        if (command == "triangle")
        {
            Vector4f points[3];
            for (int i = 0; i < 3; ++i)
            {
                fin >> points[i].v[0] >> points[i].v[1] >> points[i].v[2];
                points[i].v[3] = 1;
                points[i] = transformPoint(M, points[i]);
                f_out.precision(7);
                f_out << fixed << points[i].v[0] << " " << points[i].v[1] << " " << points[i].v[2] << endl;
                Vector4f transformedPoint = transformPoint(V, points[i]);
                f_out2.precision(7);
                f_out2 << fixed << transformedPoint.v[0] << " " << transformedPoint.v[1] << " " << transformedPoint.v[2] << endl;
                Vector4f projectedPoint = transformPoint(P, transformedPoint);
                f_out3.precision(7);
                f_out3 << fixed << projectedPoint.v[0] << " " << projectedPoint.v[1] << " " << projectedPoint.v[2] << endl;
            }
            f_out << endl;
            f_out2 << endl;
            f_out3 << endl;
        }
        else if (command == "translate")
        {
            float tx, ty, tz;
            fin >> tx >> ty >> tz;
            M = multiplyMatrix(M, generateTranslationMatrix(tx, ty, tz));
        }
        else if (command == "scale")
        {
            float sx, sy, sz;
            fin >> sx >> sy >> sz;
            M = multiplyMatrix(M, generateScalingMatrix(sx, sy, sz));
        }
        else if (command == "rotate")
        {
            float angle, ax, ay, az;
            fin >> angle >> ax >> ay >> az;
            M = multiplyMatrix(M, generateRotationMatrix(angle, ax, ay, az));
        }
        else if (command == "push")
        {
            S.push(M);
        }
        else if (command == "pop")
        {
            if (!S.empty())
            {
                M = S.top();
                S.pop();
            }
        }
        else if (command == "end")
        {
            break;
        }
    }
    fin.close();
    f_out.close();
    f_out2.close();
    f_out3.close();
    z_buffer_main();
    return 0;
}