#include <bits/stdc++.h>
#include <fstream>
#include <vector>
#include "bitmap_image.hpp"
struct Point
{
    double x, y, z;
};

struct Triangle
{
    Point points[3];
    int color[3];
};

static unsigned long int g_seed = time(0);
inline int random()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

void random_color(int color[3])
{
    color[0] = random() % 256;
    color[1] = random() % 256;
    color[2] = random() % 256;
}
void read_config(int &Screen_Width, int &Screen_Height)
{
    std::ifstream config_file("config.txt");
    config_file >> Screen_Width >> Screen_Height;
    config_file.close();
}

void read_triangles(std::vector<Triangle> &triangles)
{
    std::ifstream stage3_file("stage3.txt");
    while (stage3_file)
    {
        Triangle triangle;
        for (int i = 0; i < 3; ++i)
        {
            stage3_file >> triangle.points[i].x >> triangle.points[i].y >> triangle.points[i].z;
        }
        if (!stage3_file)
        {
            break;
        }
        random_color(triangle.color);
        triangles.push_back(triangle);
    }
    stage3_file.close();
    //std::cout << "Number of triangles: " << triangles.size() << std::endl;
}
void write_and_clear_zbuffer(std::vector<std::vector<double>> &z_buffer, int Screen_Width, int Screen_Height)
{
    std::ofstream z_buffer_file("z_buffer.txt");
    double z_max = 1;
    for (int i = 0; i < Screen_Height; ++i)
    {
        for (int j = 0; j < Screen_Width; ++j)
        {
            if (z_buffer[i][j] < z_max)
            {
                z_buffer_file << z_buffer[i][j] << "\t";
            }
        }
        z_buffer_file << "\n";
    }
    z_buffer_file.close();

    for (int i = 0; i < Screen_Height; ++i)
    {
        z_buffer[i].clear();
    }
    z_buffer.clear();
}
void z_buffer_main()
{
    int Screen_Width, Screen_Height;
    std::vector<Triangle> triangles;
    read_config(Screen_Width, Screen_Height);
    read_triangles(triangles);
    double dx = 2.0 / Screen_Width, dy = 2.0 / Screen_Height;
    double Top_Y = 1 - dy / 2.0, Left_X = -1 + dx / 2.0, Right_X = 1 - dx / 2.0, Bottom_Y = -1 + dy / 2.0;

    std::vector<std::vector<double>> z_buffer(Screen_Height, std::vector<double>(Screen_Width, 1.0));

    bitmap_image image(Screen_Width, Screen_Height);
    image.set_all_channels(0, 0, 0);

    for (auto &triangle : triangles)
    {
        double max_y = std::min(Top_Y, std::max(triangle.points[0].y, std::max(triangle.points[1].y, triangle.points[2].y)));
        double min_y = std::max(Bottom_Y, std::min(triangle.points[2].y, std::min(triangle.points[0].y, triangle.points[1].y)));
        double max_x = std::min(std::max(std::max(triangle.points[0].x, triangle.points[1].x), triangle.points[2].x), Right_X);
        double min_x = std::max(std::min(std::min(triangle.points[0].x, triangle.points[1].x), triangle.points[2].x), Left_X);

        int top_scanline = round((Top_Y - max_y) / dy);
        int bottom_scanline = round((Top_Y - min_y) / dy);

        for (int j = top_scanline; j <= bottom_scanline; j++)
        {
            double y = Top_Y - j * dy;
            std::vector<double> x_list, z_list;

            for (int k = 0; k < 3; k++)
            {
                double x1 = triangle.points[k].x, y1 = triangle.points[k].y;
                double x2 = triangle.points[(k + 1) % 3].x, y2 = triangle.points[(k + 1) % 3].y;
                double z1 = triangle.points[k].z, z2 = triangle.points[(k + 1) % 3].z;
                if (y1 == y2)
                    continue;
                if (y1 > y2)
                {
                    std::swap(x1, x2);
                    std::swap(y1, y2);
                    std::swap(z1, z2);
                }
                if (y < y1 || y > y2)
                    continue;
                double x = x1 + (1.0 * (y - y1) * (x2 - x1)) / (y2 - y1);
                x_list.push_back(x);
                double z = z1 + (1.0 * (y - y1) * (z2 - z1)) / (y2 - y1);
                z_list.push_back(z);
            }
            if (x_list.size() < 2)
                continue;
            if (x_list[0] > x_list[1])
            {
                x_list[0] = x_list[0] + x_list[1];
                x_list[1] = x_list[0] - x_list[1];
                x_list[0] = x_list[0] - x_list[1];
                z_list[0] = z_list[0] + z_list[1];
                z_list[1] = z_list[0] - z_list[1];
                z_list[0] = z_list[0] - z_list[1];
            }
            int left_scanline = round((std::max(x_list[0], min_x) - Left_X) / dx);
            int right_scanline = round((std::min(x_list[1], max_x) - Left_X) / dx);
            double x1 = std::max(x_list[0], min_x), x2 = std::min(x_list[1], max_x), z1 = z_list[0], z2 = z_list[1];

            for (int k = left_scanline; k <= right_scanline; k++)
            {
                double x = Left_X + k * dx;
                double z = z1 + (1.0 * (x - x1) * (z2 - z1)) / (x2 - x1);
                if (z < z_buffer[j][k] && z >= -1.0)
                {
                    z_buffer[j][k] = z;
                    image.set_pixel(k, j, triangle.color[0], triangle.color[1], triangle.color[2]);
                }
            }
        }
    }
    image.save_image("out.bmp");
    image.clear();
    image.setwidth_height(0, 0);
    write_and_clear_zbuffer(z_buffer, Screen_Width, Screen_Height);
}