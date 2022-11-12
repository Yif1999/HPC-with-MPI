#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <mpi.h>
using namespace std;

// 定义PPM图像类用于读取保存和像素操作
class Image
{
public:
    int width, height, max_value, n_bytes;
    char *data;

    void init(int _width, int _height)
    {
        width = _width;
        height = _height;
        n_bytes = _width * _height * 3;
        data = new char[n_bytes];
        for (int i = 0; i < n_bytes; i += 3)
        {
            data[i] = 0;
            data[i + 1] = 0;
            data[i + 2] = 0;
        }
    }

    void load(const char *filename)
    {
        ifstream file(filename, ios::binary);
        if (!file.is_open())
        {
            cout << "Error opening file" << endl;
            exit(1);
        }

        char magic_number[2];
        file.read(magic_number, 2);
        if (magic_number[0] != 'P' || magic_number[1] != '6')
        {
            cout << "Error: not a binary PPM file" << endl;
            exit(1);
        }

        file >> width >> height >> max_value;
        file.get();
        n_bytes = width * height * 3;
        data = new char[n_bytes];
        file.read(data, n_bytes);

        file.close();
        cout << "Image loaded from " << filename << endl;
    }

    void save(const char *filename)
    {
        ofstream file(filename, ios::binary);
        if (!file.is_open())
        {
            cout << "Error opening file" << endl;
            exit(1);
        }

        file << "P6" << endl;
        file << width << " " << height << endl;
        file << 255 << endl;
        file.write(data, n_bytes);

        file.close();

        cout << "Image saved to " << filename << endl;
    }

    void get_pixel(int x, int y, int &red, int &green, int &blue)
    {
        int index = get_address(x, y);
        red = (int)(unsigned char)(data[index]);
        green = (int)(unsigned char)data[index + 1];
        blue = (int)(unsigned char)data[index + 2];
    }

    void get_luminance(int x, int y, float &luminance)
    {
        int index = get_address(x, y);
        float red = (int)(unsigned char)(data[index]) / 255.0;
        float green = (int)(unsigned char)data[index + 1] / 255.0;
        float blue = (int)(unsigned char)data[index + 2] / 255.0;
        luminance = 0.2125 * red + 0.7154 * green + 0.0721 * blue;
    }

    void set_pixel(int x, int y, int red, int green, int blue)
    {
        int index = get_address(x, y);
        data[index] = red;
        data[index + 1] = green;
        data[index + 2] = blue;
    }

private:
    int get_address(int x, int y)
    {
        if (x < 0 || x >= width || y < 0 || y >= height)
        {
            cout << "Error: (" << x << "," << y << ") pixel out of bounds" << endl;
            exit(1);
        }

        return (y * width + x) * 3;
    }
};

int main(int argc, char **argv)
{
    clock_t start, end;
    start = clock();

    // 初始化 MPI 环境
    int id, p;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // 根节点加载图片并广播图片信息
    int imgSize[3];
    Image file;
    if (!id)
    {
        file.load("input.ppm");
        imgSize[0] = file.width;
        imgSize[1] = file.height;
        imgSize[2] = file.n_bytes;
    }
    MPI_Bcast(imgSize, 3, MPI_INT, 0, MPI_COMM_WORLD);

    // 行分块划分数据
    Image texture;
    texture.init(imgSize[0], imgSize[1] / p);
    MPI_Scatter(file.data, imgSize[2] / p, MPI_BYTE, texture.data, imgSize[2] / p, MPI_BYTE, 0, MPI_COMM_WORLD);

    // 从全分辨率开始
    int tw = texture.width;
    int th = texture.height;

    // 确定迭代次数
    int maxSize = min(tw, th);
    int iterations = max(1, int(log2(maxSize) - 1));
    int mipCount = min(iterations, 7);

    // 创建临时RenderTexture
    Image temporaryRT[2 * mipCount + 1];

    // 亮度钳制预处理
    temporaryRT[0].init(tw, th + 2); // 多开两行空间作为虚拟网格
    for (int x = 0; x < tw; x++)
    {
        for (int y = 0; y < th; y++)
        {
            float l;
            texture.get_luminance(x, y, l);
            if (l > 0.8)
            {
                int r, g, b;
                texture.get_pixel(x, y, r, g, b);
                temporaryRT[0].set_pixel(x, y + 1, r, g, b);
            }
        }
    }

    // 降采样
    for (int i = 1; i <= mipCount; i++)
    {
        tw = tw >> 1;
        th = th >> 1;
        temporaryRT[i].init(tw, th + 2); // 多开两行空间作为虚拟网格

        // 应用Downsample filter
        int r, g, b;
        float R, G, B;
        for (int x = 0; x < tw; x++)
        {
            for (int y = 0; y < th; y++)
            {
                temporaryRT[i - 1].get_pixel(2 * x, 2 * y + 1, r, g, b);
                R = 0.25 * r;
                G = 0.25 * g;
                B = 0.25 * b;
                temporaryRT[i - 1].get_pixel(2 * x + 1, 2 * y + 1, r, g, b);
                R += 0.25 * r;
                G += 0.25 * g;
                B += 0.25 * b;
                temporaryRT[i - 1].get_pixel(2 * x, 2 * y + 2, r, g, b);
                R += 0.25 * r;
                G += 0.25 * g;
                B += 0.25 * b;
                temporaryRT[i - 1].get_pixel(2 * x + 1, 2 * y + 2, r, g, b);
                R += 0.25 * r;
                G += 0.25 * g;
                B += 0.25 * b;
                temporaryRT[i].set_pixel(x, y + 1, R, G, B);
            }
        }

        // 每步降采样结果通信
        if (p > 1)
        {
            if (!id)
            {
                for (int x = 0; x < tw; x++)
                {
                    temporaryRT[i].get_pixel(x, 1, r, g, b);
                    temporaryRT[i].set_pixel(x, 0, r, g, b);
                }
                MPI_Sendrecv(&(temporaryRT[i].data[tw * th * 3]), tw * 3, MPI_BYTE, 1, 0,
                             &(temporaryRT[i].data[tw * (th + 1) * 3]), tw * 3, MPI_BYTE, 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (id == p - 1)
            {
                for (int x = 0; x < tw; x++)
                {
                    temporaryRT[i].get_pixel(x, th, r, g, b);
                    temporaryRT[i].set_pixel(x, th + 1, r, g, b);
                }
                MPI_Sendrecv(&(temporaryRT[i].data[tw * 3]), tw * 3, MPI_BYTE, p - 2, 0,
                             temporaryRT[i].data, tw * 3, MPI_BYTE, p - 2, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else
        {
            for (int x = 0; x < tw; x++)
            {
                temporaryRT[i].get_pixel(x, 1, r, g, b);
                temporaryRT[i].set_pixel(x, 0, r, g, b);
                temporaryRT[i].get_pixel(x, th, r, g, b);
                temporaryRT[i].set_pixel(x, th + 1, r, g, b);
            }
        }

        if (id > 0 && id < p - 1)
        {
            MPI_Sendrecv(&(temporaryRT[i].data[tw * 3]), tw * 3, MPI_BYTE, id - 1, 0,
                         temporaryRT[i].data, tw * 3, MPI_BYTE, id - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&(temporaryRT[i].data[tw * th * 3]), tw * 3, MPI_BYTE, id + 1, 0,
                         &(temporaryRT[i].data[tw * (th + 1) * 3]), tw * 3, MPI_BYTE, id + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // 升采样
    for (int i = mipCount + 1; i <= mipCount * 2; i++)
    {
        tw = tw << 1;
        th = th << 1;
        temporaryRT[i].init(tw, th + 4); // 多开四行作为虚拟网格

        // 建立双线性插值查找表
        for (int x = 0; x < tw; x++)
        {
            for (int y = 0; y < th; y++)
            {
                int interp_r[4], interp_g[4], interp_b[4];
                temporaryRT[2 * mipCount - i + 1].get_pixel(max((int)floor(x / 2.0 - 0.5), 0), floor(y / 2.0 - 0.5) + 1, interp_r[0], interp_g[0], interp_b[0]);
                temporaryRT[2 * mipCount - i + 1].get_pixel(min((int)floor(x / 2.0 - 0.5) + 1, tw / 2 - 1), floor(y / 2.0 - 0.5) + 1, interp_r[1], interp_g[1], interp_b[1]);
                temporaryRT[2 * mipCount - i + 1].get_pixel(max((int)floor(x / 2.0 - 0.5), 0), floor(y / 2.0 - 0.5) + 2, interp_r[2], interp_g[2], interp_b[2]);
                temporaryRT[2 * mipCount - i + 1].get_pixel(min((int)floor(x / 2.0 - 0.5) + 1, tw / 2 - 1), floor(y / 2.0 - 0.5) + 2, interp_r[3], interp_g[3], interp_b[3]);

                int r1, g1, b1, r2, g2, b2;
                float interp_x = 0.25 + 0.5 * (0.5 / tw + 1.0 / tw * x > 2.0 / tw * floor(x / 2.0 + 0.5));
                r1 = interp_r[0] + (interp_r[1] - interp_r[0]) * interp_x;
                r2 = interp_r[2] + (interp_r[3] - interp_r[2]) * interp_x;
                g1 = interp_g[0] + (interp_g[1] - interp_g[0]) * interp_x;
                g2 = interp_g[2] + (interp_g[3] - interp_g[2]) * interp_x;
                b1 = interp_b[0] + (interp_b[1] - interp_b[0]) * interp_x;
                b2 = interp_b[2] + (interp_b[3] - interp_b[2]) * interp_x;

                int r, g, b;
                float interp_y = 0.25 + 0.5 * (0.5 / th + 1.0 / th * y > 2.0 / th * floor(y / 2.0 + 0.5));
                r = r1 + (r2 - r1) * interp_y;
                g = g1 + (g2 - g1) * interp_y;
                b = b1 + (b2 - b1) * interp_y;

                temporaryRT[i].set_pixel(x, y + 2, r, g, b);
            }
        }

        // 查找表边界通信
        if (p > 1)
        {
            if (!id)
            {
                for (int x = 0; x < tw; x++)
                {
                    int r, g, b;
                    temporaryRT[i].get_pixel(x, 2, r, g, b);
                    temporaryRT[i].set_pixel(x, 1, r, g, b);
                    temporaryRT[i].get_pixel(x, 3, r, g, b);
                    temporaryRT[i].set_pixel(x, 0, r, g, b);
                }
                MPI_Sendrecv(&(temporaryRT[i].data[tw * th * 3]), tw * 6, MPI_BYTE, 1, 0,
                             &(temporaryRT[i].data[tw * (th + 2) * 3]), tw * 6, MPI_BYTE, 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (id == p - 1)
            {
                for (int x = 0; x < tw; x++)
                {
                    int r, g, b;
                    temporaryRT[i].get_pixel(x, th + 1, r, g, b);
                    temporaryRT[i].set_pixel(x, th + 2, r, g, b);
                    temporaryRT[i].get_pixel(x, th, r, g, b);
                    temporaryRT[i].set_pixel(x, th + 3, r, g, b);
                }
                MPI_Sendrecv(&(temporaryRT[i].data[tw * 6]), tw * 6, MPI_BYTE, p - 2, 0,
                             temporaryRT[i].data, tw * 6, MPI_BYTE, p - 2, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else
        {
            for (int x = 0; x < tw; x++)
            {
                int r, g, b;
                temporaryRT[i].get_pixel(x, 2, r, g, b);
                temporaryRT[i].set_pixel(x, 1, r, g, b);
                temporaryRT[i].get_pixel(x, 3, r, g, b);
                temporaryRT[i].set_pixel(x, 0, r, g, b);
                temporaryRT[i].get_pixel(x, th + 1, r, g, b);
                temporaryRT[i].set_pixel(x, th + 2, r, g, b);
                temporaryRT[i].get_pixel(x, th, r, g, b);
                temporaryRT[i].set_pixel(x, th + 3, r, g, b);
            }
        }

        if (id > 0 && id < p - 1)
        {
            MPI_Sendrecv(&(temporaryRT[i].data[tw * 6]), tw * 6, MPI_BYTE, id - 1, 0,
                         temporaryRT[i].data, tw * 6, MPI_BYTE, id - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&(temporaryRT[i].data[tw * th * 3]), tw * 6, MPI_BYTE, id + 1, 0,
                         &(temporaryRT[i].data[tw * (th + 2) * 3]), tw * 6, MPI_BYTE, id + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // 应用Upsample filter算子
        int r, g, b;
        float R, G, B;
        for (int x = 0; x < tw; x++)
        {
            for (int y = 2; y < th + 2; y++)
            {
                temporaryRT[i].get_pixel(abs(x - 2), y, r, g, b);
                R = 1.0 / 12 * r;
                G = 1.0 / 12 * g;
                B = 1.0 / 12 * b;
                temporaryRT[i].get_pixel(x, y - 2, r, g, b);
                R += 1.0 / 12 * r;
                G += 1.0 / 12 * g;
                B += 1.0 / 12 * b;
                temporaryRT[i].get_pixel(x + 2 < tw ? x + 2 : 2 * tw - x - 3, y, r, g, b);
                R += 1.0 / 12 * r;
                G += 1.0 / 12 * g;
                B += 1.0 / 12 * b;
                temporaryRT[i].get_pixel(x, y + 2, r, g, b);
                R += 1.0 / 12 * r;
                G += 1.0 / 12 * g;
                B += 1.0 / 12 * b;
                temporaryRT[i].get_pixel(abs(x - 1), y - 1, r, g, b);
                R += 1.0 / 6 * r;
                G += 1.0 / 6 * g;
                B += 1.0 / 6 * b;
                temporaryRT[i].get_pixel(x + 1 < tw ? x + 1 : 2 * tw - x - 2, y - 1, r, g, b);
                R += 1.0 / 6 * r;
                G += 1.0 / 6 * g;
                B += 1.0 / 6 * b;
                temporaryRT[i].get_pixel(abs(x - 1), y + 1, r, g, b);
                R += 1.0 / 6 * r;
                G += 1.0 / 6 * g;
                B += 1.0 / 6 * b;
                temporaryRT[i].get_pixel(x + 1 < tw ? x + 1 : 2 * tw - x - 2, y + 1, r, g, b);
                R += 1.0 / 6 * r;
                G += 1.0 / 6 * g;
                B += 1.0 / 6 * b;

                temporaryRT[mipCount * 2 - i].set_pixel(x, y - 1, R, G, B);
            }
        }

        // 每步升采样结果通信
        if (p > 1)
        {
            if (!id)
            {
                for (int x = 0; x < tw; x++)
                {
                    temporaryRT[mipCount * 2 - i].get_pixel(x, 1, r, g, b);
                    temporaryRT[mipCount * 2 - i].set_pixel(x, 0, r, g, b);
                }
                MPI_Sendrecv(&(temporaryRT[mipCount * 2 - i].data[tw * th * 3]), tw * 3, MPI_BYTE, 1, 0,
                             &(temporaryRT[mipCount * 2 - i].data[tw * (th + 1) * 3]), tw * 3, MPI_BYTE, 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (id == p - 1)
            {
                for (int x = 0; x < tw; x++)
                {
                    temporaryRT[mipCount * 2 - i].get_pixel(x, th, r, g, b);
                    temporaryRT[mipCount * 2 - i].set_pixel(x, th + 1, r, g, b);
                }
                MPI_Sendrecv(&(temporaryRT[mipCount * 2 - i].data[tw * 3]), tw * 3, MPI_BYTE, p - 2, 0,
                             temporaryRT[mipCount * 2 - i].data, tw * 3, MPI_BYTE, p - 2, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else
        {
            for (int x = 0; x < tw; x++)
            {
                temporaryRT[mipCount * 2 - i].get_pixel(x, 1, r, g, b);
                temporaryRT[mipCount * 2 - i].set_pixel(x, 0, r, g, b);
                temporaryRT[mipCount * 2 - i].get_pixel(x, th, r, g, b);
                temporaryRT[mipCount * 2 - i].set_pixel(x, th + 1, r, g, b);
            }
        }

        if (id > 0 && id < p - 1)
        {
            MPI_Sendrecv(&(temporaryRT[mipCount * 2 - i].data[tw * 3]), tw * 3, MPI_BYTE, id - 1, 0,
                         temporaryRT[mipCount * 2 - i].data, tw * 3, MPI_BYTE, id - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&(temporaryRT[mipCount * 2 - i].data[tw * th * 3]), tw * 3, MPI_BYTE, id + 1, 0,
                         &(temporaryRT[mipCount * 2 - i].data[tw * (th + 1) * 3]), tw * 3, MPI_BYTE, id + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // 合并结果
    for (int x = 0; x < tw; x++)
    {
        for (int y = 0; y < th; y++)
        {
            int r, g, b, R, G, B;
            temporaryRT[0].get_pixel(x, y, r, g, b);
            texture.get_pixel(x, y, R, G, B);
            texture.set_pixel(x, y, min(r + R, 255), min(g + G, 255), min(b + B, 255));
        }
    }

    // 收集子进程分布的后处理结果到根节点
    MPI_Gather(texture.data, imgSize[2] / p, MPI_BYTE, file.data, imgSize[2] / p, MPI_BYTE, 0, MPI_COMM_WORLD);

    end = clock();
    cout << "The subprocess cost " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

    if (!id)
        file.save("output_parallel.ppm");

    // 释放 MPI 资源
    MPI_Finalize();

    return 0;
}