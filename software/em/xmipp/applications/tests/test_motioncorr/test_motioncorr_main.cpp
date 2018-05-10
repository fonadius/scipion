#include <gtest/gtest.h>
#include <iostream>
#include <data/xmipp_image.h>
#include <data/multidim_array.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>

class myMotioncorrTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    { 
    // Code here will be called immediately after the constructor (right before each test).
    }

    virtual void TearDown() 
    { 
    // Code here will be called immediately after each test (right before the destructor). 
    }

    // Objects declared here can be used by all tests in the test case for Foo. 
};

void addGrid(Image<int> &img) {
    FOR_ALL_ELEMENTS_IN_ARRAY3D(img.data) {
        int val = A3D_ELEM(img.data, k, i, j);

        if (i % 45 <= 5 or j % 45 <=5) {
            val = 0;
        }
        
        A3D_ELEM(img.data, k, i, j) = val;
    }
}

/**
 * @brief Performs linear interpolation
 * @details  With the following argument naming
    y1: v11 ......... v12
        .................
        ...... p ........
        .................
        .................
    y2: v21 ......... v22
        x1             x2
    (y increases from top to down, x increases from left to right)
    Values q_i are expected to be in strict orthogonal grid, where q11 and q21 have the same x coordinates (the same
    applies for q12 and q22) and q11 and q12 have the same y coordinates (the same applies for q21 and q22).
 * 
 * @param y [description]
 * @param x [description]
 * 
 * @return [description]
 */
double linearInterpolation(double y1, double x1, double y2, double x2, double v11, double v12, double v21, double v22,
                            double p_y, double p_x)
{
    double p1, p2, p;
    // x interpolation
    if (x1 == x2)
    {
        p1 = v11;
        p2 = v21;
    }
    else 
    {
        double diffx = (x2 - x1);
        double rat1 = (x2 - p_x) / diffx;
        double rat2 = (p_x - x1) / diffx;

        p1 = v11 * rat1 + v12 * rat2;
        p2 = v21 * rat1 + v22 * rat2;
    }

    // y interpolation
    if (y1 == y2)
    {
        p = p1;
    }
    else
    {
        double diffy = (y2 - y1);
        double rat1 = (y2 - p_y) / diffy;
        double rat2 = (p_y - y1) / diffy;

        p = p1 * rat1 + p2 * rat2;
    }

    return p;
}


// void loadImages(std::vector<String> filePaths) {
//     Image(1024, 768, 1, filePaths.size());
//     for (String s : filePaths) {
//         FileName fn(s);

//     }
// }

double calculateShift(double x, double y, double t, std::vector<double>& c)
{
    double t2 = t * t;
    double t3 = t2 * t;

    return (c[0] + c[1] * x + c[2] * x * x + c[3] * y + c[4] * y * y + c[5] * x * y) *
            (c[6] * t + c[7] * t2 + c[8] * t3);
}

int getValue(MultidimArray<int>& array, int y, int x) {
    if (array.outside(y, x)) {
        return 0;
    }
    return A2D_ELEM(array, y, x);
}

void applyDeformation(MultidimArray<int>& input, MultidimArray<int>& output, std::vector<double>& cy,
    std::vector<double>& cx, double t1, double t2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(input)
    {
        int y = i;
        int x = j;

        double posy = y + calculateShift(y, x, t2, cy) - calculateShift(y, x, t1, cy);
        double posx = x + calculateShift(y, x, t2, cx) - calculateShift(y, x, t1, cx);
        // std::cout << posx << ":" << x << "   " << posy << ":" << y << std::endl;

        double x_left = floor(posx);
        double x_right = x_left + 1;
        double y_down = floor(posy);
        double y_up = y_down + 1;

        double v11 = getValue(input, y_down, x_left);
        double v12 = getValue(input, y_down, x_right);
        double v21 = getValue(input, y_up, x_left);
        double v22 = getValue(input, y_up, x_right);

        A2D_ELEM(output, i, j) = linearInterpolation(y_down, x_left, y_up, x_right, v11,
                                                            v12, v21, v22, posy, posx);
    }
}

void generateCoefficients(std::vector<double>& out) {
    for (int i = 0; i < 9; i++) {
        double random = ((double) (rand() % 100)) / 1000;
        std::cout << random << std::endl;
        out[i] = random;
    }
}

TEST_F(myMotioncorrTest, superTest)
{
    std::cout << "-----------**************------------" << std::endl;
    FileName testFile((String) "/home/fonadius/Downloads/test.tif");
    srand(time(NULL));
    
    Image<int> img1(1024, 768, 1, 1);
    img1.read(testFile);
    Image<int> img2(1024, 768, 1, 1);

    std::vector<double> cy(9, 0);
    std::vector<double> cx(9, 0);
    generateCoefficients(cy);
    generateCoefficients(cx);

    addGrid(img1);

    applyDeformation(img1(), img2(), cy, cx, 0, 1);

    FileName deformedFile((String) "/home/fonadius/Downloads/deformed.jpg");
    img2.write(deformedFile);
}

// TEST_F(myMotioncorrTest, interpolationTest) 
// {
//     // In form (y1, x1, y2, x2, v11, v12, v21, v22, y, x, expected)
//     // In form (q11, q12, q21, q22, y, x, expected)
//     std::vector<double> knownValues = {
//         // all points are one point
//         5.5, 5.5, 5.5, 5.5, 12.3, 12.3, 12.3, 12.3, 5.5, 5.5, 12.3,
//         // all points lie on the same x coordinates
//         1, 1, 2, 1, 1, 1, 3, 3, 1.5, 1.5, 2,
//         // all points lie on the same y coordinates
//         8, 8, 8, 11, 2, 4, 2, 4, 8, 10, 3 + 1.0/3.0,
//         // point in the middle
//         1, 1, 4, 4, 1, 1, 2, 2, 2.5, 2.5, 1.5,
//         // normal points
//         2.1, 1, 4.8, 3.62, 3.5, 1, 2, 4, 3.2, 2.5, 2.5072094995759118,
//         1, 3.12, 5.6, 12, 1, -2, 9.8, 5.1, 1.12, 10, -1.1291186839012923,
//         8.13, 4.82, 22.65, 4.99, 19, 17, 5, 24, 20.1, 4.87, 11.962202236266407
//     };

//     for (int i = 0; i < knownValues.size(); i+=11)
//     {
//         double* item = &knownValues.front() + i;
//         double result = linearInterpolation(item[0], item[1], item[2], item[3], item[4],
//             item[5], item[6], item[7], item[8], item[9]);
//         ASSERT_NEAR(result, item[10], 1e-12);
//     }
// }

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
