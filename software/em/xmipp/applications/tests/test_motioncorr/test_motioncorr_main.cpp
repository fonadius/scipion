#include <gtest/gtest.h>
#include <iostream>
#include <data/xmipp_image.h>
#include <data/multidim_array.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <data/filters.h>
#include <limits>
#include <algorithm>
#include <reconstruction/movie_alignment_deformation_model.h>

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

// void addGrid(Image<int> &img) {
//     FOR_ALL_ELEMENTS_IN_ARRAY3D(img.data) {
//         int val = A3D_ELEM(img.data, k, i, j);

//         if (i % 45 <= 5 or j % 45 <=5) {
//             val = 0;
//         }
        
//         A3D_ELEM(img.data, k, i, j) = val;
//     }
// }

// void generateCoefficients(std::vector<double>& out) {
//     for (int i = 0; i < 9; i++) {
//         double random = ((double) (rand() % 100)) / 10000;
//         std::cout << random << std::endl;
//         out[i] = random;
//     }
// }

void addSquare(MultidimArray<double>& input, int edgeSize, int y, int x) {
    for (int iy = y; iy < y+edgeSize; iy++) {
        for (int ix = x; ix < x+edgeSize; ix++) {
            A2D_ELEM(input, iy, ix) = 1.0;
        }
    }
}

// TEST_F(myMotioncorrTest, testDataLoading) 
// {
//     FileName starFile("/home/fonadius/Downloads/movie.xmd");
//     MetaData movie;
//     readMovie(movie, starFile);

//     MDRow  row;
//     double time;
//     Image<double> frame;
//     FileName fnFrame("/home/fonadius/Downloads/movie01.mrc");
//     FileName folder = starFile.getDir(); //removeFilename()

//     FOR_ALL_OBJECTS_IN_METADATA(movie) //loop through all lines
//     {   
//         movie.getRow(row, __iter.objId); //read line
//         if (row.getValue(MDL_TIME, time))
//         {
//               std::cout << "The sampling rate is: " <<  time << std::endl;
//         }
//         else
//         {
//             std::cout << "no value " << std::endl;
//         }
//         row.getValue(MDL_IMAGE, fnFrame);
//         std::cout << "name: " << fnFrame << std::endl;
//         frame.read(folder + fnFrame);
//     }

//     FileName res("/home/fonadius/Downloads/last.jpg");
//     frame.write(res);

// }

TEST_F(myMotioncorrTest, testMultipleShifts) {
    Image<double> img1(256, 512, 1, 1);
    Image<double> img2(256, 512, 1, 1);
    Image<double> img3(256, 512, 1, 1);
    Image<double> img4(256, 512, 1, 1);

    addSquare(img1(), 20, 68, 100);
    addSquare(img2(), 20, 200, 45);
    addSquare(img3(), 20, 172, 437);
    addSquare(img4(), 20, 20, 310);

    img1().setXmippOrigin();
    img2().setXmippOrigin();
    img3().setXmippOrigin();
    img4().setXmippOrigin();
    std::vector<MultidimArray<double> > frames = {img1(), img2(), img3(), img4()};

    std::vector<double> shiftsX = {0, 0, 0, 0};
    std::vector<double> shiftsY = {0, 0, 0, 0};
    calculateAndCorrectForShifts(frames, shiftsX, shiftsY);

    MultidimArray<double> helper;
    helper.initZeros(img1());
    MultidimArray<double> sum;
    sum.initZeros(helper);

    applyShift(img1(), helper, shiftsX[0], shiftsY[0]);
    sum += helper;
    applyShift(img2(), helper, shiftsX[1], shiftsY[1]);
    sum += helper;
    applyShift(img3(), helper, shiftsX[2], shiftsY[2]);
    sum += helper;
    applyShift(img4(), helper, shiftsX[3], shiftsY[3]);
    sum += helper;

    Image<double> sumImg(sum);

    FileName fn("/home/fonadius/Downloads/sumAll.jpg");
    sumImg.write(fn);
}

TEST_F(myMotioncorrTest, testShiftAndTranslation)
{
    size_t height = 256;
    size_t width = 512;
    size_t edge = 25;
    Image<double> img1(height, width, 1, 1);
    Image<double> img2(height, width, 1, 1);

    addSquare(img1(), edge, 20, 100);
    addSquare(img2(), edge, 58, 90);

    img1().setXmippOrigin();
    img2().setXmippOrigin();

    CorrelationAux aux;
    double shiftX;
    double shiftY;
    bestShift(img1(), img2(), shiftX, shiftY, aux);
    

    ASSERT_NEAR(shiftY, -38, 1e-10);
    ASSERT_NEAR(shiftX, 10, 1e-10);

    Image<double> img2Shifted(height, width, 1, 1);
    img2Shifted().setXmippOrigin();
    translate(2, img2Shifted(), img2(), vectorR2(shiftX, shiftY), false, 0.0);

    Image<double> result(height, width, 1, 1);
    result() = img1() + img2Shifted();
    size_t nonZeroPixels = result().countThreshold("above", 0.1, 0.0);
    ASSERT_EQ(nonZeroPixels, edge*edge);
    
}



// TEST_F(myMotioncorrTest, superTest)
// {
//     std::cout << "-----------**************------------" << std::endl;
//     FileName testFile((String) "/home/fonadius/Downloads/test.tif");
//     srand(time(NULL));
    
//     Image<int> img1(1024, 768, 1, 1);
//     img1.read(testFile);
//     Image<int> img2(1024, 768, 1, 1);

//     std::vector<double> cy(9, 0);
//     std::vector<double> cx(9, 0);
//     generateCoefficients(cy);
//     generateCoefficients(cx);

//     addGrid(img1);

//     applyDeformation(img1(), img2(), cy, cx, 0, 1);

//     FileName deformedFile((String) "/home/fonadius/Downloads/deformed.jpg");
//     img2.write(deformedFile);
// }

TEST_F(myMotioncorrTest, interpolationTest) 
{
    // In form (y1, x1, y2, x2, v11, v12, v21, v22, y, x, expected)
    // In form (q11, q12, q21, q22, y, x, expected)
    std::vector<double> knownValues = {
        // all points are one point
        5.5, 5.5, 5.5, 5.5, 12.3, 12.3, 12.3, 12.3, 5.5, 5.5, 12.3,
        // all points lie on the same x coordinates
        1, 1, 2, 1, 1, 1, 3, 3, 1.5, 1.5, 2,
        // all points lie on the same y coordinates
        8, 8, 8, 11, 2, 4, 2, 4, 8, 10, 3 + 1.0/3.0,
        // point in the middle
        1, 1, 4, 4, 1, 1, 2, 2, 2.5, 2.5, 1.5,
        // normal points
        2.1, 1, 4.8, 3.62, 3.5, 1, 2, 4, 3.2, 2.5, 2.5072094995759118,
        1, 3.12, 5.6, 12, 1, -2, 9.8, 5.1, 1.12, 10, -1.1291186839012923,
        8.13, 4.82, 22.65, 4.99, 19, 17, 5, 24, 20.1, 4.87, 11.962202236266407
    };

    ProgMovieAlignmentDeformationModel padm;

    for (int i = 0; i < knownValues.size(); i+=11)
    {
        double* item = &knownValues.front() + i;
        double result = padm.linearInterpolation(item[0], item[1], item[2], item[3], item[4],
            item[5], item[6], item[7], item[8], item[9]);
        ASSERT_NEAR(result, item[10], 1e-12);
    }
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
