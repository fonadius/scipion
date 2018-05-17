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

// TEST_F(myMotioncorrTest, testPartitioning) {
//     //prepare frame data
//     int partitionCount = 5;
//     std::vector<MultidimArray<double>> frames;
//     frames.resize(10);
//     bool first = true;
//     int width = 112;
//     int height = 98;

//     int partSizeY = height / partitionCount;
//     int partSizeX = width / partitionCount;
//     int yReminder = height - (partSizeY * partitionCount);
//     int xReminder = width - (partSizeX * partitionCount);
//     int longerPartY = yReminder * (partSizeY + 1);
//     int longerPartX = xReminder * (partSizeX + 1);
//     for (MultidimArray<double>& frame : frames) {
//         frame.resize(height, width);
//         FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(frame) {
//             int partYIndex, partXIndex;
//             if (i < longerPartY) { // still in the part where partitions are longer along y
//                 partYIndex = i / (partSizeY + 1);
//             } else {
//                 partYIndex = yReminder + (i - longerPartY) / partSizeY;
//             }
//             if (j < longerPartX) { //still in the part where partitions are longer along x
//                 partXIndex = j / (partSizeX + 1);
//             } else {
//                 partXIndex = xReminder + (j - longerPartX) / partSizeX;
//             }
//             dAij(frame, i, j) = partYIndex * partitionCount + partXIndex;
//         };
//     }

//     //prepare partition data. After this the vectors should have correct sizes
//     std::vector<std::vector<MultidimArray<double>>> partitions;
//     partitions.resize(partitionCount * partitionCount, {});
//     for (int i = 0; i < partitions.size(); i++) {
//         partitions[i].resize(frames.size());
//     }

//     ProgMovieAlignmentDeformationModel pmadm;
//     pmadm.partitionFrames(frames, partitions, partitionCount);

//     //check correct partition sizes
//     ASSERT_EQ(partitions.size(), partitionCount * partitionCount);
//     for (int i = 0; i < partitions.size(); i++) {
//         ASSERT_EQ(partitions[i].size(), frames.size());
//         int indexY = i / partitionCount;
//         int indexX = i % partitionCount;
//         for (const MultidimArray<double>& part : partitions[i]) {
//             ASSERT_EQ(1, NSIZE(part));
//             ASSERT_EQ(1, ZSIZE(part));
//             if (indexX < 2) {
//                 ASSERT_EQ(23, XSIZE(part));
//             } else {
//                 ASSERT_EQ(22, XSIZE(part));
//             }
//             if (indexY < 3) {
//                 ASSERT_EQ(20, YSIZE(part));
//             } else {
//                 ASSERT_EQ(19, YSIZE(part));
//             }
//         }
//     }

//     //check for correct partition values
//     for (int i = 0; i < partitions.size(); i++) {
//         for (MultidimArray<double>& part : partitions[i]) {
//             int countCorrect = part.countThreshold("range", i - 0.1, i + 0.1);
//             ASSERT_EQ(NZYXSIZE(part), countCorrect);
//         }
//     }
// }

//TODO:
// TEST_F(myMotioncorrTest, testMultipleShifts) {
//     Image<double> img1(256, 512, 1, 1);
//     Image<double> img2(256, 512, 1, 1);
//     Image<double> img3(256, 512, 1, 1);
//     Image<double> img4(256, 512, 1, 1);

//     addSquare(img1(), 20, 68, 100);
//     addSquare(img2(), 20, 200, 45);
//     addSquare(img3(), 20, 172, 437);
//     addSquare(img4(), 20, 20, 310);

//     img1().setXmippOrigin();
//     img2().setXmippOrigin();
//     img3().setXmippOrigin();
//     img4().setXmippOrigin();
//     std::vector<MultidimArray<double> > frames = {img1(), img2(), img3(), img4()};

//     std::vector<double> shiftsX = {0, 0, 0, 0};
//     std::vector<double> shiftsY = {0, 0, 0, 0};
//     calculateAndCorrectForShifts(frames, shiftsX, shiftsY);

//     MultidimArray<double> helper;
//     helper.initZeros(img1());
//     MultidimArray<double> sum;
//     sum.initZeros(helper);

//     applyShift(img1(), helper, shiftsX[0], shiftsY[0]);
//     sum += helper;
//     applyShift(img2(), helper, shiftsX[1], shiftsY[1]);
//     sum += helper;
//     applyShift(img3(), helper, shiftsX[2], shiftsY[2]);
//     sum += helper;
//     applyShift(img4(), helper, shiftsX[3], shiftsY[3]);
//     sum += helper;

//     Image<double> sumImg(sum);

//     FileName fn("/home/fonadius/Downloads/sumAll.jpg");
//     sumImg.write(fn);
// }


TEST_F(myMotioncorrTest, testShiftAndTranslation)
{
    size_t height = 256;
    size_t width = 512;
    size_t edge = 25;
    std::vector<double> xPos = {100, 90, 300};
    std::vector<double> yPos = {20, 58, 120};
    Image<double> img1(height, width, 1, 1);
    Image<double> img2(height, width, 1, 1);
    Image<double> img3(height, width, 1, 1);
    addSquare(img1(), edge, yPos[0], xPos[0]);
    addSquare(img2(), edge, yPos[1], xPos[1]);
    addSquare(img3(), edge, yPos[2], xPos[2]);
    img1().setXmippOrigin();
    img2().setXmippOrigin();
    img3().setXmippOrigin();
    std::vector<MultidimArray<double>> data = {img1(), img2(), img3()};
    std::vector<double> shiftsX(data.size(), 0);
    std::vector<double> shiftsY(data.size(), 0);

    ProgMovieAlignmentDeformationModel pmadm;
    pmadm.estimateShifts(data, shiftsX, shiftsY, 50, 0.1);

    // for (int i = 0; i < data.size(); i++) {
    //     std::cout << i << ": " << "x is " << (xPos[i] + shiftsX[i]) <<
    //     ", y is " << (yPos[i] + shiftsY[i]) << std::endl;
    // }

    for (int i = 1; i < data.size(); i++) {
        ASSERT_NEAR(xPos[0] + shiftsX[0], xPos[i] + shiftsX[i], 1e-10);
        ASSERT_NEAR(yPos[0] + shiftsY[0], yPos[i] + shiftsY[i], 1e-10);
    }
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

//     ProgMovieAlignmentDeformationModel padm;

//     for (int i = 0; i < knownValues.size(); i+=11)
//     {
//         double* item = &knownValues.front() + i;
//         double result = padm.linearInterpolation(item[0], item[1], item[2], item[3], item[4],
//             item[5], item[6], item[7], item[8], item[9]);
//         ASSERT_NEAR(result, item[10], 1e-12);
//     }
// }

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
