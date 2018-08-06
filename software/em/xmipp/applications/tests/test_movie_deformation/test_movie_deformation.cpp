#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>
#include <data/metadata_extension.h>
#include <data/xmipp_fftw.h>
#include <data/filters.h>
#include <iostream>
#include <gtest/gtest.h>
// MORE INFO HERE: http://code.google.com/p/googletest/wiki/AdvancedGuide
// This test is named "Size", and belongs to the "MetadataTest"
// test case.
class FuncTest : public ::testing::Test
{
};

void addSquare(MultidimArray<double>& data, int x, int y) {
    const int SIZE = 3;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(data) {
        if ((i < y or i > y + SIZE) or (j < x or j > x + SIZE)) {
            DIRECT_A2D_ELEM(data, i, j) = 0;
        } else {
            DIRECT_A2D_ELEM(data, i, j) = 1;
        }
    }
}

void printArray(MultidimArray<double>& data) {
    for (int y = 0; y < 50; y++) {
        for (int x = 0; x < 50; x++) {
            double val = DIRECT_A2D_ELEM(data, y, x); 
            if (abs(val) < 0.5) {
                val = 0;
            }
            std::cout << val << ",";
        }
        std::cout << std::endl;
    }
}

void applyShifts(std::vector<MultidimArray<double> >& data,
        std::vector<MultidimArray<double> >& out,
        const std::vector<double>& shiftsX, const std::vector<double>& shiftsY)
{
	for (int i = 0; i < data.size(); i++) {
		translate(2, out[i], data[i], vectorR2(shiftsX[i], shiftsY[i]),
                false, 0.0);
	}
}

void averageFrames(const std::vector<MultidimArray<double> >& data,
        MultidimArray<double>& out)
{
	out.initZeros(data[0]);
    out.setXmippOrigin();
    for (int i = 0; i < data.size(); i++) {
        out += data[i];
    }
}

void estimateShifts(const std::vector<MultidimArray<double> >& data,
		std::vector<double>& shiftsX, std::vector<double>& shiftsY,
        int maxIterations)
{
	// prepare sum of images
    MultidimArray<double> sum;
    averageFrames(data, sum);

    std::vector<MultidimArray<double> > shiftedData;
    shiftedData.resize(data.size());
    for (int i = 0; i < shiftedData.size(); i++) {
        shiftedData[i].initZeros(data[i]);
        shiftedData[i] += data[i];
        printArray(shiftedData[i]);
    }

    // estimate the shifts
    double shiftX, shiftY;
    CorrelationAux aux;
    for (int cycle = 0; cycle < maxIterations; cycle++) {
        for (int i = 0; i < data.size(); i++) {
            sum -= shiftedData[i];
            //sum = sum / (data.size() - 1);
            bestShift(sum, data[i], shiftX, shiftY, aux);
            //sum = sum * (data.size() - 1);
            sum += shiftedData[i];

            shiftsY[i] = shiftY;
            shiftsX[i] = shiftX;
        }

        // recalculate avrg
        sum.initZeros();
        for (int i = 0; i < data.size(); i++) {
            translate(3, shiftedData[i], data[i],
                    vectorR2(shiftsX[i], shiftsY[i]), false, 0.0);
            sum += shiftedData[i];
        }
    }
}

TEST_F( FuncTest, CompareTwoFiles)
{
    const int SIZE = 3;
    std::vector<MultidimArray<double> > data;
    data.resize(SIZE, MultidimArray<double>(50, 50));

    addSquare(data[0], 10, 10);
    addSquare(data[1], 15, 13);
    addSquare(data[2], 12, 8);
    for (int i = 0; i < data.size(); i++) {
        data[i].setXmippOrigin();
    }

    std::vector<double> shiftsX(SIZE);
    std::vector<double> shiftsY(SIZE);
    std::cout << "CALCULATE" << std::endl;

    estimateShifts(data, shiftsX, shiftsY, 10);

    std::cout << "SHIFTING" << std::endl;

    std::vector<MultidimArray<double> > res;
    res.resize(SIZE, MultidimArray<double>(50, 50));
    for (int i =0; i < res.size(); i++) {
        res[i].setXmippOrigin();
    }
    applyShifts(data, res, shiftsX, shiftsY);

    std::cout << "FINISHED" << std::endl;
    for (int i = 0; i < shiftsX.size(); i++) {
        std::cout << i << ": (" << shiftsX[i] << ", " << shiftsY[i] << ")"
            << std::endl;
    }

    MultidimArray<double> help(50,50);
    help.setXmippOrigin();
    averageFrames(res, help);
    printArray(help);
}


GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
