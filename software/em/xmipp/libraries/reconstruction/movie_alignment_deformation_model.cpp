/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
 *             David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "reconstruction/movie_alignment_deformation_model.h"

void ProgMovieAlignmentDeformationModel::readParams()
{
	fnMovie = getParam("-i");
	fnMicrograph = getParam("-o");
    initDose = getDoubleParam("--initDose");
    perFrameDose = getDoubleParam("--perFrameDose");
    maxIterations = getIntParam("--maxIterations");
    upScaling = getIntParam("--upscaling");
    fnUnaligned = getParam("--ounaligned");
    threadNumbers = getIntParam("-j");
    fnDark = getParam("--dark");
    fnGain = getParam("--gain");
    show();
}

void ProgMovieAlignmentDeformationModel::show()
{
    if (!verbose)
        return;
    std::cout 
    << "Input movie:          " << fnMovie           << std::endl
    << "Output micrograph:    " << fnMicrograph      << std::endl
    << "Initial dose:         " << initDose          << std::endl
    << "Per frame dose:       " << perFrameDose      << std::endl
    << "Max iterations:       " << maxIterations     << std::endl
    << "Up scaling coef.:     " << upScaling         << std::endl
	<< "Unaligned micrograph: " << fnUnaligned       << std::endl
    << "Threads number:       " << threadNumbers     << std::endl
    << "Dark image:           " << fnDark            << std::endl
    << "Gain image:           " << fnGain            << std::endl;
}

void ProgMovieAlignmentDeformationModel::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <metadata>               : Metadata with the list of frames to align");
    addParamsLine("   -o <fn=\"\"> 		          : Give the name of a micrograph to generate an aligned micrograph");
    addParamsLine("  [--initDose <N=0>]           : Radiation dose received before first frame is taken");
    addParamsLine("  [--perFrameDose <s=0>]       : Radiation dose received after imaging each frame");
    addParamsLine("  [--maxIterations <s=20>]	  : Number of robust least squares iterations");
    addParamsLine("  [--upscaling <N=1>]          : UpScaling coefficient for super resolution image generated from model application");
    addParamsLine("  [--ounaligned <fn=\"\">]     : Give the name of a micrograph to generate an unaligned (initial) micrograph");
    addParamsLine("  [-j <N=5>]                   : Maximum threads the program is allowed to use");
    addParamsLine("  [--dark <fn=\"\">]           : Dark correction image");
    addParamsLine("  [--gain <fn=\"\">]           : Gain correction image");
}

void ProgMovieAlignmentDeformationModel::run()
{  
    loadMovie(fnMovie, frames, timeStamps, fnDark, fnGain);
    
    if (!fnUnaligned.isEmpty()) {
        averageFrames(frames, unalignedMicrograph);
        saveMicrograph(fnUnaligned, unalignedMicrograph);
    }

    globalShiftsX.clear();
    globalShiftsX.resize(frames.size(), 0.0);
    globalShiftsY.clear();
    globalShiftsY.resize(frames.size(), 0.0);
    std::cout << "Estimating global shifts" << std::endl;
    estimateShifts(frames, globalShiftsX, globalShiftsY, maxIterations);
    std::cout << "Applying global shifts" << std::endl;
    applyShifts(frames, globalShiftsX, globalShiftsY);

    std::cout << "Partitioning" << std::endl;
    partitions.clear();
    partitions.resize(PARTITION_AXIS_COUNT * PARTITION_AXIS_COUNT);
    for (int i = 0; i < partitions.size(); i++) {
        partitions[i].resize(frames.size());
    }
    partitionFrames(frames, partitions, PARTITION_AXIS_COUNT);

    localShiftsX.clear();
    localShiftsX.resize(frames.size() * PARTITION_AXIS_COUNT * PARTITION_AXIS_COUNT, 0.0);
    localShiftsY.clear();
    localShiftsY.resize(frames.size() * PARTITION_AXIS_COUNT * PARTITION_AXIS_COUNT, 0.0);
    std::cout << "Estimating local shifts" << std::endl;
    estimateLocalShifts(partitions, localShiftsX, localShiftsY);
    
    deformationCoefficientsX.clear();
    deformationCoefficientsX.resize(9, 0.0);
    deformationCoefficientsY.clear();
    deformationCoefficientsY.resize(9, 0.0);
    std::cout << "Calculating deformation model coefficients" << std::endl;
    calculateModelCoefficients(localShiftsX, timeStamps,
            deformationCoefficientsX, frames[0].ydim, frames[0].xdim);
    calculateModelCoefficients(localShiftsY, timeStamps,
            deformationCoefficientsY, frames[0].ydim, frames[0].xdim);

    correctedFrames.clear();
    correctedFrames.resize(frames.size());
    for (int i = 0; i < correctedFrames.size(); i++) {
        correctedFrames[i].initZeros(frames[0]);
        correctedFrames[i].setXmippOrigin();
    }

    std::cout << "Applying local motion correction" << std::endl;
    motionCorrect(frames, correctedFrames, timeStamps, deformationCoefficientsX,
            deformationCoefficientsY, upScaling);

    std::cout << "Frame averaging" << std::endl;
    averageFrames(correctedFrames, correctedMicrograph);

    std::cout << "Saving end result" << std::endl;
    saveMicrograph(fnMicrograph, correctedMicrograph);
}

void ProgMovieAlignmentDeformationModel::loadMovie(FileName fnMovie,
        std::vector<MultidimArray<double> >& frames,
        std::vector<double>& timeStamps, FileName fnDark, FileName fnGain)
{
    Image<double> dark, gain;
    if (not fnDark.isEmpty()) {
        dark.read(fnDark);
    }
    if (not fnGain.isEmpty()) {
        gain.read(fnGain);
    }

    Image<double> movieStack;
    movieStack.read(fnMovie);
    size_t Xdim, Ydim, Zdim, Ndim;
    movieStack.getDimensions(Xdim, Ydim, Zdim, Ndim);
    std::cout << "Loading stack with dimensions (x, y, z, n): ";
    std::cout << Xdim << "," << Ydim << "," << Zdim <<"," << Ndim << std::endl;
    bool switched = false;
    if (Zdim == 1 and Ndim > 1) {
        Zdim = Ndim;
        switched = true;
    }
    frames.resize(Zdim, MultidimArray<double>(Ydim, Xdim));
    for (size_t z = 0; z < Zdim; z++) {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(frames[z]) {
            if (not switched) {
                DIRECT_A2D_ELEM(frames[z], i, j) = DIRECT_NZYX_ELEM(movieStack(),
                                                                0, z, i, j);
            } else {
                DIRECT_A2D_ELEM(frames[z], i, j) = DIRECT_NZYX_ELEM(movieStack(),
                                                                z, 0, i, j);
            }
        }
        frames[z].setXmippOrigin();

        if (dark().xdim > 0) {
            frames[z] -= dark();
        }
        if (gain().xdim > 0) {
            frames[z] *= gain();
        }

        timeStamps.push_back(initDose + perFrameDose * z);
    }
}

void ProgMovieAlignmentDeformationModel::saveMicrograph(FileName fnMicrograph,
        const MultidimArray<double>& micrograph)
{
	Image<double> img(micrograph);
	img.write(fnMicrograph);
}

void ProgMovieAlignmentDeformationModel::calculateShift2(
        const alglib::real_1d_array &c, const alglib::real_1d_array &dim,
        double &func, void *ptr)
{
	double y = dim[0];
	double x = dim[1];
	double t = dim[2];
	func = (c[0] + c[1]*x + c[2]*x*x + c[3]*y + c[4]*y*y + c[5]*x*y) *
			(c[6]*t + c[7]*t*t + c[8]*t*t*t);
}

double ProgMovieAlignmentDeformationModel::calculateShift(double x, double y, 
        double t, const std::vector<double>& c)
{
    return (c[0] + c[1]*x + c[2]*x*x + c[3]*y + c[4]*y*y + c[5]*x*y) * 
    		(c[6]*t + c[7]*t*t + c[8]*t*t*t);
}

void ProgMovieAlignmentDeformationModel::applyShifts(
        std::vector<MultidimArray<double> >& data,
        const std::vector<double>& shiftsX, const std::vector<double>& shiftsY)
{
	MultidimArray<double> helper;
	helper.initZeros(data[0]);
	helper.setXmippOrigin();
	for (int i = 0; i < data.size(); i++) {
		translate(3, helper, data[i], vectorR2(shiftsX[i], shiftsY[i]),
                false, 0.0);
		data[i] = helper;
	}
}

double ProgMovieAlignmentDeformationModel::linearInterpolation(double y1,
        double x1, double y2, double x2, double v11, double v12, double v21,
        double v22, double p_y, double p_x)
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

void ProgMovieAlignmentDeformationModel::partitionFrames(
        const std::vector<MultidimArray<double> >& frames,
        std::vector<std::vector<MultidimArray<double> > >& partitions,
        int edgeCount)
{
    // properly resize the individual partitions
    int partSizeY = YSIZE(frames[0]) / edgeCount;
    int partSizeX = XSIZE(frames[0]) / edgeCount;
    int yReminder = YSIZE(frames[0]) - (partSizeY * edgeCount);
    int xReminder = XSIZE(frames[0]) - (partSizeX * edgeCount);

    for (int i = 0; i < partitions.size(); i++) {
        int partY = i / edgeCount;
        int partX = i % edgeCount;
        int xSize = partSizeX + (partX < xReminder ? 1 : 0);
        int ySize = partSizeY + (partY < yReminder ? 1 : 0);
        for (int j = 0; j < partitions[i].size(); j++) {
            partitions[i][j].resize(1, 1, ySize, xSize);
            partitions[i][j].initZeros();
            partitions[i][j].setXmippOrigin();
        }
    }

    int longerPartY = yReminder * (partSizeY + 1);
    int longerPartX = xReminder * (partSizeX + 1);
    for (int fi = 0; fi < frames.size(); fi++) {
        FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(frames[fi]) {
            int partY, partX; // index of the partition
            int innerY, innerX; // index inside the partition
            
            if (i < longerPartY) { //in this part, partitions are longer along y
                partY = i / (partSizeY + 1);
                innerY = i % (partSizeY + 1);
            } else {
                partY = yReminder + (i - longerPartY) / partSizeY;
                innerY = (i - longerPartY) % partSizeY;
            }

            if (j < longerPartX) { //in this part, partitions are longer along x
                partX = j / (partSizeX + 1);
                innerX = j % (partSizeX + 1);
            } else {
                partX = xReminder + (j - longerPartX) / partSizeX;
                innerX = (j - longerPartX) % partSizeX;
            }
            dAij(partitions[partY * edgeCount + partX][fi], innerY, innerX) = dAij(frames[fi], i, j);
        }
    }
}

void ProgMovieAlignmentDeformationModel::estimateShifts(
        const std::vector<MultidimArray<double> >& data,
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
    }

    // estimate the shifts
    double shiftX, shiftY;
    CorrelationAux aux;
    for (int cycle = 0; cycle < maxIterations; cycle++) {
        std::cout << "Cycle: " << cycle << std::endl;
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

void ProgMovieAlignmentDeformationModel::estimateLocalShifts(
        const std::vector<std::vector<MultidimArray<double> > >& partitions,
        std::vector<double>& shiftsX, std::vector<double>& shiftsY)
{
    //shiftsX and shiftsY contains shifts for all partitions []
    //shifts are organized 
	int partsPerFrame = partitions.size();
	int partDepth = partitions[0].size(); //frame count
    std::vector<double> tmpXShifts(partDepth);
    std::vector<double> tmpYShifts(partDepth);
    for (int i = 0; i < partitions.size(); i++) {
        std::cout << "Local movement estimation for partition " << i
            << std::endl;
        estimateShifts(partitions[i], tmpXShifts, tmpYShifts, maxIterations);
        for (int j = 0; j < partDepth; j++) {
        	shiftsX[i + j*partsPerFrame] = tmpXShifts[j];
        	shiftsY[i + j*partsPerFrame] = tmpYShifts[j];
        }
    }
}

void ProgMovieAlignmentDeformationModel::calculateModelCoefficients(
        const std::vector<double>& shifts,
        const std::vector<double>& timeStamps, std::vector<double>& coeffs,
        int frameHeight, int frameWidth)
{
	alglib::real_1d_array c;
	c.setlength(coeffs.size());
	for (size_t i = 0; i < c.length(); i++) {
		c[i] = 0.05;
	}

	alglib::real_1d_array y;
	y.setlength(shifts.size());
	for (size_t i = 0; i < y.length(); i++) {
		y[i] = shifts[i];
	}


	alglib::real_2d_array positions;
	positions.setlength(shifts.size(), 3);
	int partsInFrame = PARTITION_AXIS_COUNT * PARTITION_AXIS_COUNT;
    int curFrame = -1;
    double cummulativeX, cummulativeY;
    int partSizeX, partSizeY;
	for (size_t i = 0; i < positions.rows(); i++) {
		int frameIndex = i / partsInFrame;
		int partIndex = i % partsInFrame;
        if (partIndex == 0) { //starting next frame
            cummulativeY = 0;
            cummulativeX = 0;
        } else if (partIndex % PARTITION_AXIS_COUNT == 0) { //new partition line
            cummulativeX = 0;
            cummulativeY += partSizeY;
        }

        calculatePartitionSize(partIndex, PARTITION_AXIS_COUNT, frameHeight,
                frameWidth, partSizeX, partSizeY);

		positions[i][0] = cummulativeY + (partSizeY / 2.0);
		positions[i][1] = cummulativeX + (partSizeX / 2.0);
		positions[i][2] = timeStamps[frameIndex];

        cummulativeX += partSizeX;
	}


    alglib::ae_int_t info;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    double epsf = 0.000001;
    double epsx = 0.000001;
    alglib::ae_int_t maxits = 0;
    double diffstep = 0.000001;


    alglib::lsfitcreatef(positions, y, c, diffstep, state);
    alglib::lsfitsetcond(state, epsf, epsx, maxits);
    alglib::lsfitfit(state, calculateShift2);
    alglib::lsfitresults(state, info, c, rep);

    for (size_t i = 0; i < c.length(); i++) {
    	coeffs[i] = c[i];
    }
}

void ProgMovieAlignmentDeformationModel::motionCorrect(
        const std::vector<MultidimArray<double> >& input,
		std::vector<MultidimArray<double> >& output,
        const std::vector<double>& timeStamps, const std::vector<double>& cx,
        const std::vector<double>& cy, int scalingFactor)
{
    MultidimArray<double> helper;
    helper.resize(input[0].ydim * scalingFactor, input[0].xdim * scalingFactor);
    helper.setXmippOrigin();
	for (int i = 0; i < input.size(); i++) {
        output[i].setXmippOrigin();
        if (scalingFactor != 1) {
    		applyDeformation(input[i], helper, cx, cy, timeStamps[i], 0);
            downsampleFrame(helper, output[i], scalingFactor);
        } else {
            applyDeformation(input[i], output[i], cx, cy, timeStamps[i], 0);
        }
    }
}

void ProgMovieAlignmentDeformationModel::downsampleFrame(
        MultidimArray<double>& in, MultidimArray<double>& out,
        int scalingFactor)
{
    int nThreads = 5;
    double omin=0.,omax=0.;
    in.computeDoubleMinMax(omin,omax);

    FourierTransformer trInput;
    trInput.setThreadsNumber(nThreads);
    MultidimArray<std::complex<double> > inputReciprocal;
    trInput.FourierTransform(in, inputReciprocal, false);

    FourierTransformer trResult;
    trResult.setThreadsNumber(nThreads);
    trResult.setReal(out);
    //TODO: is it calculated??
    MultidimArray<std::complex<double> > resultReciprocal;
    trResult.getFourierAlias(resultReciprocal);

    for (int i = 0; i < resultReciprocal.ydim; i++) {
        for (int j = 0; j < resultReciprocal.xdim; j++) {
            dAij(resultReciprocal, i, j) = dAij(inputReciprocal,
                    i*scalingFactor, j*scalingFactor);
        }
    }
    trResult.inverseFourierTransform();
    out.rangeAdjust(omin, omax);
}

void ProgMovieAlignmentDeformationModel::applyDeformation(
        const MultidimArray<double>& input, MultidimArray<double>& output,
        const std::vector<double>& cx, const std::vector<double>& cy,
		double t1, double t2)
{
    int scalingFactor = upScaling;
    int maxX = input.xdim;
    int maxY = input.ydim;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(input)
    {
        int y = i;
        int x = j;

        int realY = floor(y / scalingFactor);
        int realX = floor(x / scalingFactor);

        double posy = (y + calculateShift(realX, realY, t2, cy) -
                calculateShift(realX, realY, t1, cy));
        double posx = (x + calculateShift(realX, realY, t2, cx) -
                calculateShift(realX, realY, t1, cx));

        int x_left = floor(posx);
        int x_right = x_left + 1;
        int y_down = floor(posy);
        int y_up = y_down + 1;
        //TODO: scaling
        double val;
        if ((x_right < 0 or x_right >= maxX) or (x_left < 0 or x_left >= maxX)
                or (y_up < 0 or y_up >=maxY) or (y_down < 0 or y_down >= maxY)) {
            val = 0;
        } else {
            double v11 = dAij(input, y_down, x_left);
            double v12 = dAij(input, y_down, x_right);
            double v21 = dAij(input, y_up, x_left);
            double v22 = dAij(input, y_up, x_right);
            val = linearInterpolation(y_down, x_left, y_up, x_right, v11, v12,
                    v21, v22, posy, posx);
        }

        DIRECT_A2D_ELEM(output, i, j) = val;
    }	
}

void ProgMovieAlignmentDeformationModel::averageFrames(
        const std::vector<MultidimArray<double> >& data,
        MultidimArray<double>& out)
{
	out.initZeros(data[0]);
    out.setXmippOrigin();
    for (int i = 0; i < data.size(); i++) {
        out += data[i];
    }
}

void ProgMovieAlignmentDeformationModel::calculatePartitionSize(int partIndex,
        int edgeCount, int frameHeight, int frameWidth, int& partXSize,
        int& partYSize)
{
	partYSize = frameHeight / edgeCount;
    partXSize = frameWidth / edgeCount;
    int yReminder = frameHeight - (partYSize * edgeCount);
    int xReminder = frameWidth - (partXSize * edgeCount);

    int partIndexY = partIndex / edgeCount;
    int partIndexX = partIndex % edgeCount;
    partYSize = partYSize + (partIndexY < yReminder ? 1 : 0);
    partXSize = partXSize + (partIndexX < xReminder ? 1 : 0);
}
