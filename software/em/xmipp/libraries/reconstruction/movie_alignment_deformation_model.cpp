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
    maxIterations = getIntParam("--maxIterations");
    upScaling = getIntParam("--upscaling");
    fnUnaligned = getParam("--ounaligned");
    show();
}

void ProgMovieAlignmentDeformationModel::show()
{
    if (!verbose)
        return;
    std::cout 
    << "Input movie:          " << fnMovie           << std::endl
    << "Output micrograph:    " << fnMicrograph      << std::endl
    << "Max iterations:       " << maxIterations     << std::endl
    << "Up scaling coef.:     " << upScaling         << std::endl
	<< "Unaligned micrograph: " << fnUnaligned       << std::endl;
}

void ProgMovieAlignmentDeformationModel::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <metadata>               : Metadata with the list of frames to align");
    addParamsLine("   -o <fn=\"\"> 		          : Give the name of a micrograph to generate an aligned micrograph");
    addParamsLine("  [--maxIterations <N=5>]	  : Number of robust least squares iterations");
    addParamsLine("  [--upscaling <N=1>]         : UpScaling coefficient for super resolution image generated from model application");
    addParamsLine("  [--ounaligned <fn=\"\">]     : Give the name of a micrograph to generate an unaligned (initial) micrograph");
}

void ProgMovieAlignmentDeformationModel::run()
{
    std::cout << "UPPPP" << upScaling << std::endl;
	loadMovie(this->fnMovie, this->frames, this->timeStamps);

    if (!fnUnaligned.isEmpty()) {
        averageFrames(this->frames, this->unalignedMicrograph);
        saveMicrograph(this->fnUnaligned, this->unalignedMicrograph);
    }

    this->globalShiftsX.clear();
    this->globalShiftsX.resize(this->frames.size(), 0.0);
    this->globalShiftsY.clear();
    this->globalShiftsY.resize(this->frames.size(), 0.0);
	estimateShifts(this->frames, this->globalShiftsX, this->globalShiftsY, maxIterations);

	applyShifts(this->frames, this->globalShiftsX, this->globalShiftsY);

    this->partitions.clear();
    this->partitions.resize(this->PARTITION_AXIS_COUNT * this->PARTITION_AXIS_COUNT, {});
    for (int i = 0; i < this->partitions.size(); i++) {
        this->partitions[i].resize(this->frames.size());
    }
    partitionFrames(this->frames, this->partitions, this->PARTITION_AXIS_COUNT);

    this->localShiftsX.clear();
    this->localShiftsX.resize(this->frames.size() * this->PARTITION_AXIS_COUNT * this->PARTITION_AXIS_COUNT, 0.0);
    this->localShiftsY.clear();
    this->localShiftsY.resize(this->frames.size() * this->PARTITION_AXIS_COUNT * this->PARTITION_AXIS_COUNT, 0.0);
	estimateLocalShifts(this->partitions, this->localShiftsX, this->localShiftsY);

    this->deformationCoefficientsX.clear();
    this->deformationCoefficientsX.resize(9, 0.0);
    this->deformationCoefficientsY.clear();
    this->deformationCoefficientsY.resize(9, 0.0);
	calculateModelCoefficients(this->localShiftsX, this->timeStamps,
			this->deformationCoefficientsX, this->frames[0].ydim, this->frames[0].xdim);
	calculateModelCoefficients(this->localShiftsY, this->timeStamps,
			this->deformationCoefficientsY, this->frames[0].ydim, this->frames[0].xdim);

    this->correctedFrames.clear();
    this->correctedFrames.resize(this->frames.size());
    for (MultidimArray<double>& ma : this->correctedFrames) {
        ma.resize(1, 1, this->frames[0].ydim * this->upScaling, this->frames[0].xdim * this->upScaling);
        ma.initZeros();
        ma.setXmippOrigin();
    }
	motionCorrect(this->frames, this->correctedFrames, this->timeStamps, this->deformationCoefficientsX,
			this->deformationCoefficientsY);

    for (int i = 0; i < correctedFrames.size(); i++) {
        Image<double> img;
        img() = correctedFrames[i];
        img.write("/home/fonadius/Downloads/" + std::to_string(i) + ".jpg");
        // for (int j = 0; j < 25; j++) {
        //     std::cout << this->localShiftsX[j] << ",";
        // }
        // std::cout << std::endl;
    }

	averageFrames(this->correctedFrames, this->correctedMicrograph);

	saveMicrograph(this->fnMicrograph, this->correctedMicrograph);
}

void ProgMovieAlignmentDeformationModel::loadMovie(FileName fnMovie, std::vector<MultidimArray<double>>& frames,
		std::vector<double>& timeStamps)
{
    MetaData movie;
    movie.read(fnMovie, NULL, "movie_stack");

    MDRow  row;
    double time;
    FileName fnFrame;
    FileName fnFolder = fnMovie.getDir();
    FOR_ALL_OBJECTS_IN_METADATA(movie)
    {
        movie.getRow(row, __iter.objId);
        // time stamp
        if (row.getValue(MDL_TIME, time)){
            timeStamps.push_back(time);
        } else {
            std::cerr << "Missing value of _time" << std::endl;
        }
        // frame image data
        if (row.getValue(MDL_IMAGE, fnFrame)){
            Image<double> frame;
            frame.read(fnFolder + fnFrame);
            frame().setXmippOrigin();
            frames.push_back(frame());
        } else {
            std::cerr << "Missing value of _image" << std::endl;
        }
    }
}

void ProgMovieAlignmentDeformationModel::saveMicrograph(FileName fnMicrograph, const MultidimArray<double>& micrograph)
{
	Image<double> img(micrograph);
	img.write(fnMicrograph);
}

void ProgMovieAlignmentDeformationModel::calculateShift2(const alglib::real_1d_array &c,
		const alglib::real_1d_array &dim, double &func, void *ptr)
{
	double y = dim[0];
	double x = dim[1];
	double t = dim[2];
	func = (c[0] + c[1]*x + c[2]*x*x + c[3]*y + c[4]*y*y + c[5]*x*y) *
			(c[6]*t + c[7]*t*t + c[8]*t*t*t);
}

double ProgMovieAlignmentDeformationModel::calculateShift(double x, double y, double t, const std::vector<double>& c)
{
    return (c[0] + c[1]*x + c[2]*x*x + c[3]*y + c[4]*y*y + c[5]*x*y) * 
    		(c[6]*t + c[7]*t*t + c[8]*t*t*t);
}

void ProgMovieAlignmentDeformationModel::applyShifts(std::vector<MultidimArray<double>>& data,
	const std::vector<double>& shiftsX, const std::vector<double>& shiftsY)
{
	MultidimArray<double> helper;
	helper.initZeros(data[0]);
	helper.setXmippOrigin();
	for (int i = 0; i < data.size(); i++) {
		translate(2, helper, data[i], vectorR2(shiftsX[i], shiftsY[i]), false, 0.0);
		data[i] = helper;
	}
}

double ProgMovieAlignmentDeformationModel::linearInterpolation(double y1, double x1, double y2, double x2, double v11,
		double v12, double v21, double v22, double p_y, double p_x)
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

void ProgMovieAlignmentDeformationModel::partitionFrames(const std::vector<MultidimArray<double>>& frames,
        std::vector<std::vector<MultidimArray<double>>>& partitions, int edgeCount)
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
        for (MultidimArray<double>& part : partitions[i]) {
            part.resize(1, 1, ySize, xSize);
            part.initZeros();
            part.setXmippOrigin();
        }
    }


    int longerPartY = yReminder * (partSizeY + 1);
    int longerPartX = xReminder * (partSizeX + 1);
    for (int fi = 0; fi < frames.size(); fi++) {
        FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(frames[fi]) {
            int partY, partX; // index of the partition
            int innerY, innerX; // index inside the partition
            
            if (i < longerPartY) { // still in the part where partitions are longer along y
                partY = i / (partSizeY + 1);
                innerY = i % (partSizeY + 1);
            } else {
                partY = yReminder + (i - longerPartY) / partSizeY;
                innerY = (i - longerPartY) % partSizeY;
            }

            if (j < longerPartX) { //still in the part where partitions are longer along x
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

void ProgMovieAlignmentDeformationModel::estimateShifts(const std::vector<MultidimArray<double>>& data,
		std::vector<double>& shiftsX, std::vector<double>& shiftsY, int maxIterations, double minImprovement)
{
	// std::vector<double> yPos = {20, 58, 120};
 //    std::vector<double> xPos = {100, 90, 300};
	
	// prepare sum of images
    MultidimArray<double> sum;
    sum.initZeros(data[0]);
    sum.setXmippOrigin();
    // std::cout << "size: " << data.size() << std::endl;
    // std::cout << "x:" << sum.xdim << std::endl;
    // std::cout << "y:" << sum.ydim << std::endl;
    for (const MultidimArray<double>& ma : data) {
        // std::cout << "x:";
        // std::cout << ma.xdim;
        // std::cout << "y:";
        // std::cout << ma.ydim << std::endl;
        sum = sum + ma;
    }    
    // std::cout << "----end---" << std::endl;
    // estimate the shifts
    double shiftX;
    double shiftY;
    CorrelationAux aux;
    int cycle = 0;
    double maxDiff = std::numeric_limits<double>::max();
    MultidimArray<double> helper;
    // Image<double> img(256,512,1,1);
    while (cycle < maxIterations) {
        maxDiff = 0;
        cycle++;
     //    for (int i = 0; i < data.size(); i++) {
	    //     std::cout << i << ": " << "x is " << (xPos[i] + shiftsX[i]) <<
	    //     	", y is " << (yPos[i] + shiftsY[i]) << std::endl;
    	// }
        for (int i = 0; i < data.size(); i++) {
            sum -= data[i];
            sum = sum / (data.size() - 1);
            // FileName fn("/home/fonadius/Downloads/a" + std::to_string(cycle) + "-" +  std::to_string(i) + ".jpg");
            // img() = sum;
            // img.write(fn);
            bestShift(sum, data[i], shiftX, shiftY, aux);
            sum = sum * (data.size() - 1);
            sum += data[i];

            double diff = std::max(std::abs(shiftsX[i] - shiftX), 
                            std::abs(shiftsY[i] - shiftY));
            maxDiff = std::max(maxDiff, diff);

            shiftsY[i] = shiftY;
            shiftsX[i] = shiftX;

        }

        // recalculate avrg
        sum.initZeros();
        for (int i = 0; i < data.size(); i++) {
            translate(2, helper, data[i], vectorR2(shiftsX[i], shiftsY[i]), false, 0.0);
            sum += helper;
        }

        // std::cout << cycle << " : " << maxDiff << std::endl;
    }
}

void ProgMovieAlignmentDeformationModel::estimateLocalShifts(
        const std::vector<std::vector<MultidimArray<double>>>& partitions,
        std::vector<double>& shiftsX, std::vector<double>& shiftsY)
{
    //shiftsX and shiftsY contains shifts for all partitions []
    //shifts are organized 
	int partsPerFrame = partitions.size();
	int partDepth = partitions[0].size(); //frame count
    std::vector<double> tmpXShifts(partDepth);
    std::vector<double> tmpYShifts(partDepth);
    for (int i = 0; i < partitions.size(); i++) {
    	// std::cout << partitions[i].size() << std::endl;
        estimateShifts(partitions[i], tmpXShifts, tmpYShifts, this->maxIterations);
        for (int j = 0; j < partDepth; j++) {
        	shiftsX[i + j*partsPerFrame] = tmpXShifts[j];
        	shiftsY[i + j*partsPerFrame] = tmpYShifts[j];
        }
    }
}

void ProgMovieAlignmentDeformationModel::calculateModelCoefficients(const std::vector<double>& shifts,
        const std::vector<double>& timeStamps, std::vector<double>& coeffs, int frameHeight, int frameWidth)
{
	alglib::real_1d_array c;
	c.setlength(coeffs.size());
	for (size_t i = 0; i < c.length(); i++) {
		c[i] = 1;
	}

	alglib::real_1d_array y;
	y.setlength(shifts.size());
	for (size_t i = 0; i < y.length(); i++) {
		y[i] = shifts[i];
	}


	alglib::real_2d_array positions;
	positions.setlength(shifts.size(), 3);
	int partsInFrame = this->PARTITION_AXIS_COUNT * this->PARTITION_AXIS_COUNT;
    int curFrame = -1;
    double cummulativeX, cummulativeY;
    int partSizeX, partSizeY;
	for (size_t i = 0; i < positions.rows(); i++) {
		int frameIndex = i / partsInFrame;
		int partIndex = i % partsInFrame;
        if (partIndex == 0) { //starting next frame
            cummulativeY = 0;
            cummulativeX = 0;
        } else if (partIndex % PARTITION_AXIS_COUNT == 0) { // new partition line
            cummulativeX = 0;
            cummulativeY += partSizeY;
        }

        calculatePartitionSize(partIndex, this->PARTITION_AXIS_COUNT, frameHeight, frameWidth, partSizeX, partSizeY);

		positions[i][0] = cummulativeY + (partSizeY / 2.0);
		positions[i][1] = cummulativeX + (partSizeX / 2.0);
		positions[i][2] = timeStamps[frameIndex];

        cummulativeX += partSizeX;
	}


    alglib::ae_int_t info;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    double epsx = 0.000001;
    double epsf = 0.000001;
    alglib::ae_int_t maxits = 0;
    double diffstep = 0.0001;

    alglib::lsfitcreatef(positions, y, c, diffstep, state);
    alglib::lsfitsetcond(state, epsf, epsx, maxits);
    alglib::lsfitfit(state, calculateShift2);
    alglib::lsfitresults(state, info, c, rep);

    for (size_t i = 0; i < c.length(); i++) {
    	coeffs[i] = c[i];
    }
}

void ProgMovieAlignmentDeformationModel::motionCorrect(const std::vector<MultidimArray<double>>& input,
		std::vector<MultidimArray<double>>& output, const std::vector<double>& timeStamps,
		const std::vector<double>& cx, const std::vector<double>& cy)
{
	for (int i = 0; i < input.size(); i++) {
		applyDeformation(input[i], output[i], cx, cy, timeStamps[i], 0);
	}
}

void ProgMovieAlignmentDeformationModel::applyDeformation(const MultidimArray<double>& input,
		MultidimArray<double>& output, const std::vector<double>& cx, const std::vector<double>& cy,
		double t1, double t2)
{
    int maxX = input.xdim;
    int maxY = input.ydim;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(input)
    {
        int y = i;
        int x = j;

        double posy = y + calculateShift(y, x, t2, cy) - calculateShift(y, x, t1, cy);
        double posx = x + calculateShift(y, x, t2, cx) - calculateShift(y, x, t1, cx);

        int x_left = floor(posx);
        int x_right = x_left + 1;
        int y_down = floor(posy);
        int y_up = y_down + 1;

        double val;
        if (x_right <= 0 || y_up <= 0 || x_left >= maxX -1 || y_down >= maxY - 1) {
            val = 0;
        } else {
            double v11 = dAij(input, y_down, x_left);
            double v12 = dAij(input, y_down, x_right);
            double v21 = dAij(input, y_up, x_left);
            double v22 = dAij(input, y_up, x_right);
            val = linearInterpolation(y_down, x_left, y_up, x_right, v11, v12, v21, v22, posy, posx);
        }

        dAij(output, i, j) = val;
    }	
}

void ProgMovieAlignmentDeformationModel::averageFrames(const std::vector<MultidimArray<double>>& data,
		MultidimArray<double>& out)
{
	out.initZeros(data[0]);
    out.setXmippOrigin();
    for (int i = 0; i < data.size(); i++) {
        out += data[i];
    }
}

void ProgMovieAlignmentDeformationModel::calculatePartitionSize(int partIndex, int edgeCount, int frameHeight,
		int frameWidth, int& partXSize, int& partYSize)
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