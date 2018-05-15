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
	<< "Unaligned micrograph: " << fnUnaligned       << std::endl;
}

void ProgMovieAlignmentDeformationModel::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <metadata>               : Metadata with the list of frames to align");
    addParamsLine("   -o <fn=\"\"> 		          : Give the name of a micrograph to generate an aligned micrograph");
    addParamsLine("  [--maxIterations <N=5>]	  : Number of robust least squares iterations");
    addParamsLine("  [--ounaligned <fn=\"\">]     : Give the name of a micrograph to generate an unaligned (initial) micrograph");
}

void ProgMovieAlignmentDeformationModel::run()
{
	loadMovie(this->fnMovie, this->frames, this->timeStamps);
	this->globalShiftsX.resize(this->frames.size(), 0.0);
	this->globalShiftsY.resize(this->frames.size(), 0.0);

	if (!fnUnaligned.isEmpty()) {
		averageFrames(this->frames, this->unalignedMicrograph);
		saveMicrograph(this->fnUnaligned, this->unalignedMicrograph);
	}

	estimateShifts(this->frames, this->globalShiftsX, this->globalShiftsY, maxIterations);

	applyShifts(this->frames, this->globalShiftsX, this->globalShiftsY);

    std::vector<std::vector<MultidimArray<double>>> partitions;
    partitions.resize(this->PARTITION_AXIS_COUNT * this->PARTITION_AXIS_COUNT, {});
    for (int i = 0; i < partitions.size(); i++) {
        partitions[i].resize(this->frames.size());
    }
    partitionFrames(this->frames, partitions, this->PARTITION_AXIS_COUNT);

	estimateLocalShifts(partitions, this->localShiftsX, this->localShiftsY);

	calculateModelCoefficients(this->timeStamps, this->deformationCoefficientsX, this->deformationCoefficientsY);

	motionCorrect(this->frames, this->correctedFrames, this->timeStamps, this->deformationCoefficientsX,
			this->deformationCoefficientsY);

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

void ProgMovieAlignmentDeformationModel::estimateShifts(const std::vector<MultidimArray<double>>& data,
		std::vector<double>& shiftsX, std::vector<double>& shiftsY, int maxIterations, double minImprovement)
{
	// prepare sum of images
    MultidimArray<double> sum;
    sum.initZeros(data[0]);
    sum.setXmippOrigin();
    for (const MultidimArray<double>& ma : data) {
        sum += ma;
    }

    // prepare shift arrays
    //TODO: mozna radeji error pri nespravne velikosti?
    shiftsY.clear();
    shiftsY.resize(data.size(), 0.0);
    shiftsY.clear();
    shiftsY.resize(data.size(), 0.0);

    // estimate the shifts
    double shiftX;
    double shiftY;
    CorrelationAux aux;
    int cycle = 0;
    double maxDiff = std::numeric_limits<double>::max();
    MultidimArray<double> helper;
    while (cycle < maxIterations && maxDiff > minImprovement) {
        maxDiff = 0;
        cycle++;

        for (int i = 0; i < data.size(); i++) {
            bestShift(sum, data[i], shiftX, shiftY, aux);

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
    }
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

    for (std::vector<MultidimArray<double>>& partStack : partitions) {
        for (int i = 0; i < partStack.size(); i++) {
            int partY = i / edgeCount;
            int partX = i % edgeCount;

            int xSize = partSizeX + (partX < xReminder ? 1 : 0);
            int ySize = partSizeY + (partY < yReminder ? 1 : 0);
            part.resize(1, 1, ySize, xSize);
        }
    }


    int longerPartY = yReminder * (partSizeY + 1);
    int longerPartX = xReminder * (partSizeX + 1);
    for (int fi = 0; fi < frames.size(); fi++) {
        FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(frames[fi]) {
            int partY, partX; // index of the partition
            int innerY, innerX; // index inside the partition
            
            if (i < longerPartY) { // still in the part where partitions are longer along y
                partY = i / (partSizeY + 1)
                innerY = i % (partSizeY + 1)
            } else {
                partY = yReminder + (i - longerPartY) / partSizeY;
                innerY = (i - longerPartY) % partSizeY;
            }

            if (j < longerPartX) { //still in the part where partitions are longer along x
                partX = j / (partSizeX + 1)
                innerX = j % (partSizeX + 1)
            } else {
                partX = xReminder + (j - longerPartX) / partSizeX;
                innerX = (j - longerPartX) % partSizeX;
            }
            dAij(partitions[partY * edgeCount + partX][fi], innerY, innerX) = dAij(frames[fi], i, j);
        }
    }
}

void ProgMovieAlignmentDeformationModel::estimateLocalShifts(
        const std::vector<std::vector<MultidimArray<double>>>& partitions,
        std::vector<double>& shiftsX, std::vector<double>& shiftsY)
{
    double* currX, currY;
    int partitionDepth = partitions[0].size();
    for (int i = 0; i < partitions.size(); i++) {
        currX = shiftsX.front() + partitionDepth * i;
        currY = shiftsY.front() + partitionDepth * i;
        estimateShifts(partitions[i], currX, currY, this->maxIterations);
    }
}

void ProgMovieAlignmentDeformationModel::calculateModelCoefficients(const std::vector<double>& shiftsX,
        const std::vector<double>& shiftsY, const std::vector<double>& timeStamps, std::vector<double>& coeffsX,
        std::vector<double>& coeffsY)
{

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
    size_t maxX = input.colNumber();
    size_t maxY = input.rowNumber();
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(input)
    {
        int y = i;
        int x = j;

        double posy = y + calculateShift(y, x, t2, cy) - calculateShift(y, x, t1, cy);
        double posx = x + calculateShift(y, x, t2, cx) - calculateShift(y, x, t1, cx);
        // std::cout << posx << ":" << x << "   " << posy << ":" << y << std::endl;

        int x_left = floor(posx);
        int x_right = x_left + 1;
        int y_down = floor(posy);
        int y_up = y_down + 1;

        double val;
        if (x_right <= 0 || y_up <= 0 || x_left >= maxX -1 || y_down >= maxY - 1) {
            val = 0;
        } else {
            double v11 = A2D_ELEM(input, y_down, x_left);
            double v12 = A2D_ELEM(input, y_down, x_right);
            double v21 = A2D_ELEM(input, y_up, x_left);
            double v22 = A2D_ELEM(input, y_up, x_right);
            val = linearInterpolation(y_down, x_left, y_up, x_right, v11, v12, v21, v22, posy, posx);
        }

        A2D_ELEM(output, i, j) = val;
    }	
}

void ProgMovieAlignmentDeformationModel::averageFrames(const std::vector<MultidimArray<double>>& data,
		MultidimArray<double>& out)
{
	out.initZeros(data[0]);
    for (int i = 0; i < data.size(); i++) {
        out += data[i];
    }
}