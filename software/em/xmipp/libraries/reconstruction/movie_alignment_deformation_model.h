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

#ifndef _PROG_MOVIE_ALIGNMENT_DEFORMATION_MODEL
#define _PROG_MOVIE_ALIGNMENT_DEFORMATION_MODEL

#include <alglib/src/interpolation.h>
#include <alglib/src/stdafx.h>
#include <float.h>
#include <stdlib.h>
#include "data/xmipp_program.h"
#include "data/metadata_extension.h"
#include "data/xmipp_fftw.h"
#include "data/filters.h"


class ProgMovieAlignmentDeformationModel: public XmippProgram
{
private:
	//----Input arguments----
	FileName fnMovie;		// Input movie stack
	FileName fnMicrograph;	// Output micrograph
    double initDose;        // radiation dose recieved before the first frame
    double perFrameDose;    // radiation dose recieved for each imagined frame 

	int maxIterations;      // max number of iterations for shift calculation
	FileName fnUnaligned;	// Micrograph calculated from unaligned frames
	int upScaling;
    int threadNumbers;
    //correction images
	FileName fnGain;
	FileName fnDark;

	//----Internal data-----
	std::vector<MultidimArray<double> > frames;
	std::vector<double> timeStamps;

	std::vector<double> globalShiftsX;
	std::vector<double> globalShiftsY;

	std::vector<double> localShiftsX;
	std::vector<double> localShiftsY;
	std::vector<std::vector<MultidimArray<double> > > partitions;

    alglib::real_1d_array deformationCoeffsX;
    alglib::real_1d_array deformationCoeffsY;

	const static int PARTITION_COUNT = 5;	//partition count on each axis 
public:
    void readParams();
    void show();
    void defineParams();
    void run();
private:
	void loadMovie(FileName fnMovie,
            std::vector<MultidimArray<double> >& frames,
            std::vector<double>& timeStamps, FileName fnDark, FileName fnGain);

	void estimateShifts(const std::vector<MultidimArray<double> >& data,
            std::vector<double>& shiftsX, std::vector<double>& shiftsY,
            int maxIterations);
	void estimateLocalShifts(
            const std::vector<std::vector<MultidimArray<double> > >& partitions,
			std::vector<double>& shiftsX, std::vector<double>& shiftsY,
            int maxIterations);
	void applyShifts(std::vector<MultidimArray<double> >& data,
            const std::vector<double>& shiftsX,
            const std::vector<double>& shiftsY);

	void partitionFrames(const std::vector<MultidimArray<double> >& frames,
			std::vector<std::vector<MultidimArray<double> > >& partitions,
            int edgeCount);
	void calculatePartitionSize(int partIndex, int edgeCount, int frameHeight,
            int frameWidth, int& partXSize, int& partYSize);

	void calculateModelCoefficients(const std::vector<double>& shifts,
            const std::vector<double>& timeStamps,
            alglib::real_1d_array& coeffs, int frameHeight, int frameWidth);
	static double pixelShift(double x, double y, double t,
            const alglib::real_1d_array& c);
	static void pixelShiftAlg(const alglib::real_1d_array &c,
            const alglib::real_1d_array &dim, double &func, void *ptr);

	void motionCorrect(std::vector<MultidimArray<double> >& data,
            const std::vector<double>& timeStamps,
            const alglib::real_1d_array& cx, const alglib::real_1d_array& cy,
            int scaling);
	void revertDeformation(MultidimArray<double>& input,
            MultidimArray<double>& output, const alglib::real_1d_array& cx,
            const alglib::real_1d_array& cy, double t, int scalingFactor);

	void averageFrames(const std::vector<MultidimArray<double> >& data,
            MultidimArray<double>& out);
	void saveMicrograph(const FileName fnMicrograph,
            const MultidimArray<double>& micrograph);
};

#endif
