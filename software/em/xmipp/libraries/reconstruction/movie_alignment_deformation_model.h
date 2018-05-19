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

#include "data/xmipp_program.h"
#include "data/metadata_extension.h"
#include "data/xmipp_fftw.h"
#include "data/filters.h"
#include <alglib/src/interpolation.h>
#include <alglib/src/stdafx.h>


/** Movie alignment correlation Parameters. */
class ProgMovieAlignmentDeformationModel: public XmippProgram
{
private:
	//----Input arguments----
	FileName fnMovie;		// Input movie data
	FileName fnMicrograph;	// Output micrograph
	int maxIterations = 10;	// Limits the number of iterations for shift calculation (Default value taken from Unblur)
	FileName fnUnaligned;	// Micrograph calculated from unaligned frames
	int upScaling = 1;

	//----Internal data-----
	std::vector<MultidimArray<double>> frames;
	std::vector<double> timeStamps;

	std::vector<double> globalShiftsX;
	std::vector<double> globalShiftsY;

	std::vector<double> localShiftsX;
	std::vector<double> localShiftsY;
	std::vector<std::vector<MultidimArray<double>>> partitions;

	std::vector<double> deformationCoefficientsX;
	std::vector<double> deformationCoefficientsY;

	std::vector<MultidimArray<double>> correctedFrames;

	MultidimArray<double> correctedMicrograph;
	MultidimArray<double> unalignedMicrograph;

	const int PARTITION_AXIS_COUNT = 5;	// how many divisions into partitions should be along each frame axis
public:
    /// Read argument from command line
    void readParams();
    /// Show
    void show();
    /// Define parameters
    void defineParams();
    /// Run
    void run();

public:
	void loadMovie(FileName fnMovie, std::vector<MultidimArray<double>>& frames, std::vector<double>& timeStamps);
	void saveMicrograph(FileName fnMicrograph, const MultidimArray<double>& micrograph);
	void saveCalculatedMetadata();

	void estimateShifts(const std::vector<MultidimArray<double>>& data, std::vector<double>& shiftsX,
			std::vector<double>& shiftsY, int maxIterations=50, double minImprovement=0.1);
	/**
	 * @brief [brief description]
	 * @details [long description]
	 * 
	 * @param c [description]
	 * @param dim dim[0] - y, dim[1] - x, dim[2] - t
	 * @param func [description]
	 * @param ptr [description]
	 */
	static void calculateShift2(const alglib::real_1d_array &c, const alglib::real_1d_array &dim, double &func, void *ptr);
	double calculateShift(double x, double y, double t, const std::vector<double>& c);
	void applyShifts(std::vector<MultidimArray<double>>& data, const std::vector<double>& shiftsX,
			const std::vector<double>& shiftsY);

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
	 * TODO: maybe move somewhere else?
	 */
	double linearInterpolation(double y1, double x1, double y2, double x2, double v11, double v12, double v21,
			double v22, double p_y, double p_x);

	void partitionFrames(const std::vector<MultidimArray<double>>& frames,
			std::vector<std::vector<MultidimArray<double>>>& partitions, int edgeCount);

	void estimateLocalShifts(const std::vector<std::vector<MultidimArray<double>>>& partitions,
			std::vector<double>& shiftsX, std::vector<double>& shiftsY);

	void calculateModelCoefficients(const std::vector<double>& shifts, const std::vector<double>& timeStamps,
			std::vector<double>& coeffs, int frameHeight, int frameWidth);

	void motionCorrect(const std::vector<MultidimArray<double>>& input, std::vector<MultidimArray<double>>& output,
			const std::vector<double>& timeStamps, const std::vector<double>& cx, const std::vector<double>& cy);
	void applyDeformation(const MultidimArray<double>& input, MultidimArray<double>& output,
			const std::vector<double>& cx, const std::vector<double>& cy, double t1, double t2);

	void averageFrames(const std::vector<MultidimArray<double>>& data, MultidimArray<double>& out);

	void calculatePartitionSize(int partIndex, int edgeCount, int frameHeight, int frameWidth, int& partXSize,
			int& partYSize);
};

#endif
