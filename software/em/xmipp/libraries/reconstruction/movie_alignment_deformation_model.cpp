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
    fnMaxIterations = getIntParam("--maxIterations");
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
    << "Max iterations:       " << fnMaxIterations   << std::endl
	<< "Unaligned micrograph: " << fnInitialAvg      << std::endl;
}

void ProgMovieAlignmentDeformationModel::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <metadata>               : Metadata with the list of frames to align");
    addParamsLine("   -o <fn=\"\"> 		          : Give the name of a micrograph to generate an aligned micrograph");
    addParamsLine("  [--maxIterations <N=5>]	  : Number of robust least squares iterations");
    addParamsLine("  [--ounaligned]    			  : Give the name of a micrograph to generate an unaligned (initial) micrograph");
}

void ProgMovieAlignmentDeformationModel::run()
{
    // preprocess input data
    MetaData movie;
	readMovie(movie);
	correctLoopIndices(movie);

	Image<double> dark, gain;
	loadDarkCorrection(dark);
	loadGainCorrection(gain);

    int bestIref;
    if (useInputShifts)
    {
    	if (!movie.containsLabel(MDL_SHIFT_X)) { // FIXME seems suspicious
    		setZeroShift(movie);
    	}
    } else {
		bestIref = findShiftsAndStore(movie, dark, gain);
    }

	size_t N, Ninitial;
	Image<double> initialMic, averageMicrograph;
    // Apply shifts and compute average
	applyShiftsComputeAverage(movie, dark, gain, initialMic, Ninitial,
			averageMicrograph, N);

	storeResults(initialMic, Ninitial, averageMicrograph, N, movie, bestIref);
}

