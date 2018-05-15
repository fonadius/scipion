#include <reconstruction/movie_alignment_deformation_model.h>

// RUN_XMIPP_PROGRAM(ProgMovieAlignmentDeformationModel);

int main(int argc, char **argv)
{
    ProgMovieAlignmentDeformationModel padm;
    padm.read(argc, argv);
    padm.run();
}