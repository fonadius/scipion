/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2001 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef _FOURIER_FILTER_HH
   #define _FOURIER_FILTER_HH

#ifdef _HAVE_VTK

#include "../CTF.hh"

/**@name Fourier Masks */
//@{
/** Filter class for Fourier space.

   Example of use for highpass filtering
   \begin{verbatim}
      ImageXmipp I("image.xmp");
      FourierMask Filter;
      Filter.FilterShape=RAISED_COSINE;
      Filter.FilterBand=HIGHPASS;
      Filter.w1=w_cutoff;
      Filter.raised_w=slope;
      Filter.apply_mask(I());
      I.write("filtered_image.xmp");
   \end{verbatim}
*/
class FourierMask {
public:
   #define RAISED_COSINE 1
   #define FROM_FILE     6
   /** Shape of the decay in the filter.
      Valid types are RAISED_COSINE, FROM_FILE. */
   int FilterShape;
   
   #define LOWPASS       1
   #define HIGHPASS      2
   #define BANDPASS      3
   #define STOPBAND      4
   #define CTF           5
   /** Pass band. LOWPASS, HIGHPASS, BANDPASS, STOPBAND, CTF, FROM_FILE */
   int FilterBand;
   
   /** Cut frequency for Low and High pass filters, first freq for bandpass.
       Normalized to 1/2*/
   double w1;
   
   /** Second frequency for bandpass and stopband. Normalized to 1/2 */
   double w2;
   
   /** Pixels around the central frequency for the raised cosine */
   double raised_w;

   /** File from which the mask is read */
   FileName fn_mask;

   /** CTF parameters. */
   XmippCTF ctf;
   
   /** Mask */
   vtkImageData *mask;   
public:
   /** Empty constructor */
   FourierMask() {mask=NULL; clear();}

   /** Destructor */
   ~FourierMask() {clear();}

   /** Assignment */
   FourierMask & operator = (const FourierMask &F);

   /** Clear */
   void clear();

   /** Read parameters from command line.
       If a CTF description file is provided it is read. */
   void read(int argc, char **argv) _THROW;

   /** Show. */
   void show();

   /** Usage. */
   void usage();

   /** Generate mask for a resized image.
   It is supposed that the image is already resized and with its logical
   origin set. If the filter is a CTF it must be already read and prepared
   in the ctf variable.
   
   An exception is thrown if you try to apply a CTF to a volume. */
   void generate_mask(vtkImageData *v) _THROW;

   /** Save mask as a Fourier ImageXmipp or VolumeXmipp. */
   void write_mask(const FileName &fn);

   /** Save amplitude as a ImageXmipp or VolumeXmipp */
   void write_amplitude(const FileName &fn, bool do_not_center=FALSE);

   /** Apply mask.
       It should have been already generated. The given image is modified.
       An exception is thrown if the mask do not fit the size and shape of
       the */
   void apply_mask(vtkImageData *v) _THROW;

   /** Apply mask to image.
       It doesn't need to have a mask already generated. If the mask is equal
       to NULL, the apropiate mask is generated if not, nothing is done with
       the mask except filtering. The given image is modified. */
   void apply_mask(matrix2D<double> &v);

   /** Apply mask to volume.
       The same as the previous one but for volumes*/
   void apply_mask(matrix3D<double> &v);

   /** Resize fourier mask to a desired scale. */
   void resize_mask(int Ydim, int Xdim);
   
   /** Mask power. Return the power of the Fourier Image contained in mask 
   within the given frequencies. */
   double mask_power(double wmin=0, double wmax=0.5);
};
//@}
#endif
#endif
