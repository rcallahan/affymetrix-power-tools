#
# sdk/normalization/AFFY-README.txt ---
#
# $Id: AFFY-README.txt,v 1.1 2005-11-02 00:25:33 harley Exp $
#

   A short discussion of the Normalization SDK in which the
hows and whys of this SDK are explained.

----------

   In order to compare the probe intensities across chips,
we need to control for cross chip variation.  This is what
is called 'normalization.  This SDK provides for Quantile
and Median each with an in-core and out-of-core version.

   Most of the functions in the SDK written as templates --
this allows the numeric datatype to be changed for different
applications. (ExACT uses this to provide four datatypes
which can be used.)

   The in-core functions require that all the data be loaded
into memory.  The inputs to the sdk is a vector which
contains a vector of intensity data for each celfile.  This
is specficed as the start and end of the first vector.  The
length of the intensity vector is defined by passing in the
end of the first intensity vector.

Vector of Intensity Vectors
     |
     V
    [0]  =>  ( I[0][0], 1,  2, 3, ..., I[0][M] )
    [1]  =>  ( ... )
    ...
    [N]  =>  ( I[N][0],           ..., I[N][M] )

   In the current scheme, the input data is expected to be
of a uniform type.  If it isnt C++ will give a compile time
error.  This could be hidden by creating a fancy iterator
which "hid" the underling types and did the required
conversions on the fly.  ExACT reformats the celfile
intensity data into a uniform vector before invoking the
normalization code.  The preserves the seperation of the
file SDK and the normalization SDK.

   The in-core median just applies a multiplicative factor
to each intensity value.  This value is supplied by the
caller.

   The in-core quantile normalization makes a mapping and
maps the values based on rank order.  The "vectormap" data
struture was written to reduce the runtime and memory usage
compared to the normal STL "map".  We dont require the data
kept in sorted order as data is added.  Only that it be
sorted before use.  (Only the methods required for the norm
code were written. It is not a complete STL container.)


   Both these in-core methods will use lots of memory. This
makes them unsuitable and unable to normalize many celfiles.
To resolve this, we have written code to do normalizations
out-of-core.  Each celfile is read into memory and a
"Sketch" of the data is taken.  The size of the sketch is
chosen by the user and are normally several thousand points
in size.  The sketch is defined to include the min value as
the first value and the max value as the last value of the
sketch. (This is much smaller than the celfile.)

  The average sketch is computed by averaging the values
across that rank.  Once the average is computed, it is used
with the sketch of the celfile to create a mapping.  Each
celfile is then read into memory again

   This process reduces the peak memory requirement to
((2*size of a celfile)+ ((number of celfiles + 1) * size of
sketch)) at the cost of reading each celfile an extra time
to generate the sketches.


CHANGES SINCE ExACT 1.0 --------------------

   The first release of ExACT did not include the min and
max endpoints when the sketch was generated.  Instead the
Sketch Interpolate function would extrapolate from the end
of the sketch to the required values.  This could lead to
negative values which would cause an assertion error in the
File SDK.  (No negative values allowed.)  By including the
min and max in the sketch the Interpolate function does only
interpolation and negative values are only generated if they
actually exist in the inputs.

   As a result of this change, there will be small
differences in the extreme ends of the outputs when
comparing the output of the last version of ExACT to the
later versions.
