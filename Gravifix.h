#ifndef _GFX_OBJ_H
#define	_GFX_OBJ_H


#include "imatistl.h"

using namespace IMATI_STL;
 
	class GFX_Object : public TriMesh
	{
	public:
	    GFX_Object(): TriMesh() {}
		int GFX_removeSmallestComponents(int *bsize);
		int GFX_strongDegeneracyRemoval(/*int * bsize, */int max_iters=10);
		int GFX_strongIntersectionRemoval(/*int * bsize, */int max_iters=10);
		int GFX_meshclean(/*int minSize, */int max_iters = 10, int inner_loops = 3);
		//int GFX_meshclean2(int max_iters = 10, int inner_loops = 3);
		int GFX_tagNativeBoundaries(int *minSize); 
		//int GFX_fillSmallBoundaries(/*int nbe, */bool refine_patches);
		int GFX_fillSmallBoundaries2(/*int nbe, */bool refine_patches);
		int GFX_simplifyByError(double max_error, int optimal, int edgelen, int check, int SIZE, double *err);
		void GFX_edgeStats(double *avg, double *min, double *max);
		int savePLYwithNormals(const char *fname, bool ascii);
		int checkDegenerate(double minArea);
		double minArea();
	};
 //namespace IMATI_STL
#endif