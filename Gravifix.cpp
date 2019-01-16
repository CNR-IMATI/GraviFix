#include <iostream>
#include <fstream>
#include <limits.h>
#include <time.h>
#include <string.h>

#include "Gravifix.h"


#define SMALL_HOLE_EDGES 10			//holes in the input model with boundary smaller than 10 edges will be automatically patched

//change these numbers to achieve output models with different resolutions:
#define PROCESSING_SIZE 1000000 	//target number of vertices in the model simplified for efficient processing
#define MATING_SIZE 100000			//target number of vertices in the model simplified for efficient mating
#define VISUALIZATION_SIZE 50000	//target number of vertices in the model simplified for efficient visualization


#define TVI1(a) (TMESH_TO_INT(((Triangle *)a->data)->v1()->x))
#define TVI2(a) (TMESH_TO_INT(((Triangle *)a->data)->v2()->x))
#define TVI3(a) (TMESH_TO_INT(((Triangle *)a->data)->v3()->x))

#define PIGNOLO 1
#define NORMALS 0  		//switch to 1 to save the output model with normals 
#define CHECK_VERTICES 1

using namespace T_MESH;
using namespace std;
int debug = 0;
bool save_ascii=false;


inline void PRINT_PLY_COMMENT(FILE *f)
{
 if (TMesh::app_name != NULL)
 {
  fprintf(f, "comment File created by %s",TMesh::app_name);
  if (TMesh::app_version != NULL)
  {
   fprintf(f, " v%s",TMesh::app_version);
   if (TMesh::app_year != NULL) fprintf(f, " (%s)",TMesh::app_year);
  }
  fprintf(f, "\n");
  if (TMesh::app_url != NULL) fprintf(f, "comment %s\n",TMesh::app_url);
 }
}




int GFX_Object::savePLYwithNormals(const char *fname, bool ascii)
{
	FILE *fp;
 int i, ii[3];
 float fc[6];  //3 coords, 3 normal coords
 char triname[256];
 unsigned char ii0 = 3;
 Node *n;
 coord *ocds;
 Vertex *v;

 strcpy(triname,fname);
 
 if ((fp = fopen(triname,"w")) == NULL)
 {
  TMesh::warning("Can't open '%s' for output !\n",triname);
  return 1;
 }

 fprintf(fp,"ply\n");
 if (ascii) fprintf(fp,"format ascii 1.0\n");
 else fprintf(fp,"format binary_little_endian 1.0\n");
 PRINT_PLY_COMMENT(fp);
 
 fprintf(fp,"element vertex %d\n",V.numels());
 fprintf(fp,"property float x\n");
 fprintf(fp,"property float y\n");
 fprintf(fp,"property float z\n");
 fprintf(fp,"property float nx\n"); //vertex normal
 fprintf(fp,"property float ny\n");
 fprintf(fp,"property float nz\n");
 
 fprintf(fp,"element face %d\n",T.numels());
 fprintf(fp,"property list uchar int vertex_indices\n");
 fprintf(fp,"end_header\n");
 
 if (ascii) FOREACHVERTEX(v, n) 
	 fprintf(fp, "%f %f %f %f %f %f\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z), TMESH_TO_FLOAT(v->getNormal().x), TMESH_TO_FLOAT(v->getNormal().y), TMESH_TO_FLOAT(v->getNormal().z) );
 else FOREACHVERTEX(v, n)
 {
  fc[0] = TMESH_TO_FLOAT(v->x); fc[1] = TMESH_TO_FLOAT(v->y); fc[2] = TMESH_TO_FLOAT(v->z); 
  fc[3] = TMESH_TO_FLOAT(v->getNormal().x); fc[4] = TMESH_TO_FLOAT(v->getNormal().y); fc[5] = TMESH_TO_FLOAT(v->getNormal().z); //ADD NORMAL!!!
  fwrite(fc, sizeof(float), 6, fp);
 }

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
 i=0; FOREACHVERTEX(v, n) v->x = i++;

 if (ascii) FOREACHNODE(T, n) fprintf(fp,"3 %d %d %d\n",TVI1(n),TVI2(n),TVI3(n));
 else FOREACHNODE(T, n)
 {
  ii[0]=TVI1(n); ii[1]=TVI2(n); ii[2]=TVI3(n);
  fwrite(&ii0, sizeof(unsigned char), 1, fp);
  fwrite(ii, sizeof(int), 3, fp);
 }

 fclose(fp);
 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}


void GFX_Object::GFX_edgeStats(double *avg, double *min, double *max)
{
	Edge *e;
	Node *n;
	double lsum = 0.0, lmin=DBL_MAX, lmax=0.0;
	double l;
	
	FOREACHEDGE(e, n) 
	{
		l= e->length();
		lsum += l;
		if(l<lmin) lmin=l;
		if(l>lmax) lmax=l;
	}
	
	*min=lmin;
	*avg=(lsum / E.numels());
	*max=lmax;
	return;
}



bool GFX_remints_appendCubeToList(Triangle *t0, List& l)
{
	if (!IS_VISITED(t0) || IS_BIT(t0, 6)) return false;

	Triangle *t, *s;
	Vertex *v;
	List triList(t0);
	MARK_BIT(t0, 6);
	coord minx = DBL_MAX, maxx = -DBL_MAX, miny = DBL_MAX, maxy = -DBL_MAX, minz = DBL_MAX, maxz = -DBL_MAX;

	while (triList.numels())
	{
		t = (Triangle *)triList.popHead();
		v = t->v1();
		minx = MIN(minx, v->x); miny = MIN(miny, v->y); minz = MIN(minz, v->z);
		maxx = MAX(maxx, v->x); maxy = MAX(maxy, v->y); maxz = MAX(maxz, v->z);
		v = t->v2();
		minx = MIN(minx, v->x); miny = MIN(miny, v->y); minz = MIN(minz, v->z);
		maxx = MAX(maxx, v->x); maxy = MAX(maxy, v->y); maxz = MAX(maxz, v->z);
		v = t->v3();
		minx = MIN(minx, v->x); miny = MIN(miny, v->y); minz = MIN(minz, v->z);
		maxx = MAX(maxx, v->x); maxy = MAX(maxy, v->y); maxz = MAX(maxz, v->z);
		if ((s = t->t1()) != NULL && !IS_BIT(s, 6) && IS_VISITED(s)) { triList.appendHead(s); MARK_BIT(s, 6); }
		if ((s = t->t2()) != NULL && !IS_BIT(s, 6) && IS_VISITED(s)) { triList.appendHead(s); MARK_BIT(s, 6); }
		if ((s = t->t3()) != NULL && !IS_BIT(s, 6) && IS_VISITED(s)) { triList.appendHead(s); MARK_BIT(s, 6); }
	}

	l.appendTail(new Point(minx, miny, minz));
	l.appendTail(new Point(maxx, maxy, maxz));
	return true;
}

bool GFX_remints_isVertexInCube(Vertex *v, List& loc)
{
	Node *n;
	Point *p1, *p2;
	FOREACHNODE(loc, n)
	{
		p1 = (Point *)n->data; n = n->next(); p2 = (Point *)n->data;
		if (!(v->x < p1->x || v->y < p1->y || v->z < p1->z ||
			v->x > p2->x || v->y > p2->y || v->z > p2->z)) return true;
	}

	return false;
}

void GFX_remints_selectTrianglesInCubes(Basic_TMesh *tin)
{
	Triangle *t;
	Vertex *v;
	Node *n;
	List loc;
	FOREACHVTTRIANGLE((&(tin->T)), t, n) GFX_remints_appendCubeToList(t, loc);
	FOREACHVVVERTEX((&(tin->V)), v, n) if (GFX_remints_isVertexInCube(v, loc)) MARK_BIT(v, 5);
	FOREACHVTTRIANGLE((&(tin->T)), t, n)
	{
		UNMARK_BIT(t, 6);
		if (IS_BIT(t->v1(), 5) || IS_BIT(t->v2(), 5) || IS_BIT(t->v3(), 5)) MARK_VISIT(t);
	}
	FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_BIT(v, 5);
	loc.freeNodes();
}




// basically the same as Basic_TMesh::removeSmallestComponents, but it also returns
// the dimension of the biggest hole created. Necessary for the following fillSmallBoundaries as a parameter
int GFX_Object::GFX_removeSmallestComponents(int *the_max=NULL)
{
	Node *n, *m;
	List todo;
	List components;
	List *component, *biggest = NULL;
	Triangle *t, *t1, *t2, *t3;
	int nt = 0, gnt = 0;

	
	FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);

	t = ((Triangle *)T.head()->data);
	n = T.head();
	do //divides the connected components
	{
		component = new List;
		components.appendHead(component);
		todo.appendHead(t);
		while (todo.numels())
		{
			t = (Triangle *)todo.head()->data;
			todo.removeCell(todo.head());
			if (!IS_BIT(t, 5))
			{
				t1 = t->t1();
				t2 = t->t2();
				t3 = t->t3();

				if (t1 != NULL && !IS_BIT(t1, 5)) todo.appendHead(t1);
				if (t2 != NULL && !IS_BIT(t2, 5)) todo.appendHead(t2);
				if (t3 != NULL && !IS_BIT(t3, 5)) todo.appendHead(t3);

				MARK_BIT(t, 5);
				component->appendTail(t);
			}
		}
		todo.removeNodes();
		for (; n != NULL; n = n->next()) { t = ((Triangle *)n->data); if (!IS_BIT(t, 5)) break; }
	} while (n != NULL);

	int num_comps = components.numels();
     //cout<<"computed "<<num_comps<<" components\n"<<flush;
	 
	FOREACHNODE(components, n)
		if ((nt = ((List *)n->data)->numels()) > gnt) { gnt = nt; biggest = (List *)n->data; } //selects component with highest number of triangles
	
	
	//computing max boundary of the components
	int max_nbc = 0;
	int nbc = 0;
	FOREACHNODE(components, n)
	if ((nt = ((List *)n->data)->numels()) < gnt)  //for smaller components only!
	{
		nbc = 0;
		//std::cout << "new component - " << std::flush;
		component = (List *)n->data;
		//std::cout << "num triangles in this component:" <<component->numels()<<"\n"<<std::flush;
		for (m = component->head(); m != NULL; m = m->next())
		{
			
			// here begin to compute the size of the boundary of each component (in terms of boundary edges)
			t = (Triangle *)m->data;
						
			if(t->t1()==NULL) nbc++; 
			else if (!(IS_BIT(t->t1(), 5))) nbc++;
			
			if (t->t2() == NULL) nbc++;
			else if (!(IS_BIT(t->t2(), 5))) nbc++;
			
			if (t->t3() == NULL) nbc++;
			else if (!(IS_BIT(t->t3(), 5))) nbc++;
			
		}
		//std::cout << "nbc for this component is " << nbc << "\n" << std::flush;
		if (nbc > max_nbc) max_nbc = nbc;
	}
	
	FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5); 

	
	nt = 0;
	FOREACHNODE(components, n)
	{
		
		if (((List *)n->data) != biggest) //removes all components but the biggest
			FOREACHVTTRIANGLE(((List *)n->data), t, m)
		{
			
			if (t->e1->v1 != NULL) t->e1->v1->e0 = NULL;
			if (t->e1->v2 != NULL) t->e1->v2->e0 = NULL;
			if (t->e2->v1 != NULL) t->e2->v1->e0 = NULL;
			if (t->e2->v2 != NULL) t->e2->v2->e0 = NULL;
			if (t->e3->v1 != NULL) t->e3->v1->e0 = NULL;
			if (t->e3->v2 != NULL) t->e3->v2->e0 = NULL;
			t->e1->v1 = t->e1->v2 = t->e2->v1 = t->e2->v2 = t->e3->v1 = t->e3->v2 = NULL;
			t->e1 = t->e2 = t->e3 = NULL;
			nt++;
		}
		
	}
	

	FOREACHNODE(components, n) delete((List *)n->data);

	if (nt)
	{
		d_boundaries = d_handles = d_shells = 1;
		removeUnlinkedElements();
		
		 *the_max= max_nbc;
		
		
		return num_comps - 1;
		
	}

	return 0;
}

int GFX_Object::GFX_strongDegeneracyRemoval(int max_iters)
{
	int n, iter_count = 0;
		bool qstatus = TMesh::quiet;

		TMesh::info("Removing degeneracies...\n");
		while ((++iter_count) <= max_iters && removeDegenerateTriangles()<0)
		{
			std::cout << "there are degeneracies to remove!\n" << std::flush;
			for (n = 1; n<iter_count; n++) growSelection();
			removeSelectedTriangles();
			removeSmallestComponents();
			//GFX_removeSmallestComponents(BSIZE);
			
			//TMesh::quiet = true; fillSmallBoundaries(/*E.numels()*/BSIZE, false); TMesh::quiet = qstatus;
			TMesh::quiet = true; GFX_fillSmallBoundaries2(false); TMesh::quiet = qstatus;
			coordBackApproximation();
		}

		if (iter_count > max_iters) return false;
		return true;
}










void GFX_remints_selectTrianglesInCubes(GFX_Object *tin)
{
	Triangle *t;
	Vertex *v;
	Node *n;
	List loc;
	FOREACHVTTRIANGLE((&(tin->T)), t, n) GFX_remints_appendCubeToList(t, loc);
	FOREACHVVVERTEX((&(tin->V)), v, n) if (GFX_remints_isVertexInCube(v, loc)) MARK_BIT(v, 5);
	FOREACHVTTRIANGLE((&(tin->T)), t, n)
	{
		UNMARK_BIT(t, 6);
		if (IS_BIT(t->v1(), 5) || IS_BIT(t->v2(), 5) || IS_BIT(t->v3(), 5)) MARK_VISIT(t);
	}
	FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_BIT(v, 5);
	loc.freeNodes();
}




int GFX_Object::GFX_strongIntersectionRemoval(int max_iters)
{
	int n, iter_count = 0;
	bool qstatus = TMesh::quiet;

	TMesh::info("Removing self-intersections...\n");

	while ((++iter_count) <= max_iters && selectIntersectingTriangles())
	{
		std::cout << "there are intersections to remove!\n" << std::flush;
		for (n = 1; n<iter_count; n++) growSelection();
		removeSelectedTriangles();
		removeSmallestComponents();
		//GFX_removeSmallestComponents(BSIZE);
		TMesh::quiet = true; GFX_fillSmallBoundaries2(false); TMesh::quiet = qstatus;
		
		coordBackApproximation();
		GFX_remints_selectTrianglesInCubes(this);
	}

	if (iter_count > max_iters) return false;
	return true;
}


//fills all the holes exept for the native ones, tagged with bit 4
int GFX_Object::GFX_meshclean(int max_iters, int inner_loops)
{
bool ni, nd;
Triangle *t;
Node *m;
deselectTriangles();
invertSelection();

for (int n = 0; n<max_iters; n++)
{
	TMesh::info("********* ITERATION %d *********\n", n);
	nd = GFX_strongDegeneracyRemoval(inner_loops);
	
	deselectTriangles(); invertSelection();
	ni = GFX_strongIntersectionRemoval(inner_loops);
	
	if (ni && nd)
	{
		FOREACHTRIANGLE(t, m) if (t->isExactlyDegenerate()) ni = false;
		if (ni) return true;
	}
}

return false;
}





double closestPair(List *bl1, List *bl2, Vertex **closest_on_bl1, Vertex **closest_on_bl2)
{
	Node *n, *m;
	Vertex *v, *w;
	double adist, mindist = DBL_MAX;

	FOREACHVVVERTEX(bl1, v, n)
		FOREACHVVVERTEX(bl2, w, m)
		if ((adist = w->squaredDistance(v))<mindist)
		{
			mindist = adist;
			*closest_on_bl1 = v;
			*closest_on_bl2 = w;
		}

	return mindist;
}

bool joinClosestComponents(Basic_TMesh *tin)
{
	Vertex *v, *w, *gv, *gw=NULL; 
	Triangle *t, *s;
	Node *n;
	List triList, boundary_loops, *one_loop;
	List **bloops_array;
	int i, j, numloops;

	// Mark triangles with connected component's unique ID 
	i = 0;
	FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
	FOREACHVTTRIANGLE((&(tin->T)), t, n) if (t->info == NULL)
	{
		i++;
		triList.appendHead(t);
		t->info = (void *)i;

		while (triList.numels())
		{
			t = (Triangle *)triList.popHead();
			if ((s = t->t1()) != NULL && s->info == NULL) { triList.appendHead(s); s->info = (void *)i; }
			if ((s = t->t2()) != NULL && s->info == NULL) { triList.appendHead(s); s->info = (void *)i; }
			if ((s = t->t3()) != NULL && s->info == NULL) { triList.appendHead(s); s->info = (void *)i; }
		}
	}

	if (i<2)
	{
		FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
		//   JMesh::info("Mesh is a single component. Nothing done.");
		return false;
	}

	FOREACHVTTRIANGLE((&(tin->T)), t, n)
	{
		t->v1()->info = t->v2()->info = t->v3()->info = t->info;
	}

	FOREACHVVVERTEX((&(tin->V)), v, n) if (!IS_VISITED2(v) && v->isOnBoundary())
	{
		w = v;
		one_loop = new List;
		do
		{
			one_loop->appendHead(w); MARK_VISIT2(w);
			w = w->nextOnBoundary();
		} while (w != v);
		boundary_loops.appendHead(one_loop);
	}
	FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_VISIT2(v);

	bloops_array = (List **)boundary_loops.toArray();
	numloops = boundary_loops.numels();

	int numtris = tin->T.numels();
	double adist, mindist = DBL_MAX;

	gv = NULL;
	for (i = 0; i<numloops; i++)
		for (j = 0; j<numloops; j++)
			if (((Vertex *)bloops_array[i]->head()->data)->info != ((Vertex *)bloops_array[j]->head()->data)->info)
			{
				adist = closestPair(bloops_array[i], bloops_array[j], &v, &w);
				if (adist<mindist) { mindist = adist; gv = v; gw = w; }
			}

	if (gv != NULL) tin->joinBoundaryLoops(gv, gw, 1, 0);

	FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
	FOREACHVVVERTEX((&(tin->V)), v, n) v->info = NULL;

	free(bloops_array);
	while ((one_loop = (List *)boundary_loops.popHead()) != NULL) delete one_loop;

	return (gv != NULL);
}
 

char *createFilename(const char *iname,const char *subext, char *oname, const char *newextension, int numver=0)
{
	static char tname[2048];
	static char modelname[2048];
	const char *outdir="./processed/";  //PATH OF THE OUTPUT DIRECTORY - YOU CAN CHANGE THIS
	int n,m;
	char name[2048];
	 
	strcpy(tname, iname); 
	for (n = strlen(tname) - 1; n>0; n--) if (tname[n] == '/') break; 
	int startname=n;
	for(m=n;m<strlen(tname);m++)	if(tname[m]=='.')break;	else modelname[m-startname]=tname[m];
	
	tname[m] = '\0';	
	
	if(numver==0) sprintf(oname, "%s%s%s%s", outdir, modelname, subext, newextension);
	else   if(numver==50000)	sprintf(oname, "%s%s%s_50K%s", outdir, modelname, subext, newextension);
           else if (numver==100000)	sprintf(oname, "%s%s%s_100K%s", outdir, modelname, subext, newextension);
	            else if (numver==1000000)	sprintf(oname, "%s%s%s_1M%s", outdir, modelname, subext, newextension);
		              else sprintf(oname, "%s%s%s_%d%s", outdir, modelname, subext, numver,newextension);
	  
	
	return oname;
}



int GFX_Object::GFX_tagNativeBoundaries(int *minSize)
{
  Node *n;
  Edge *e;
  int ne=0;
  FOREACHEDGE(e, n) if (e->isOnBoundary()) {MARK_BIT(e, 4); ne++;}
  cout<<"taggati "<<ne<<" edges\n"<<flush;
  
  //test count of the number of boundaries
  int nc=0;
  int size=0;
  int minsize=E.numels();
  Edge *nexte;
  Edge *e1, *e2;
  List comp;
  FOREACHEDGE(e, n) UNMARK_VISIT(e);
  
  FOREACHEDGE(e, n) 
  if ((IS_BIT(e,4))&&(!IS_VISITED(e)))
	{
		nc++;
		comp.appendTail(e);
		size=1;
		
		MARK_VISIT(e);
		
		while (comp.numels() > 0)
		{
			nexte = ((Edge *)comp.popHead());
			e1=nexte->v1->nextBoundaryEdge();
			e2=nexte->v2->nextBoundaryEdge();
			if((IS_BIT(e1,4))&&(!IS_VISITED(e1)))
			{
				comp.appendTail(e1);
				size++;
				MARK_VISIT(e1);
			}
			if((IS_BIT(e2,4))&&(!IS_VISITED(e2)))
			{
				comp.appendTail(e2);
				size++;
				MARK_VISIT(e2);
			}
		}
		//end of a component
		if(size<minsize) minsize=size;
	}
	FOREACHEDGE(e, n) UNMARK_VISIT(e);
	
	*minSize=minsize;
	return nc;
	
}








int GFX_Object::GFX_fillSmallBoundaries2(bool refine_patches)
{
  Node *n;
  Edge *e;
  int n_filled=0;
  int newt=0;
  Triangle *t;
  
  FOREACHEDGE(e, n) 
    if(e->isOnBoundary())
	{
		if(!IS_BIT(e,4))  //not native
	    {
		 n_filled++;  
		 
		 newt=TriMesh::TriangulateHole(e);
		 //cout<<"introduced "<<newt<<" new triangles\n"<<flush;
		 if (newt && refine_patches)
		 {			
			t = (Triangle *)T.head()->data;
			refineSelectedHolePatches(t); 
			//NB: this is the version for AMF_Object 			
		 }
		 	    }
	    
	}

  cout<<"filled "<<n_filled<<" holes\n"<<flush;
}

int GFX_Object::checkDegenerate(double ASOGLIA)
{
	int nttt=0;
	Node *n; Triangle *t;
	
	
	FOREACHTRIANGLE(t,n) 
	
	if((t->area())<ASOGLIA)
		  nttt++; 
	
	
	return nttt;
}


double GFX_Object::minArea()
{
 Triangle *t;
 Node *n;
 
 double AMIN=DBL_MAX;
 double averageA=0.0;
 
 int verodeg=0;
 FOREACHTRIANGLE(t,n)  
 {
	 if(t->isExactlyDegenerate()) verodeg++;
	 averageA+=t->area();
	 
	 if((t->area()<AMIN)&&(t->area()>0.0)) 
	  {
		  //printf("AMIN= %.10f\n", t->area());
		  AMIN=t->area();
	  }
 }
 
 averageA/=T.numels();
 printf("average area= %f\n",averageA);
  printf("AMIN= %.16f\n", AMIN);
 
  
 printf("degenerate triangles: %d\n",verodeg);
 
 return AMIN;
}

int main(int argc, char ** argv)
{
	if(argc < 2)
	{
		std::cout <<"Usage: test model.ply"<<endl;
		return 10;
	}

	
	FILE *fp;
	short int IS_CLOSED=-1, G2FIX=-1, C2FIX=-1, SMALL_HOLES_PATCHED=0, COMP_JOINED=0, COMP_REMOVED=0, SIMPL1=0,SIMPL2=0,SIMPL3=0, TAGGED=-1 ; 
	int NVS1,NVS2,NVS3;
	int b_tagged=0;
	double SOGLIA1=0.0,SOGLIA2=0.0,SOGLIA3=0.0,ERR1=0.0,ERR2=0.0, ERR3=0.0;
	double max=0.0,min=0.0,avg=0.0;
	double AMIN, ASOGLIA;
	TMesh::init();
     
	clock_t beginning = clock();

	// Uncomment the following to prevent message reporting
	// TMesh::quiet = true;

	bool stl_output = false;
	bool skip_if_fixed = false;
	bool join_multiple_components = true;
	char infilename[2048], closed_outfilename[2048], clean_outfilename[2048], cleanclosed_outfilename[2048],patched_outfilename[2048], simplifiedP_outfilename[2048], simplifiedM_outfilename[2048],simplifiedV_outfilename[2048], extension[] = ".ply";
	char logfilename[2048];
	char label_start[2048], label_mid[2048], label_mid2[2048],label_end[2048], label_out[2048];
	
	const char *labelfilename="./processed/labels.txt"; //MUST BE CHANGED CONSISTENTLY WITH THE OUTPUT DIRECTORY NAME
	
	int nsh = 0;
	int bsize = 0; 
	int rc = 0;
	char *ms;
	
	sprintf(infilename, "%s", argv[1]);
	//std::cout<<"opening file "<<argv[1]<<"\n"<<std::flush; 
	createFilename(infilename, "", logfilename, ".log");
	fp = fopen(logfilename, "w");
	if (!fp) TMesh::warning("could not open output log file %s \n",logfilename );
	fprintf(fp, "[%d]-Processing started ([clock time] in ms)\n", clock());
	
	
	//LOAD
	
    GFX_Object *obj = new GFX_Object();
	if(obj->loadPLY(argv[1]) == -1) 
	{ 
      fprintf(fp, "[%d]- could not load input file %s- exit \n", clock(), argv[1]);
	  fclose(fp);
      TMesh::error("could not load input file %s- exit", argv[1]); 
	  exit(1); 
    }
	std::cout<<"from input  file "<<argv[1]<<"\n"<<std::flush; 
	TMesh::info("Num. of shells: %d  Num. of boundaries: %d\n", obj->shells(), obj->boundaries());
	fprintf(fp, "[%d]-File loaded\n*** object has:  \n\t%d shells \n\t%d vertices \n\t%d triangles \n\t%d (native) boundaries\n***\n", clock(), obj->shells(),obj->V.numels(), obj->T.numels(),obj->boundaries());
	
	
	
	//CHECK NORMALS
	
	printf("\n\nCHECKING NORMALS\n\n");
	
	 Triangle *t0=(Triangle *)obj-> T.head()->data;
     Triangle *top = obj->topTriangle(t0);
	
     if (top->getNormal().z < 0) {TMesh::warning("Wrong normals, flipping..."); obj->flipNormals(t0);}
	 else TMesh::info("Normals are correctly oriented");
   
		
	sprintf(label_start, "*%d_%d_%d_%d*", obj->V.numels(), obj->T.numels(), obj->shells(), obj->boundaries());	
	if(obj->boundaries()==0) IS_CLOSED=1; else IS_CLOSED=0; 
	fclose(fp);
	fp = fopen(logfilename, "a");
    

	//REDUCTION TO SINGLE SHELL
	
	if (obj->shells() > 1)
	{
			fprintf(fp, "[%d]- more than one shell, trying reduction \n", clock());
			if (join_multiple_components)
			{
				TMesh::info("Joining input components ...\n");
				fprintf(fp, "[%d]-Joining multiple components...\n", clock());
				TMesh::begin_progress();

				while (joinClosestComponents(obj)) TMesh::report_progress("Num. components: %d       ", obj->shells());
				TMesh::end_progress();
				fprintf(fp, "[%d]-Done (%d components left)\n", clock(), obj->shells());
				obj->deselectTriangles();
				COMP_JOINED=1;
			}

			// Keep only the largest component (i.e. with most triangles)
			std::cout << "Removing smallest components...\n" << std::flush;
			fprintf(fp, "[%d]-Removing smallest components...\n", clock());

			rc = obj->GFX_removeSmallestComponents(&bsize); //bsize is initialized inside!
			 
			
			if (rc) {TMesh::warning("Removed %d small component(s); BSIZE= %d\n", rc,bsize); COMP_REMOVED=1;}
			else  TMesh::info("No component removed\n");
			fprintf(fp, "[%d]-Done (%d components removed, %d components left)\n", clock(), rc, obj->shells());
			
	}
    else fprintf(fp, "[%d]-single shell, no need to join or remove components\n", clock());
		
	fclose(fp);
	fp = fopen(logfilename, "a");

    
	
	//CHECK GEOMETRY AND TOPOLOGY
	
			
	TMesh::info("\n***Starting Connectivity and Geometry Check***\n\n");
	fprintf(fp, "[%d]-Starting Connectivity and Geometry check...\n", clock());
	//these just check, do not modify
	ms=obj->checkConnectivity();
	if(ms) printf("***checkConnectivity: %s\n",ms);
	if((ms=obj->checkConnectivity())!=NULL) 
    {
		C2FIX=1;
		//TMesh::info("Connectivity has problems\n");
	    fprintf(fp, "[%d]-there are connectivity issues: %s\n", clock(),ms);
	}
    else
	{
		C2FIX=0;
        TMesh::info("CheckConnectivity passed\n");
	    fprintf(fp, "[%d]-checkConnectivity passed\n", clock());
	}
	
	if(obj->checkGeometry()!=NULL)       //returns the first problematic element
	{
		G2FIX=1;
		//TMesh::info("Geometry has problems\n");
	    fprintf(fp, "[%d]-there are geometry issues\n", clock());
	}
    else
	{
		G2FIX=0;
        TMesh::info("CheckGeometry passed\n");
	    fprintf(fp, "[%d]-checkGeometry passed\n", clock());
	}
	
	
    sprintf(label_mid, "%d_%d_%d_%d_%d_%d_%d_%d*", COMP_JOINED, COMP_REMOVED, obj->V.numels(), obj->T.numels(),obj->shells(), obj->boundaries(),G2FIX, C2FIX  );
	fclose(fp);
	fp = fopen(logfilename, "a");	
	
	
	
	//PATCHING SMALL HOLES
	
	if (obj->boundaries())  
		{
			TMesh::info("Patching small holes\n");
			fprintf(fp, "[%d]-There are boundaries\n", clock());
			fprintf(fp, "[%d]-bsize= %d, SMALL_HOLE_EDGES=%d\n", clock(),bsize,SMALL_HOLE_EDGES);
			fprintf(fp, "[%d]-Patching small holes with threshold %d boundary edges\n", clock(),(bsize>SMALL_HOLE_EDGES)?bsize:SMALL_HOLE_EDGES);
		    nsh=obj->fillSmallBoundaries(((bsize>SMALL_HOLE_EDGES)?bsize:SMALL_HOLE_EDGES), false);
			fprintf(fp, "[%d]-Done, patched %d small holes, %d boundaries remaining\n", clock(), nsh, obj->boundaries());
			
			if (nsh > 0)
			{
				SMALL_HOLES_PATCHED=1;
				createFilename(infilename, "_patched", patched_outfilename, extension);
				TMesh::info("Saving partial output mesh in %s ...\n",patched_outfilename);
				fprintf(fp, "[%d]-Saving patched mesh in %s...\n", clock(),patched_outfilename);
				obj->savePLY(patched_outfilename, save_ascii);
				fprintf(fp, "[%d]-Done\n", clock());
			}
		
		}
	else fprintf(fp, "[%d]-no boundaries to patch\n", clock());
	int nb=obj->boundaries();
	if (nb) TMesh::info("There are still holes (%d boundary loops)\n", nb);
	 		
	
	
	
	//TAGGING NATIVE BOUNDARIES
	
	int minSize=0;
	if (nb) 
	{
		TMesh::info("Tagging native boundaries \n");
		b_tagged=obj->GFX_tagNativeBoundaries(&minSize);	
		cout<<"minSize="<<minSize<<"\n"<<flush;
		if(!b_tagged) TMesh::warning("no native boundary out of %d was actually tagged!");
		if(b_tagged<nb) TMesh::warning("%d native boundary could not be tagged! remaining: %d", nb-b_tagged);
		TAGGED=1; 
	    cout<<"Tagged "<<b_tagged<<" native boundaries\n"<<flush;
	    fprintf(fp, "[%d]-Tagged %d native boundaries, minSize=\n", clock(),b_tagged, minSize);
		
	}
	else {TAGGED=0;fprintf(fp, "[%d]-no native boundaries to tag\n", clock());}
	
	sprintf(label_mid2, "%d_%d_%d_%d_%d*", SMALL_HOLES_PATCHED, nsh, nb, TAGGED,b_tagged);
	fclose(fp);
	fp = fopen(logfilename, "a");
	
	bsize=(minSize-SMALL_HOLE_EDGES>0)?(minSize-SMALL_HOLE_EDGES):(SMALL_HOLE_EDGES);
	
	
	
	// RUN GEOMETRY CORRECTION
	
	if(G2FIX||C2FIX||PIGNOLO)
	{
		TMesh::info("Fixing degeneracies and/or intersections...\n");
		fprintf(fp, "[%d]-Fixing degeneracies and intersections, bsize=%d\n", clock(),bsize);
		if (!(obj->GFX_meshclean()))
		//if (!(obj->GFX_meshclean2()))	
		{
		  TMesh::warning("GraviFix could not fix everything.\n");
		  fprintf(fp, "[%d]-GraviFix could not fix everything. Skipping.\n", clock());	
		}
		
		else 
		{
			fprintf(fp, "[%d]-Done\n", clock());
			
			obj->GFX_strongIntersectionRemoval();
			TMesh::info("Saving clean output mesh ...\n");
			createFilename(infilename,"_clean", clean_outfilename, extension);
			obj->savePLY(clean_outfilename, save_ascii);
			
            if (obj->boundaries()) 
			{
				fprintf(fp, "[%d]-Saving clean mesh in %s, still %d holes\n", clock(), clean_outfilename, obj->boundaries());
				TMesh::info("%d holes still to be filled\n", obj->boundaries());
			}
			else { 
			       fprintf(fp, "[%d]-Saving clean (and closed) mesh in %s...\n", clock(), cleanclosed_outfilename);
			       TMesh::info("no holes to be filled\n");
			     }
			fprintf(fp, "[%d]-Done\n", clock());
			TMesh::info("Elapsed time: %d ms\n", clock() - beginning);
		}
	}
	else fprintf(fp, "[%d]-no degenaracies or intersections to fix\n", clock());

	
	//obj->GFX_strongIntersectionRemoval(); QUESTO non credo che serva, appena fatto in meshclean... però controllare il _clean che non abbia intersezioni!!!!
	
	TMesh::info("Saving anyway clean output mesh ...\n");
	createFilename(infilename,"_clean", clean_outfilename, extension);
	if(NORMALS) obj->savePLYwithNormals(clean_outfilename, save_ascii);
	else  obj->savePLY(clean_outfilename, save_ascii);
	
    
	
	fclose(fp);
	fp = fopen(logfilename, "a");
 
 
 
    ////////////////////////////////////
	//CHECK the MIN AREA - ADDED BY ICCHIA
    AMIN=obj->minArea();
    //printf("dopo open, MIN AREA OF THE ORIGINAL MESH= %.16f\n",AMIN);
	ASOGLIA=AMIN;
	
	///////////////////////////////////////
 
 
	//OPTIONAL: Closes all the holes
	/*	
	
	//CLOSE THE MODEL
	
	if (obj->boundaries())
	{
			TMesh::info("Patching all the remaining holes\n");
			//obj->fillSmallBoundaries(0,true);
			fprintf(fp, "[%d]-Patching all the remaining %d holes...\n", clock(), obj->boundaries());
			obj->GFX_fillSmallBoundaries(0, true);
			
            if (obj->boundaries()) 
			 { TMesh::warning("Gravifix could not close all the holes\n"); 
			  fprintf(fp, "[%d]-WARNING: GraviFix could not close all the holes. %d boundaries remaining\n", clock(), obj->boundaries()); 
			 }
			else
			{ 
		        createFilename(infilename, "_closed", closed_outfilename, extension);
				fprintf(fp, "[%d]-Done\n", clock());
				TMesh::info("Saving closed output mesh ...\n");
				fprintf(fp, "[%d]-Saving closed mesh in %s...\n", clock(), closed_outfilename);
				obj->savePLY(closed_outfilename, save_ascii);
				fprintf(fp, "[%d]-Done\n", clock());
				TMesh::info("Elapsed time: %d ms\n", clock() - beginning);
			}
	}
    */
	
	
	//SIMPLIFICATION
	
	//printf("ASOGLIA=%.16f\n",ASOGLIA);	
	//printf("triangloli iniziali con area minore di ASOGLIA= %d\n",obj->checkDegenerate(ASOGLIA));
	
	int nvprima=obj->V.numels();
	NVS1=NVS2=NVS3=obj->V.numels();
	obj->GFX_edgeStats(&avg,&min,&max);
	double unit=avg;
	int nforced=-1;		//! forceNormalConsistence returns: 0 if mesh was already oriented; 1 if
						//! the mesh could be oriented without cuts; >1 if cuts were necessary.
		
	fprintf(fp, "before simplification: avg_edge= %f, min_edge= %f, max_edge=%f\n", avg,min,max);
	
	if(obj->V.numels()>PROCESSING_SIZE)
	{
		TMesh::info("Starting Simplification for PROCESSING...\n");
		
		fprintf(fp, "[%d]-First simplification round, original vertex count= %d\n", clock(), nvprima);
		SOGLIA1=0.05*unit;
		if (obj->simplify(PROCESSING_SIZE, 1, 0, 2, ASOGLIA)) 
		//if (obj->GFX_simplifyByError(SOGLIA1, 1, 0, 2, PROCESSING_SIZE, &ERR1)) //using quadratic error function, optimal point, best check_collapse
		{
			fprintf(fp, "[%d]-Done, ", clock());
			
			NVS1 = obj->V.numels();
			SIMPL1=1;
			
			std::cout << "Number of Vertices: " << NVS1 << " " << std::flush;
			std::cout << "Simplification rate= " << 100 - (100 * (float)NVS1 / (float)nvprima) << "% \n" << std::flush;
			fprintf(fp, "Simplification rate = %f percent, Number of Vertices = %d \n", 100 - (100 * (float)NVS1 / (float)nvprima),NVS1);
			TMesh::info("Saving simplified output mesh for PROCESSING...\n");
			
			/**/
			if (PIGNOLO) obj->GFX_strongIntersectionRemoval();
			
										  
			createFilename(infilename,"", simplifiedP_outfilename, extension, NVS1);
			fprintf(fp, "[%d]-Saving simplified mesh in %s...\n", clock(), simplifiedP_outfilename);
			if(NORMALS) obj->savePLYwithNormals(simplifiedP_outfilename, save_ascii);
			else obj->savePLY(simplifiedP_outfilename, save_ascii);
			fprintf(fp, "[%d]-Done\n", clock());
			TMesh::info("Elapsed time: %d ms\n", clock() - beginning);
		}
		
		
		ms=obj->checkConnectivity();
	    if(ms) printf("***checkConnectivity: %s\n",ms);
		 
		obj->GFX_edgeStats(&avg,&min,&max);
		fprintf(fp, "after first round: avg_edge= %f, min_edge= %f, max_edge=%f, error reached= %f\n", avg, min,max, ERR1);
	}
	else fprintf(fp, "[%d]-no need to simplify for processing: %d vertices\n", clock(),obj->V.numels());
	
	if(obj->V.numels()>MATING_SIZE)	
	{
		TMesh::info("Further Simplification ...\n");
		fprintf(fp, "[%d]-Second simplification round, vertex count= %d\n", clock(), obj->V.numels());
		SOGLIA2=0.09*unit;
		if (obj->simplify(MATING_SIZE, 1, 0, 2,ASOGLIA)) //using quadratic error function, no optimal point, best check_collapse
		//if (obj->GFX_simplifyByError(SOGLIA2, 1, 0, 2, MATING_SIZE,&ERR2)) //using quadratic error function, optimal point, best check_collapse
		{ 
			fprintf(fp, "[%d]-Done, ", clock());
			NVS2 = obj->V.numels();
			SIMPL2=1;
			std::cout << "Number of Vertices: " << NVS2 << " " << std::flush;
			std::cout << "Simplification rate= " << 100 - (100 * (float)NVS2 / (float)nvprima) << "% \n" << std::flush;
			fprintf(fp, "Simplification rate = %f percent, Number of Vertices = %d \n", 100 - (100 * (float)NVS2 / (float)nvprima), NVS2);
			TMesh::info("Saving simplified output mesh for MATING...\n");
			
            /**/if(PIGNOLO) obj->GFX_strongIntersectionRemoval();
           
		  			
			createFilename(infilename, "", simplifiedM_outfilename, extension, NVS2);
			fprintf(fp, "[%d]-Saving simplified mesh in %s...\n", clock(), simplifiedM_outfilename);
			if(NORMALS) obj->savePLYwithNormals(simplifiedM_outfilename, save_ascii);
			else obj->savePLY(simplifiedM_outfilename, save_ascii);
			fprintf(fp, "[%d]-Done\n", clock());
			TMesh::info("Elapsed time: %d ms\n", clock() - beginning);
		}
		
		
		ms=obj->checkConnectivity();
	    if(ms) printf("***checkConnectivity: %s\n",ms);
		
		
		obj->GFX_edgeStats(&avg,&min,&max);
		fprintf(fp, "after second round: avg_edge= %f, min_edge= %f, max_edge=%f, error reached=%f\n", avg, min,max,ERR2);
	}
	else fprintf(fp, "[%d]-no need to simplify for mating: %d vertices\n", clock(),obj->V.numels());	
	
		
	if(obj->V.numels()>VISUALIZATION_SIZE)
	{
		TMesh::info("Further Simplification ...\n");
		fprintf(fp, "[%d]-Third simplification round, vertex count= %d\n", clock(), obj->V.numels());
		SOGLIA3=0.5*unit;
		if (obj->simplify(VISUALIZATION_SIZE, 1, 0, 2,ASOGLIA)) 
		//if (obj->GFX_simplifyByError(SOGLIA3, 1, 0, 2, VISUALIZATION_SIZE,&ERR3)) //using quadratic error function, optimal point, best check_collapse
		{ 
			fprintf(fp, "[%d]-Done, ", clock());
			NVS3 = obj->V.numels();
			SIMPL3=1;
			std::cout << "Number of Vertices: " << NVS3 << " " << std::flush;
			std::cout << "Simplification rate= " << 100 - (100 * (float)NVS3 / (float)nvprima) << "% \n" << std::flush;
			fprintf(fp, "Simplification rate = %f percent, Number of Vertices = %d \n", 100 - (100 * (float)NVS3 / (float)nvprima), NVS3);
			TMesh::info("Saving simplified output mesh for VISUALIZATION...\n");
			
            /**/if (PIGNOLO) obj->GFX_strongIntersectionRemoval();	
			
									  
			createFilename(infilename,"",simplifiedV_outfilename, extension, NVS3);
			fprintf(fp, "[%d]-Saving simplified mesh in %s...\n", clock(), simplifiedV_outfilename);
			if(NORMALS) obj->savePLYwithNormals(simplifiedV_outfilename, save_ascii);
			else obj->savePLY(simplifiedV_outfilename, save_ascii);
			fprintf(fp, "[%d]-Done\n", clock());
			TMesh::info("Elapsed time: %d ms\n", clock() - beginning);
		}
		
		
		obj->GFX_edgeStats(&avg,&min,&max);
		fprintf(fp, "after third round: avg_edge= %f, min_edge= %f, max_edge=%f, error reached=\n", avg, min,max, ERR3);
	}
	else fprintf(fp, "[%d]-no need to simplify for visualization: %d vertices\n", clock(),obj->V.numels());	
	
	
	fclose(fp);
	fp = fopen(logfilename, "a");
	
	
	
	//FINAL CHECK	
	
	TMesh::info("\n***Starting Final Connectivity and Geometry Check***\n\n");
	fprintf(fp, "[%d]-Starting Final Connectivity and Geometry check...\n", clock());
	//these just check, do not modify
	ms=obj->checkConnectivity();
	if(ms)printf("%s\n",ms);
	if((obj->checkConnectivity())!=NULL) 
    {
		C2FIX=1;
		TMesh::info("Connectivity has problems\n");
	    fprintf(fp, "[%d]-there are connectivity issues\n", clock());
	}
    else
	{
		C2FIX=0;
        TMesh::info("CheckConnectivity passed\n");
	    fprintf(fp, "[%d]-checkConnectivity passed\n", clock());
	}
	obj->checkGeometry();
	if(obj->checkGeometry()!=NULL)       //returns the first problematic element
	{
		G2FIX=1;
		TMesh::info("Geometry has problems\n");
	    fprintf(fp, "[%d]-there are geometry issues\n", clock());
	}
    else
	{
		G2FIX=0; 
        TMesh::info("CheckGeometry passed\n");
	    fprintf(fp, "[%d]-checkGeometry passed\n", clock());
	}
	if (PIGNOLO) obj->GFX_strongIntersectionRemoval();		
	
	
	fclose(fp);
	
	fp = fopen(logfilename, "a");
	
	
	fprintf(fp, "[%d]-Processing Finished. Goodbye!\n", clock());	
	fclose(fp);	
   
	sprintf(label_end, "%d_%d_%d_%d_%d_%d*%d_%d_%d_%d_%d_%d*", SIMPL1,NVS1, SIMPL2,NVS2,SIMPL3,NVS3,obj->V.numels(), obj->T.numels(), obj->shells(), obj->boundaries(),G2FIX, C2FIX);
	
	
	
	//MERGE of all labels and append the new string to the labels.txt file 
	
	sprintf(label_out,"%s%s%s%s",label_start,label_mid,label_mid2,label_end);
	
	fp = fopen(labelfilename, "a");
	fprintf(fp,"%s %s\n",argv[1],label_out);
	fclose(fp);
	
	printf("\nCHECK DEGENERACIES on the final mesh\n\n");
	printf("triangles with area less than ASOGLIA= %d\n",obj->checkDegenerate(ASOGLIA));
	
	printf("Mission accomplished! Goodbye.\n");
	delete obj;
	
		
	
    return 0;  
}
