*Gravifix*

GraviFix is pre-processing procedure for 3D models represented as triangle meshes in PLY format, either in asci or binary encoding. GraviFix performs on the input model 2 main actions:
1.	Cleans the model from geometrical and topological defects that are typically introduced in the triangle mesh by the acquisition and reconstruction process;
2.	Simplifies the model at three different resolutions in order to make processing afterward efficient.
For each input model (filename.ply), GraviFix produces:
-	A patched model (filename_patched.ply) only if small holes have been filled. The patched model has the same number of vertices of the original model but has more triangles, introduced to fill the holes.
-	A clean model (filename_clean.ply) of nearly the same size of the input model, where geometrical and connectivity errors have been fixed. If no correction is needed, the algorithm saves a clean version anayway (will be the same as the input model)
-	at most 3 simplified models (filename_1M.ply, filename_100K.ply, filename_50K.ply) having 1 million points, 100K points and 50K points, respectively. If the clean model is smaller, simplification might not be needed (at one or all the scales) so not all the simplified files might be produced.
-	One log file (filename.log) reporting all the operation performed on the input model and any error or warning raised by GraviFix
Overall, GraviFix updates a labels.txt file keeping track of statistics on all the processed models. This was needed for research purposes and might be of no interest for you.
Output models are saved in a “processed” subfolder, with respect to the position in the file system of the GraviFix executable. A file naming convention is used to associate input and output models and the clean or simplified versions (see above).
The program is launched through command line on a single model. GraviFix has no graphical user interface nor parameters. To iteratively process a set of 3D models, a script should be used (providing an input text file containing the list of models).  
The application was developed in C++ on Linux CentOS 7 on a 64-bit architecture. It is based on the ImatiSTL library https://sourceforge.net/projects/imatistl/files/ 

The software related to this manual is licensed under the GNU GPL v3.
IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

The GraviFix module and the information contained in GraviFix including, but not limited to, the ideas, concepts and know-how are proprietary, confidential and trade secret to:
	
 IMATI-GE / CNR
 Consiglio Nazionale delle Ricerche
 Istituto di Matematica Applicata e Tecnologie Informatiche “E. Magenes”, Sezione di Genova
 Via de Marini 6, Torre di Francia, 16149, Genova - ITALY

If you use GraviFix for your research paper, please cite:
Michela Mortara, Corrado Pizzi, Michela Spagnuolo: Streamlining the Preparation of Scanned 3D Artifacts to Support Digital Analysis and Processing: the GRAVITATE Case Study. GCH 2017: 165-169

* Build Gravifix
Download and unpack the GraviFix files in a Linux directory. In the same directory, create the folder “ImatiSTL-4.2-2” and download the ImatiSTL library there. Check that the processed directory has been created too.
Download the ImatiSTL library in a subfolder “ImatiSTL-4.2-2”

To build the source code the user needs C++11 or higher. Simply type:
make

This will produce your gravifix executable.
The Makefile expects a 64-bit architecture. You’ll have to edit it if this is not the case.

 
* System Usage
Within GRAVITATE, this component is a pre-processing step performed on the model collections before loading them into the GRAVITATE platform. It could be integrated in the system in case new models or collections are considered. However, due to the complexity of the algorithm, it cannot give results in real-time, but it is better deployed in batch.
Since it is a stand-alone algorithm, GraviFix can easily apply to any model that needs cleaning and simplification; likely it will need to be adapted or customized according to any specific requirements (e.g., number of simplified models, simplification rate, cleaning extent).
The input model must be a triangle mesh, possibly having holes. Small holes will be patched, while bigger holes will be maintained. If more connected components are present, GraviFix will discard the smaller components and work on the biggest component only. 
The expected file format is ply. 
To launch GraviFix simply type the input file name:
 ./gravifix  filename.ply

GraviFix will start analysing the input file, providing statistics (e.g., number of vertices) and a report on the quality of the input mesh. If there are geometrical or topological flaws, GraviFix will proceed to the cleaning operations. 
At the end of the cleaning stage, a filename_clean.ply is created in the  ./processed folder. 
NB: the processed folder must exist before launching GraviFix.
In case small holes were patched in the input model, also a filename_patched.ply will be created.
Then, GraviFix proceeds to simplify the model adaptively, i.e., maintaining the small surface features as much as possible.
If the model is greater than 1 million points, GraviFix will simplify it until it size is less than 1M, and save the simplified model in  ./processed/filename_1M.ply
Simplification continues producing a filename_100K.ply model having less than 100.000 vertices and a filename_50K.ply having less than 50.000 vertices.
If the input model is greater than 1M points, all the three resolutions (filename_1M.ply, filename_100K.ply and filename_50K.ply) will be created in ./processed. If the input model has between 999.999 and 101.000 vertices, the filename_100k.ply and filename_50K.ply will be created, and so on. If the input file is less than 50K vertices, no simplification is performed.

Besides the output models, GraviFix also produces a filename.log text file, reporting details about the performed operations. Furthermore, a cumulative asci file labels.txt is enriched with statistics about all the processed files (do not worry about it, it was used for research purposes). 
