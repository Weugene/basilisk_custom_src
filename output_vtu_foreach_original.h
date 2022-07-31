#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/**
# output_pvtu_ascii
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_ascii_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_ascii (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}

/**
# output_vtu_asci_foreach
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).

Bug correction: %g turns into scientific notation for high integer values. This is 
not supported by paraview. Hence a fix was needed.
Oystein Lande 2017
*/
void output_vtu_ascii_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(serial, noauto){
    marker[] = _k;
    no_points += 1;
  }
  foreach(serial, noauto){
    no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
    foreach(serial, noauto){
      fprintf (fp, "%g\n", val(s));
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    foreach(serial, noauto){
#if dimension == 2
      fprintf (fp, "%g %g 0.\n", val(v.x), val(v.y));
#endif
#if dimension == 3
      fprintf (fp, "%g %g %g\n", val(v.x), val(v.y), val(v.z));
#endif
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  foreach_vertex(serial, noauto){
#if dimension == 2
    fprintf (fp, "%g %g 0\n", x, y);
#endif
#if dimension == 3
    fprintf (fp, "%g %g %g\n", x, y, z);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(serial, noauto){
#if dimension == 2
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0];
    int ape3 = marker[1,1];
    int ape4 = marker[0,1];
    fprintf (fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4);
#endif
#if dimension == 3
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0,0];
    int ape3 = marker[1,1,0];
    int ape4 = marker[0,1,0];
    int ape5 = marker[0,0,1];
    int ape6 = marker[1,0,1];
    int ape7 = marker[1,1,1];
    int ape8 = marker[0,1,1];
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);

#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "%d \n", i*4);
#endif
#if dimension == 3
    fprintf (fp, "%d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(serial, noauto){
#if dimension == 2
    fputs ("9 \n", fp);
#endif
#if dimension == 3
    fputs ("12 \n", fp);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

/**
# output_pvtu_bin
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_bin_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_bin (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}

/**
# output_vtu_bin_foreach
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on binary format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu_bin() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).

Bug correction: %g turns into scientific notation for high integer values. This is 
not supported by paraview. Hence a fix was needed.
Oystein Lande 2017
*/
void output_vtu_bin_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(serial, noauto){
    marker[] = _k;
    no_points += 1;
  }
  foreach(serial, noauto){
    no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  int count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name,count);
    count += ((no_cells)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name,count);
    count += ((no_cells*3)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n",count);
  count += ((no_points*3)+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(serial, noauto){
#if dimension == 2
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0];
    int ape3 = marker[1,1];
    int ape4 = marker[0,1];
    fprintf (fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4);
#endif
#if dimension == 3
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0,0];
    int ape3 = marker[1,1,0];
    int ape4 = marker[0,1,0];
    int ape5 = marker[0,0,1];
    int ape6 = marker[1,0,1];
    int ape7 = marker[1,1,1];
    int ape8 = marker[0,1,1];
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);

#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "%d \n", i*4);
#endif
#if dimension == 3
    fprintf (fp, "%d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(serial, noauto){
#if dimension == 2
    fputs ("9 \n", fp);
#endif
#if dimension == 3
    fputs ("12 \n", fp);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  unsigned long long block_len=no_cells*8;
#if dimension == 2
  double z=0, vz=0;
#endif
  for (scalar s in list) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach(serial, noauto)
      fwrite (&val(s), sizeof (double), 1, fp);
  }
  block_len=no_cells*8*3;
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach(serial, noauto){
      fwrite (&val(v.x), sizeof (double), 1, fp);
      fwrite (&val(v.y), sizeof (double), 1, fp);
#if dimension == 2
      fwrite (&vz, sizeof (double), 1, fp);
#endif
#if dimension == 3
      fwrite (&val(v.z), sizeof (double), 1, fp);
#endif
    }
  }
  block_len=no_points*8*3;
  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
  foreach_vertex(serial, noauto){
    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

void output_pvd_file(FILE * fp, int nf, float * file_timesteps, char * subname){
//    fputs("<?xml version=\"1.0\"?>\n",fp);
    fputs ("<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
           "\t <Collection>\n", fp);
    for (int i=0; i<=nf; i++) {
        fprintf(fp, "\t\t<DataSet timestep=\"%g\" part=\"0\" file=\"res/%s_0_%4.4d.pvtu\"/>\n", file_timesteps[i], subname, i);
    }//12.9g
    fputs ("\t </Collection>\n "
           "</VTKFile>\n", fp);
}
/* output_vtu_MPI produces *.pvtu files and *.vtu files. The user needs to specify list of scalars and vectors and subname.
*/

struct PVD_output {
    char * subname;
    double myt;
    scalar * list;
    vector * vlist;
};

struct stat st = {0};
static int iter_fp=0;
static float file_timesteps[9999];
void output_vtu_MPI(struct PVD_output o){
    scalar * list = o.list;
    vector * vlist = o.vlist;
    char * subname = o.subname;
    double myt = fabs(o.myt);
    int nf = iter_fp;
    char name_vtu[1000];
    if (iter_fp == 0) {
        if (stat("res", &st) == -1) {
            mkdir("res", 0755);
        }
    }
    FILE *fp;
    if (nf>9999) { fprintf(stderr, "too many files, more than 9999"); exit(1); }
    sprintf(name_vtu, "res/%s_%4.4d_n%3.3d.vtu", subname, nf, pid());
    fp = fopen(name_vtu, "w");
    output_vtu_bin_foreach(list, vlist, 64, fp, true);//64 and true is useless. It needs to support the interface
    fclose(fp);
    if (pid() == 0) {
        //pvtu file
        char name_pvtu[1000], tmp[1000];
	    sprintf(name_pvtu, "res/%s_0_%4.4d.pvtu", subname, nf);
        sprintf(tmp, "%s_%4.4d", subname, nf);
        fprintf(ferr, "vtk: iter_fp = %d myt = %g\n"
                      "<DataSet timestep=\"%g\" part=\"0\" file=\"%s\"/>\n", nf, myt, myt, name_pvtu);
        fp = fopen(name_pvtu, "w");
        output_pvtu_bin(list, vlist, 64, fp, tmp);
        fclose(fp);
        //pvd file with timesteps
        char name_pvd[1000];
        sprintf(name_pvd, "%s.pvd", subname);
        fp = fopen(name_pvd, "w");
        file_timesteps[nf] = myt;
        output_pvd_file(fp, nf, file_timesteps, subname);
        fclose(fp);
    }
    @if _MPI
        MPI_Barrier(MPI_COMM_WORLD);
    @endif
#ifdef DEBUG_OUTPUT_VTU_MPI
    fprintf (ferr, "iter_fp: %d t=%g dt=%g\n", nf, t, dt);
#endif
    iter_fp++;
}


void face_vector2vector(face vector fv, vector mapped_data_lower, vector mapped_data_upper){
//    face vector kappa;
//    kappa = some_face_data;
//    vector mapped_data_lower, mapped_data_upper;
    foreach()
    foreach_dimension(){
        mapped_data_lower.x[] = fv.x[];
        mapped_data_upper.x[] = fv.x[1];
    }
}
