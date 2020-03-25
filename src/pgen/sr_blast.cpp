//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sr_blast.cpp
//  \brief Problem generator for spherical blast wave in flat spacetime.

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <cstring>    // strcmp()
#include <sys/stat.h>
// Athena++ headers
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"          // ParameterInput

#define MAX_STRING         512

void CartesianRotate( double x[], double theta, double phi, bool inverse ); 
bool Aux_CheckFileExist( const char *FileName );
int  Aux_CountRow( const char *FileName );
template <typename T> int  Aux_LoadTable( T *&Data, const char *FileName, const int NCol_Target, const int TCol[],                       
                                          const bool RowMajor, const bool AllocMem );



static double  Prim_BG[5]                 = { 1.0, 0.0, 0.0, 0.0, 1e-5 };
static double  Blast_Dens_Src             = 1.0;
static double  Blast_Radius_x             = 0.02;
static double  Blast_Radius_y             = 0.02;
static double  Blast_Radius_z             = 0.01;
static int     NumSource                  = 64;
static char    Random_File[MAX_STRING]    = "RandomBlastWaves";
static double *Random_Data                = NULL;
static int     Random_Table_NBin;
static double *Table_x, *Table_y, *Table_z, *Table_theta,  *Table_phi, *Table_Temp;



static Real threshold = 0.01;
 
int RefinementCondition(MeshBlock *pmb);
 
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
  return;
}



void SetTable()
{
   const bool RowMajor_No  = false;                // load data into the column-major order
   const bool AllocMem_Yes = true;                 // allocate memory for Merger_Prof1/2
   const int  NCol         = 6;                    // total number of columns to load
   const int  Col[NCol]    = {1, 2, 3, 4, 5, 6};   // target columns: (radius, density, temperature)
             
   Random_Table_NBin = Aux_LoadTable( Random_Data, Random_File, NCol, Col, RowMajor_No, AllocMem_Yes );
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  AthenaArray<Real> b;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);

  // Initialize hydro variables
  SetTable();

  double Theta, Phi;
  double Cons_BG[5], Prim_EXP[5], Cons_EXP[5];
  double Blast_Center[3], RotatedCartesian[3];

  Table_x     = Random_Data + 0*Random_Table_NBin;
  Table_y     = Random_Data + 1*Random_Table_NBin;
  Table_z     = Random_Data + 2*Random_Table_NBin;
  Table_theta = Random_Data + 3*Random_Table_NBin;
  Table_phi   = Random_Data + 4*Random_Table_NBin;
  Table_Temp  = Random_Data + 5*Random_Table_NBin;

  for (int k=kl; k<=ku; ++k) 
  {
    for (int j=jl; j<=ju; ++j) 
    {
      for (int i=il; i<=iu; ++i)  
      {
        
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        bool InsideEllipsoid[ NumSource ];

        for ( int iBlast = 0; iBlast < NumSource; iBlast++ )
        {
           Blast_Center[0] = Table_x[iBlast];
           Blast_Center[1] = Table_y[iBlast];
           Blast_Center[2] = Table_z[iBlast];
           Theta           = Table_theta[iBlast];
           Phi             = Table_phi[iBlast];


           //1. rotate ellipsoid
           RotatedCartesian[0] = x - Blast_Center[0];
           RotatedCartesian[1] = y - Blast_Center[1];
           RotatedCartesian[2] = z - Blast_Center[2];

           CartesianRotate( RotatedCartesian, Theta, Phi, false );



           InsideEllipsoid[iBlast] = SQR ( RotatedCartesian[0] / Blast_Radius_x ) 
                                   + SQR ( RotatedCartesian[1] / Blast_Radius_y ) 
                                   + SQR ( RotatedCartesian[2] / Blast_Radius_z ) < 1.0;

        }

        bool Inside = false;
        
        for ( int iBlast = 0; iBlast < NumSource; iBlast++ )
        {
           Prim_EXP[4] = Blast_Dens_Src*Table_Temp[iBlast];

           Inside |= InsideEllipsoid[iBlast];

           if ( InsideEllipsoid[iBlast] )
           {
              phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = Prim_EXP[0];
              phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = Prim_EXP[4];
              phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = 0.0;
              phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = 0.0;
              phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = 0.0;
              continue;
           }
        }

        if ( !Inside )
        {
           phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = Prim_BG[0];
           phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = Prim_BG[4];
           phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = 0.0;
           phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = 0.0;
           phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = 0.0;
        }


      }
    }
  }

  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  return;
}

// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &u = pmb->phydro->u;

  Real maxeps = 0.0;

  int level = pmb->loc.level;
  bool Flag=false;
  bool UnFlag=false;

  double Threshold[20];

  Threshold[0] = 0.00001;
  Threshold[1] = 0.00005;
  Threshold[2] = 0.00010;
  Threshold[3] = 0.00050;
  Threshold[4] = 0.00100;

  //double Threshold_Lv = Threshold[level];
  double Threshold_Lv = 0.01;
  double Threshold_Lv_Derefine = 0.01*Threshold_Lv;
  double Diff_i, Diff_j,Diff_k;

    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {

          Diff_i = fabs( u(IEN,k  ,j  ,i+1) - u(IEN,k  ,j  ,i-1) ) * 0.5  / u(IEN,k,j,i);
          Diff_j = fabs( u(IEN,k  ,j+1,i  ) - u(IEN,k  ,j-1,i  ) ) * 0.5  / u(IEN,k,j,i);
          Diff_k = fabs( u(IEN,k+1,j  ,i  ) - u(IEN,k-1,j  ,i  ) ) * 0.5  / u(IEN,k,j,i);
     

          Flag |= Diff_i > Threshold_Lv;
          Flag |= Diff_j > Threshold_Lv;
          Flag |= Diff_k > Threshold_Lv;
  
          UnFlag |=  Diff_i  < Threshold_Lv_Derefine;
          UnFlag |=  Diff_j  < Threshold_Lv_Derefine;
          UnFlag |=  Diff_k  < Threshold_Lv_Derefine;
        }
      }   
    }   

  if ( Flag )   return 1;
  if ( UnFlag ) return -1;

  return 0;
}


#define COMMENT_SYM     '#'      // comment symbol
#define DELIMITER       " \t"    // delimiter characters used by strtok()



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_LoadTable
// Description :  Load the target columns from the table
//
// Note        :  1. Overloaded with different types
//                2. Put the target columns in "TCol[]", which must be sorted into ascending numerical order
//                   in advance
//                3. Allocate memory for the pointer "Data" if AllocMem == true
//                   --> Must be freed manually
//                4. Delimiter characters for strtok() are defined by DELIMITER
//
// Parameter   :  Data        : Pointer to be allocated (if AllocMem == true) and to store the data
//                              --> call-by-reference
//                FileName    : Filename of the target table
//                NCol_Target : Total number of target columns
//                TCol        : Target columns (must be sorted into ascending numerical order in advance)
//                RowMajor    : true/false --> store data into "Data" in the row-/column-major order
//                              -->    Row-major: Data[Row][Column]
//                                  Column-major: Data[Column][Row]
//                AllocMem    : true/false --> allocate/do not allocate memory for the pointer "Data"
//
// Return      :  Total number of matched rows
//-------------------------------------------------------------------------------------------------------
template <typename T>
int Aux_LoadTable( T *&Data, const char *FileName, const int NCol_Target, const int TCol[], const bool RowMajor,
                   const bool AllocMem )
{

// TCol[] must be sorted into ascending numerical order in advance
   for (int t=1; t<NCol_Target; t++)
   {
      if ( TCol[t] <= TCol[t-1] )
         printf( "TCol is not in ascending numerical order ([%d]=%d, [%d]=%d) !!\n",
                    t-1, TCol[t-1], t, TCol[t] );
   }


// count the number of rows
   const int NRow_Target = Aux_CountRow( FileName );


// allocate memory
   if ( AllocMem )   Data = new T [NCol_Target*NRow_Target];


// load data
   char  FirstItem[MAX_STRING];
   int   NCol_Check, NthCol, NCol_Match, NthRow=0, NRow_Match=0;
   char *Token=NULL;

   char *Line = new char [MAX_STRING];
   FILE *File = fopen( FileName, "r" );

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      NCol_Check = sscanf( Line, "%s", FirstItem );
      NthRow ++;

//    skip empty lines or lines starting with the comment symbol
//    --> must check NCheck < 0 as well since EOF is negative
      if ( NCol_Check <= 0  ||  FirstItem[0] == COMMENT_SYM )    continue;

//    loop over all tokens seperated by the delimiter characters
      NthCol     = 0;
      NCol_Match = 0;
      while (  NthCol <= TCol[NCol_Target-1]  &&  ( Token = (NthCol==0)?strtok(Line,DELIMITER):strtok(NULL,DELIMITER) ) != NULL  )
      {
         if ( NthCol == TCol[NCol_Match] )
         {
            if ( RowMajor )   Data[ NRow_Match*NCol_Target + NCol_Match ] = atof( Token );
            else              Data[ NCol_Match*NRow_Target + NRow_Match ] = atof( Token );

            NCol_Match ++;
         }

         NthCol ++;
      }

//    check if we find all target columns
      if ( NCol_Match != NCol_Target )
         printf( "Number of matched columns (%d) != expect (%d) at row %d !!\n",
                    NCol_Match, NCol_Target, NthRow-1 );

      NRow_Match ++;
   } // while ( fgets(Line, MAX_STRING, File) != NULL )

// check if we find all target rows
   if ( NRow_Match != NRow_Target )
      printf( "Number of matched rows (%d) != expect (%d) !!\n", NRow_Match, NRow_Target );

   fclose( File );
   delete [] Line;


// return the total number of matched rows
   return NRow_Match;

} // FUNCTION : Aux_LoadTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CountRow
// Description :  Count the total number of data rows in the target file
//
// Note        :  1. Empty lines and lines starting with the comment symbol will be skipped
//                   --> The comment symbol is defined by COMMENT_SYM
//
// Parameter   :  FileName : Filename of the target table
//
// Return      :  Total number of matched rows
//-------------------------------------------------------------------------------------------------------
int Aux_CountRow( const char *FileName )
{

   if ( !Aux_CheckFileExist(FileName) )
      printf( "table \"%s\" does not exist !!\n", FileName );


   char FirstItem[MAX_STRING];
   int  NRow=0, NItem;

   char *Line = new char [MAX_STRING];
   FILE *File = fopen( FileName, "r" );

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      NItem = sscanf( Line, "%s", FirstItem );

//    skip empty lines and lines starting with the comment symbol
      if ( NItem > 0  &&  FirstItem[0] != COMMENT_SYM )  NRow ++;
   }

   fclose( File );
   delete [] Line;


   return NRow;

} // FUNCTION : Aux_CountRow

bool Aux_CheckFileExist( const char *FileName )
{
 
   struct stat Buf;   
   return ( stat(FileName,&Buf) == 0 );
 
} // FUNCTION : Aux_CheckFileExist



void CartesianRotate( double x[], double theta, double phi, bool inverse )
{
  double xp[3];
 
  if ( inverse )
  {
     xp[0] = -            sin(phi)*x[0] - cos(theta)*cos(phi)*x[1] + sin(theta)*cos(phi)*x[2];
     xp[1] = +            cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*sin(phi)*x[2];
     xp[2] =                            + sin(theta)*         x[1] + cos(theta)*         x[2];
  }
  else
  {
     xp[0] = -            sin(phi)*x[0] +            cos(phi)*x[1];
     xp[1] = - cos(theta)*cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*         x[2];
     xp[2] = + sin(theta)*cos(phi)*x[0] + sin(theta)*sin(phi)*x[1] + cos(theta)*         x[2];
  }
 
  for (int i=0;i<3;i++) x[i] = xp[i];
}                 


// explicit template instantiation
template int Aux_LoadTable <float > ( float  *&, const char *, const int, const int [], const bool, const bool );
template int Aux_LoadTable <double> ( double *&, const char *, const int, const int [], const bool, const bool );
template int Aux_LoadTable <int   > ( int    *&, const char *, const int, const int [], const bool, const bool );
template int Aux_LoadTable <long  > ( long   *&, const char *, const int, const int [], const bool, const bool );

